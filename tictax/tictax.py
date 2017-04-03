import os
import sys
import gzip
import json
import tqdm
import asyncio
import aiohttp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from tictax import tictax


def config():
    conf = {}
    conf_path = os.path.join(os.path.expanduser('~'), '.tictax')
    if not os.path.exists(conf_path):
        print('----------------\n- TICTAX SETUP -\n----------------')
        print('Tictax needs a One Codex API key to do streaming sequence classification')
        print('1) Sign up for a One Codex account at https://www.onecodex.com')
        print('2) Paste your API key below and it will be saved for future sessions')
        conf['one_codex_api_key'] = input('API key: ').strip()
        with open(conf_path, 'w') as conf_fh:
            json.dump(conf, conf_fh)
        print('Config saved to {}'.format(conf_path))
        return conf
    else:
        with open(conf_path, 'r') as conf_fh:
            return json.load(conf_fh)


def parse_seqs(seqs_path, fastq):
    is_fastq = True if seqs_path.endswith(('.fastq','.fastq.gz','.fq','.fq.gz')) or fastq else False
    if seqs_path.endswith('.gz'):
        with gzip.open(seqs_path, 'rt') as gzip_fh:
            if is_fastq:  # Compressed fastq
                return list(SeqIO.parse(gzip_fh, 'fastq'))
            else:  # Compressed fasta
                return list(SeqIO.parse(gzip_fh, 'fasta'))
    elif is_fastq:  # Uncompressed fastq
        return list(SeqIO.parse(seqs_path, 'fastq')) 
    else:  # Assume uncompressed fasta
        return list(SeqIO.parse(seqs_path, 'fasta'))


# --------------------------------------------------------------------------------------------------


def build_record(id, classification):
    description = (classification['classifier']
                   + '|' + str(classification['taxid'])
                   + '|' + classification['sciname']
                   + '|' + classification['rank']
                   + '|' + ';'.join(classification['lineage']))
    record = SeqRecord(Seq(classification['sequence'], IUPAC.ambiguous_dna),
                       id=id, description=description)
    return record


async def oc_classify_single(session, sequence_id, sequence):
    url = 'https://app.onecodex.com/api/v0/search'
    payload = {'sequence': str(sequence)[:25000]}
    try:
        async with session.post(url, data=payload, timeout=600) as response:
            r = await response.json()
    except (asyncio.TimeoutError, aiohttp.errors.ClientOSError, json.decoder.JSONDecodeError):
        r = {}
    return {'k': r.get('k', ''),
            'sequence': sequence,
            'classifier': 'one_codex_rt',
            'taxid': r.get('tax_id', 0),
            'support': round(r.get('n_hits', 0)/max(r.get('n_lookups', 1), 1), 4)}


async def taxify_ebi(session, sequence_id, taxid):
    url = 'http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{}'
    template = {'sciname': '', 'lineage': '', 'rank': ''}
    if taxid == 0:
        return template
    elif taxid == 1:
        template['sciname'] = 'root'
        return template
    try:
        async with session.get(url.format(taxid), timeout=600) as response:
            r = await response.json()
            template['sciname'] = r.get('scientificName', '')
            template['lineage'] = r['lineage'].strip(' ;').split('; ') if 'lineage' in r else []
            template['rank'] = r.get('rank', '')
            return template
    except (asyncio.TimeoutError, aiohttp.errors.ClientOSError, json.decoder.JSONDecodeError):
        return template


async def classify_taxify(oc_session, ebi_session, sequence_id, sequence):
    classification = await oc_classify_single(oc_session, sequence_id, sequence)
    taxid = classification.get('taxid')
    taxification = await taxify_ebi(ebi_session, sequence_id, taxid)
    return sequence_id, {**classification, **taxification}  # Merge dicts


async def oc_classify(records, one_codex_api_key, progress=False, stdout=False):
    oc_auth = aiohttp.BasicAuth(one_codex_api_key)
    conn = aiohttp.TCPConnector(limit=10)
    with aiohttp.ClientSession(auth=oc_auth, connector=conn) as oc_session:
        with aiohttp.ClientSession(connector=conn) as ebi_session:
            tasks = [classify_taxify(oc_session, ebi_session, r.id, str(r.seq)) for r in records]
            # No async generators in 3.5... :'(
            # return [await f for f in tqdm.tqdm(asyncio.as_completed(tasks), total=len(tasks))]
            records = []
            for f in tqdm.tqdm(asyncio.as_completed(tasks),
                               disable=not progress,
                               total=len(tasks)):
                response = await f
                record = build_record(response[0], response[1])
                if stdout:
                    print(record.format('fasta'), end='')
                records.append(record)
            return records


# --------------------------------------------------------------------------------------------------


def kmer_lca_records(seqs_path,
                     one_codex_api_key: 'One Codex API key' = None,
                     fastq: 'input is fastq; disable autodetection' = False,
                     progress: 'show progress bar (sent to stderr)' = False):
    '''
    Parallel lowest common ancestor sequence classification of fasta/q using the One Codex API.
    Returns Biopython SeqRecords with tictax annotations as the `description` attribute.
    LCAs are assigned using an LCA index of 31mers from the One Codex database.
    '''
    records = parse_seqs(seqs_path, fastq)
    one_codex_api_key = one_codex_api_key if one_codex_api_key else config()['one_codex_api_key']
    print('Classifying sequences…', file=sys.stderr)
    records = asyncio.get_event_loop().run_until_complete(oc_classify(records,
                                                                      one_codex_api_key,
                                                                      progress,
                                                                      False))
    return records


if __name__ == '__main__':
    tictax.cli.main()
