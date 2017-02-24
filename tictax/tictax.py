import os
import sys
import json
import tqdm
import asyncio
import aiohttp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

import tictax


def config():
    conf = {}
    conf_path = os.path.join(os.path.expanduser('~'), '.tictax')
    if not os.path.exists(conf_path):
        print('----------------\n- TICTAX SETUP -\n----------------')
        print('Tictax needs a One Codex API key to do streaming sequence classification')
        print('1) Sign up for a One Codex account at https://www.onecodex.com')
        print('2) Paste your API key below and it will be saved for future sessions')
        conf['one_codex_api_key'] = input('API key: ').strip()
        with open(conf_path, 'w') as conf_fh: # Needs storing in a sensible place... dotfile?
            json.dump(conf, conf_fh)
        print('Config saved to {}'.format(conf_path))
        return conf
    else:
        with open(conf_path, 'r') as conf_fh:
            return json.load(conf_fh)


def parse_fasta(fasta_path):
    return list(SeqIO.parse(fasta_path, 'fasta'))


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
    # claxification = {**classification, **taxification}
    return sequence_id, {**classification, **taxification} # merge dicts


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


#---------------------------------------------------------------------------------------------------


def kmer_lca_records(fasta_path,
                     one_codex_api_key=None,
                     progress: 'show progress bar (sent to stderr)' = False):
    '''
    Parallel lowest common ancestor sequence (LCA) classification using the One Codex API
    Returns Biopython SeqRecords with tictax annotations as the `description` attribute
    LCAs are assigned using an LCA index of 31mers from the One Codex database
    '''
    records = parse_fasta(fasta_path)
    one_codex_api_key = one_codex_api_key if one_codex_api_key else config()['one_codex_api_key']
    print('Classifying sequences…', file=sys.stderr)
    records = asyncio.get_event_loop().run_until_complete(oc_classify(records,
                                                            one_codex_api_key,
                                                            progress,
                                                            False))
    return records


def blast_lca(fasta_path, progress: 'show progress bar (sent to stderr)' = False):
    pass


def sort():
    '''stub'''
    pass


def rank_abundance():
    '''stub'''
    pass


def blast_lca_offline(fasta_path, classifications_path, acc2tax_path):
    '''stub'''
    pass


def plot(fasta_path):
    '''stub'''
    pass


if __name__ == '__main__':
    tictax.cli.main()
