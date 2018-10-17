import os
import sys
import gzip
import json
import tqdm
import ete3
import asyncio
import aiohttp
import warnings
import pandas as pd


from collections import defaultdict

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


def parse_seqs(seqs_path, fastq=False):
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


# --------------------------------------------------------------------------------------------------


def annotate_diamond(records, diamond_path):
    '''
    Retrieve scientific names and lineages for taxon IDs in Diamond output
    Returns taxonomically annotated SeqRecords with modified description attributes
    '''
    contigs_metadata = {}
    with open(diamond_path) as diamond_tax_fh:
        for line in diamond_tax_fh:
            contig, taxid, evalue = line.strip().split('\t')
            contigs_metadata[contig] = dict(taxid=int(taxid), evalue=float(evalue))

    ncbi = ete3.NCBITaxa()
    taxids = {m['taxid'] for m in contigs_metadata.values()} # set of taxids
    taxids_lineages = ncbi.get_lineage_translator(taxids)  # dict of taxid lists
    taxids_with_children = {x for v in taxids_lineages.values() for x in v}  # flatten to set
    taxids_names = ncbi.get_taxid_translator(taxids_with_children)
    taxids_ranks = ncbi.get_rank(taxids_with_children)

    for contig, md in contigs_metadata.items():
        md['sciname'] = taxids_names.get(md['taxid'])
        md['rank'] = taxids_ranks.get(md['taxid'], '')
        md['lineage_fmt'] = (':'.join([taxids_names.get(t, '')
                                       for t in taxids_lineages.get(md['taxid'], [])])
                                       if md['taxid'] else None)

    for r in records:
        md = contigs_metadata[r.id]
        r.description = f"{md['taxid']}|{md['rank']}|{md['sciname']}|{md['lineage_fmt']}"
    return records


def filter_taxa(records, taxids, unclassified=False, discard=False):
    '''
    Selectively include or discard specified taxon IDs from tictax annotated FASTA/Qs
    Filters all children of specified taxon IDs
    Returns subset of input SeqRecords
    Taxon IDs of 1 and 2 are considered unclassified 
    '''
    taxids = set(taxids)
    kept_records = []
    ncbi = ete3.NCBITaxa()
    unclassified_taxids = {0,1}
    for r in records:
        taxid, rank, sciname, lineage = r.description.strip().partition(' ')[2].split('|')
        taxid = int(taxid)
        lineage = set(ncbi.get_lineage(taxid) or [0])  # lineage defined twice?
        intersection = lineage & taxids
        if taxid in unclassified_taxids and unclassified:
            kept_records.append(r)
        elif intersection and not discard:
            kept_records.append(r)
        elif not intersection and discard and taxid not in unclassified_taxids:
            kept_records.append(r)

    return kept_records


# Single record version
# def matrix(records, scafstats_path):
#     '''
#     Generate taxonomic count matrix from tictax classified contigs
#     '''
#     ncbi = ete3.NCBITaxa()
#     contigs_lineages = {r.id: r.description.strip().partition(' ')[2].split('|')[2] for r in records}
#     df = pd.read_csv(scafstats_path, sep='\t').set_index('#name')
#     contigs_counts = pd.Series(df.assignedReads.values, index=df.index).to_dict()
#     lineages_counts = defaultdict(int)
    
#     for contig, lineage in contigs_lineages.items():
#         lineages_counts[lineage] += contigs_counts.get(contig, 0)

#     lineages_counts_sorted = dict(reversed(sorted(lineages_counts.items(), key=lambda x: x[1])))
#     print(pd.DataFrame(lineages_counts_sorted, index=[scafstats_path]).transpose().to_csv())


def matrix(records, scafstats_path):
    '''
    Generate taxonomic count matrix from BBMap scafstats output to tictax classified contigs
    '''
    ncbi = ete3.NCBITaxa()
    
    # # Single file version 
    # contigs_lineages = {r.id: r.description.strip().partition(' ')[2].split('|')[2] for r in records}
    # df = pd.read_csv(scafstats_path, sep='\t').set_index('#name')
    # contigs_counts = pd.Series(df.assignedReads.values, index=df.index).to_dict()
    # lineages_counts = defaultdict(int)
    
    # for contig, lineage in contigs_lineages.items():
    #     lineages_counts[lineage] += contigs_counts.get(contig, 0)

    # lineages_counts_sorted = dict(reversed(sorted(lineages_counts.items(), key=lambda x: x[1])))
    # print(pd.DataFrame(lineages_counts_sorted, index=[scafstats_path]).transpose().to_csv())

    scafstats_paths = {fn.replace('.scafstats', ''): f'{scafstats_path}/{fn}'
                       for fn in os.listdir(scafstats_path)
                       if fn.endswith('.scafstats')}
    contigs_lineages = {r.id: r.description.strip().partition(' ')[2].split('|')[3] for r in records}
    samples = []
    
    for scafstat, path in scafstats_paths.items():
        raw_df = pd.read_csv(path, sep='\t').set_index('#name')
        contigs_counts = pd.Series(raw_df.assignedReads.values, index=raw_df.index).to_dict()
        lineages_counts = defaultdict(int)
        for contig, lineage in contigs_lineages.items():
            lineages_counts[lineage] += contigs_counts.get(contig, 0)
        counts_df = pd.DataFrame(lineages_counts, index=[scafstat]).transpose()
        samples.append(counts_df)

    return pd.concat(samples, axis=1, sort=False).rename_axis('lineage')



if __name__ == '__main__':
    tictax.cli.main()