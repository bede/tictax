import os
import sys
# import ete3
import json
import tqdm
import argh
import asyncio
import aiohttp
import sqlite3

import pandas as pd

import plotly.offline as py
import plotly.graph_objs as go

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


tick_tacks = '✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌'

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
    description = (str(classification['taxid'])
                   + '|' + classification['sciname']
                   + '|' + classification['rank']
                   + '|' + ':'.join(classification['lineage']))
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
    return {'sequence': sequence,
            'k': r.get('k', ''),
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
    claxification = {**classification, **taxification}
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


# async def oc_classify(records, one_codex_api_key, progress=False, stdout=False):
#     oc_auth = aiohttp.BasicAuth(one_codex_api_key)
#     conn = aiohttp.TCPConnector(limit=10)
#     with aiohttp.ClientSession(auth=oc_auth, connector=conn) as oc_session:
#         with aiohttp.ClientSession(connector=conn) as ebi_session:
#             tasks = [classify_taxify(oc_session, ebi_session, r.id, str(r.seq)) for r in records]
#             responses = await asyncio.gather(*tasks)
#     return responses


def kmer_lca(fasta_path,
             progress: 'show progress bar (sent to stderr)' = False):
    '''
    Streaming lowest common ancestor sequence (LCA) classification using the One Codex API
    Streams LCA-annotated records to stdout in fasta format
    LCAs are assigned using an LCA index of 31mers from the One Codex database
    '''
    conf = config()
    records = parse_fasta(fasta_path)
    print('Classifying sequences…', file=sys.stderr)
    asyncio.get_event_loop().run_until_complete(oc_classify(records,
                                                            conf['one_codex_api_key'],
                                                            progress,
                                                            True))
    print('✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌', file=sys.stderr)


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

#---------------------------------------------------------------------------------------------------

def filter(fasta_path: 'path to tictax annotated fasta file',
           taxids: 'single or multiple comma delimited taxon IDs' = '1',
           exclude: 'exclude matched records' = False):
    '''filter a tictax annotated fasta file by taxon ID(s)'''
    records = SeqIO.parse(fasta_path, 'fasta')
    filter_taxids = set(map(int, taxids.split(',')) if ',' in taxids else [int(taxids)]) # As set
    ncbi = ete3.NCBITaxa()
    for r in records:
        match = False
        description = r.description.partition(' ')[2]
        fields = tuple(description.split('|')[1:-1]) # remove bookend pipes
        # print(fields)
        taxid = int(fields[0])
        if taxid >= 0: # ETE can produce negative taxids
            lineage = set(ncbi.get_lineage(int(taxid))) if taxid else set([0])
            if set.intersection(lineage, filter_taxids):
                match = True
            if exclude: # exclude matched records
                if match:
                    continue
                else:
                    print(r.format('fasta'), end='')
            else: # include matched records
                if match:
                    print(r.format('fasta'), end='')
                else:
                    continue



def filter_old(fasta_path, taxid=None, unclassified=False):
    ncbi = ete3.NCBITaxa()
    bad_taxids = set([9606, 131567, 2759])
    bad_superkingdoms = set([2, 2157, 2759])


    if unclassified:
        records = SeqIO.parse(fasta_path, 'fasta')
        for r in records:
            desc = r.description.partition(' ')[2]
            if desc:
                taxid = max(int(desc.strip('|').split('|')[0]), 1)
                lineage_taxids = ncbi.get_lineage(taxid)
                if len(lineage_taxids) >= 3:
                    if lineage_taxids[2] not in bad_superkingdoms:
                        if taxid not in bad_taxids:
                            print(r.format('fasta'), end='')
                else:
                    print(r.format('fasta'), end='')


def sort():
    pass

def kmer_lca_online():
    pass

def blast_lca_offline(fasta_path, classifications_path, acc2tax_path):

    pass

def kmer_lca_offline(fasta_path, classifications_path, lineage=False):
    ncbi = ete3.NCBITaxa()
    records = SeqIO.parse(fasta_path, 'fasta')
    lcas = pd.read_csv(classifications_path, sep='\t', usecols=[0,1])
    lcas['Header'] = lcas['Header'].map(lambda x: x.lstrip('>'))
    ids_lcas = dict(lcas.to_records(index=False))
    # amended_records = []
    for r in records:
        taxid = int(ids_lcas[r.id])
        if taxid:
            sciname = ncbi.translate_to_names([taxid])[0]
            lineage_taxids = ncbi.get_lineage(taxid)
            lineage_taxids_names = ncbi.get_taxid_translator(lineage_taxids)
            lineage_fmt = ':'.join([lineage_taxids_names[t] for t in lineage_taxids])
        else:
            sciname, lineage_fmt = '', ''
        r.description = '|{}|{}|{}|'.format(taxid, sciname, lineage_fmt)
        # amended_records.append(r)
        print(r.format('fasta'), end='')
    # SeqIO.write(amended_records, sys.stdout, 'fasta')

def plot(fasta_path):
    pass

def annotate_megan_taxids(fasta_path, megan_csv_path):
    ncbi = ete3.NCBITaxa()
    records = SeqIO.parse(fasta_path, 'fasta')
    ids_taxids = dict(pd.read_csv(megan_csv_path, sep='\t').to_records(index=False))
    for r in records:
        if r.id in ids_taxids:
            taxid = ids_taxids[r.id]
            if taxid > 0:
                sciname = ncbi.translate_to_names([taxid])[0]
                lineage_taxids = ncbi.get_lineage(taxid)
                rank = ncbi.get_rank([taxid]).get(taxid, 'unknown')
                lineage_taxids_names = ncbi.get_taxid_translator(lineage_taxids)
                lineage_names_fmt = ':'.join([lineage_taxids_names.get(t, 'unknown') for t in lineage_taxids])
                r.description = '|{}|{}|{}|{}|'.format(taxid, sciname, rank, lineage_names_fmt)
            else:
                r.description = '|0|unknown|unknown|unknown|'
            print(r.format('fasta'), end='')


parser = argh.ArghParser()
parser.add_commands([plot,
                     kmer_lca,
                     kmer_lca_records,
                     sort,
                     filter,
                     filter_old, 
                     annotate_megan_taxids])


if __name__ == '__main__':
    parser.dispatch()
