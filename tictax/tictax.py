#!/usr/bin/env python3

import os
import sys
import ete3
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



def prerequisites():
    if not os.path.exists('one_codex_api_key'):
        print('Tictax requires a One Codex API key')
        print('1) Sign up for a One Codex account at https://www.onecodex.com')
        print('2) Copy and paste your API key below')
        api_key = input('API key: ').strip()
        with open('one_codex_api_key', 'w') as api_key_file: # Needs storing in a sensible place... dotfile?
            api_key_file.write(api_key)
        print('Key saved')


def fasta_seqrecords(fasta_path):
    return list(SeqIO.parse(fasta_path, 'fasta'))


def prepare_record(sequence, id, classification):
    description = (str(classification['taxid'])
                   + '|' + classification['sciname']
                   + '|' + classification['rank']
                   + '|' + ':'.join(classification['lineage']))
    record = SeqRecord(Seq(sequence, IUPAC.ambiguous_dna), id=id, description=description)
    return record.format('fasta')


async def classify_oc(session, sequence_id, sequence):
    url = 'https://app.onecodex.com/api/v0/search'
    payload = {'sequence': str(sequence)[:25000]}
    try:
        async with session.post(url, data=payload, timeout=600) as response:
            r = await response.json()
    except (asyncio.TimeoutError, aiohttp.errors.ClientOSError, json.decoder.JSONDecodeError):
        r = {}
    return {'taxid': r.get('tax_id', 0),
            'k': r.get('k', ''),
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
    classification = await classify_oc(oc_session, sequence_id, sequence)
    taxid = classification.get('taxid')
    taxification = await taxify_ebi(ebi_session, sequence_id, taxid)
    claxification = {**classification, **taxification}
    print(prepare_record(sequence, sequence_id, claxification), end='')
    return sequence_id, {**classification, **taxification} # merge dicts


async def classify_taxify_records(records):
    oc_auth = aiohttp.BasicAuth('0f63476dfb8d4b5d95c96bc96af70d7d')
    conn = aiohttp.TCPConnector(limit=10)
    with aiohttp.ClientSession(auth=oc_auth, connector=conn) as oc_session:
        with aiohttp.ClientSession(connector=conn) as ebi_session:
            tasks = [classify_taxify(oc_session, ebi_session, r.id, str(r.seq)) for r in records]
            # responses = await asyncio.gather(*tasks)
    # return dict(responses)
            return [await f for f in tqdm.tqdm(asyncio.as_completed(tasks), total=len(tasks))]
            # return [await f for f in asyncio.as_completed(tasks)]


def kmer_lca(fasta_path, one_codex=True):
    prerequisites()
    records = fasta_seqrecords(fasta_path)
    print('Classifying sequences…', file=sys.stderr)
    asyncio.get_event_loop().run_until_complete(classify_taxify_records(records))
    print('✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌', file=sys.stderr)


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
        taxid = int(fields[0])
        if taxid > 0: # ETE can produce negative taxids
            lineage = set(ncbi.get_lineage(int(taxid)))
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
                r.description = '|{}|{}|{}|{}|'.format(taxid, sciname, rank, lineage_names_fmt) # names
            else:

                rank, sciname, lineage_fmt = '', '', ''
                r.description = '|{}|{}|{}|{}|'.format(taxid, sciname, rank, lineage_names_fmt) # names
            print(r.format('fasta'), end='')


parser = argh.ArghParser()
parser.add_commands([plot, kmer_lca, kmer_lca_offline, sort, filter, filter_old, annotate_megan_taxids])


if __name__ == '__main__':
    parser.dispatch()
