#!/usr/bin/env python3
# Author: Bede Constantinides - b|at|bede|dot|im
# Usage: onecodex_filter_superkingdom.py path/to/file.fasta

import sys
import json
import logging
import concurrent.futures
from collections import OrderedDict

import argh
import requests
from Bio import SeqIO
from Bio import SeqUtils

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)

def onecodex_lca(seq, api_key):
    '''
    Returns dict of OneCodex real-time API k-mer hits for a given sequence
    e.g. {'elapsed_secs':'0.0003','k': 31,'n_hits': 97,'n_lookups': 128,'tax_id': 9606}
    '''
    url = 'https://app.onecodex.com/api/v0/search'
    payload = {'sequence':str(seq)}
    auth = requests.auth.HTTPBasicAuth(api_key, '')
    response = requests.post(url, payload, auth=auth, timeout=5)
    result = json.loads(response.text)
    return result

def ebi_taxid_to_lineage(tax_id):
    '''
    Returns scientific name and lineage for a given taxid using EBI's taxonomy API
    e.g.('Retroviridae', ['Viruses', 'Retro-transcribing viruses'])
    '''
    url = 'http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{}'
    if tax_id == 0 or tax_id == 1:
        return None, None
    response = requests.get(url.format(tax_id), timeout=5)
    result = json.loads(response.text)
    sciname = result['scientificName']
    taxonomy = [x for x in result['lineage'].split('; ') if x]
    return sciname, taxonomy

def onecodex_lca_taxa(seqrecord, api_key):
    '''
    Returns a scientific name and lineage for a SeqRecord using OneCodex and EBI APIs
    e.g. ('NODE_3_length_4481_cov_46.6129_ID_7947',
     ('Hepatitis C virus genotype 3',
      ['Viruses',
       'ssRNA viruses',
       'ssRNA positive-strand viruses, no DNA stage',
       'Flaviviridae',
       'Hepacivirus']))
    '''
    hits = onecodex_lca(str(seqrecord.seq), api_key)
    sciname, taxonomy = ebi_taxid_to_lineage(hits['tax_id'])
    result = (sciname, taxonomy, hits)
    return result

def fasta_onecodex_lca_taxa(fasta_path, api_key):
    '''
    Executes onecodex_lca_taxa() in parallel for a multifasta file
    '''
    seqrecords = SeqIO.parse(fasta_path, 'fasta')
    taxa = {}
    with concurrent.futures.ThreadPoolExecutor(50) as executor:
        futures = {executor.submit(onecodex_lca_taxa, seqrecord, api_key): seqrecord for seqrecord in seqrecords}
        for future in concurrent.futures.as_completed(futures):
            seqrecord = futures[future]
            try:
                data = future.result()
            except Exception as exception:
                taxa[seqrecord.id] = (None, None, None)
                logger.info('Skipping '.format(seqrecord.id))
            else:
                taxa[seqrecord.id] = future.result()
    return taxa

def seqrecords(fasta_path):
    '''
    Accepts path to multifasta, returns list of Biopython SeqRecords
    '''
    return SeqIO.parse(fasta_path, 'fasta')

def filter_superkingdom(seq_ids, taxa, superkingdom='Viruses'):
    '''
    Returns list of sequence IDs matching a specified superkingdom
    '''
    filtered_seq_ids = []
    for seq_id in seq_ids:
        if taxa[seq_id][1] and taxa[seq_id][1][0] == superkingdom or not taxa[seq_id][1]:
            filtered_seq_ids.append(seq_id)
    return filtered_seq_ids

    '''
    Accepts Biopython SeqRecords and OneCodex taxon assignments and returns RGB colour strings
    Colours entities belonging to the viral superkingdom
    '''
    colours = OrderedDict()
    for record in seqrecords:
        if taxa[record.id] and taxa[record.id][0]:
            if taxa[record.id][1][0] == 'Eukaryota':
                colours[record.id] = 'rgba(0, 255, 0, 0.75)'
            elif taxa[record.id][1][0] == 'Bacteria':
                colours[record.id] = 'rgba(0, 0, 255, 0.75)'
            elif taxa[record.id][1][0] == 'Viruses':
                colours[record.id] = 'rgba(255, 0, 0, 0.75)'
            else:
                colours[record.id] = 'rgba(0, 0, 0, 0.75)'
        else:
            colours[record.id] = 'rgba(0, 0, 0, 0.25)'
    return colours

def subset_to_annotated_fasta(subset_ids, seqrecords, taxa):
    '''
    Writes fasta file containing a subset of Biopython SeqRecords
    '''
    subset_seqrecords = []
    for record in seqrecords:
        if record.id in subset_ids:
            if taxa[record.id][0]:
                record.description = '| ' + taxa[record.id][0]
            subset_seqrecords.append(record)
    SeqIO.write(subset_seqrecords, sys.stdout, 'fasta')


def main(api_key='0f63476dfb8d4b5d95c96bc96af70d7d',
         superkingdom='Viruses',
         fasta_path):
    contigs = list(seqrecords(fasta_path))
    contig_ids = [r.id for r in contigs]
    contig_taxa = fasta_onecodex_lca_taxa(fasta_path, api_key)
    viral_ids = filter_superkingdom(contig_ids, contig_taxa, superkingdom='Viruses',)
    subset_to_annotated_fasta(viral_ids, contigs, contig_taxa)

argh.dispatch_command(main)