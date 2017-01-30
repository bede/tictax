#!/usr/bin/env python3

import sys
import json
import asyncio
import aiohttp

from Bio import SeqIO
from pprint import pprint
from collections import OrderedDict


def fasta_seqrecords(fasta_path):
    return list(SeqIO.parse(fasta_path, 'fasta'))


def async_onecodex_classification(records):
    async def fetch(session, name, sequence):
        url = 'https://app.onecodex.com/api/v0/search'
        payload = {'sequence': str(sequence)}
        async with session.post(url, data=payload) as response:
            return name, await response.json()
    
    async def loop():
        auth = aiohttp.BasicAuth('0f63476dfb8d4b5d95c96bc96af70d7d')
        conn = aiohttp.TCPConnector(limit=100)
        with aiohttp.ClientSession(auth=auth, connector=conn) as session:
            tasks = [fetch(session, record.id, record.seq) for record in records]
            responses = await asyncio.gather(*tasks)
        return OrderedDict(responses)

    return asyncio.get_event_loop().run_until_complete(loop())


def async_ebi_taxonomy(names_taxids):
    names_taxids_nz = OrderedDict((n,t if t else 1) for n,t in names_taxids.items()) # taxids of 0 become 1
    async def fetch(session, record_id, taxid):
        url = 'http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{}'
        async with session.get(url.format(taxid)) as response:
            json = await response.json() if 'json' in response.headers.get('content-type') else []
            return record_id, json

    async def loop():
        conn = aiohttp.TCPConnector(limit=30)
        with aiohttp.ClientSession(connector=conn) as session:
            tasks = [fetch(session, n, t) for n, t in names_taxids_nz.items()] 
            responses = await asyncio.gather(*tasks)
        return OrderedDict(responses)

    raw_results = asyncio.get_event_loop().run_until_complete(loop())
    names_taxonomy = OrderedDict((name, {'division':'', 'taxid':'', 'rank':'', 'sciname': '', 'lineage': ''}) for name in names_taxids_nz)
    for name, raw_result in raw_results.items():
        names_taxonomy[name]['division'] = raw_result['division'].strip() if 'division' in raw_result else ''
        names_taxonomy[name]['taxid'] = raw_result['taxId'].strip() if 'taxId' in raw_result else ''
        names_taxonomy[name]['rank'] = raw_result['rank'].strip() if 'rank' in raw_result and raw_result['rank'] != 'no rank' else ''
        names_taxonomy[name]['sciname'] = raw_result['scientificName'].strip() if 'scientificName' in raw_result and raw_result['scientificName'] != 'root' else ''
        names_taxonomy[name]['lineage'] = raw_result['lineage'].strip('; ').split('; ') if 'lineage' in raw_result else []
    return names_taxonomy


def annotated_fasta(records, assignments, taxonomies):
    for record in records:
        record.description = '|' + '|'.join((taxonomies[record.id]['division'],
                                             taxonomies[record.id]['taxid'],
                                             taxonomies[record.id]['rank'],
                                             taxonomies[record.id]['sciname'],
                                             ':'.join(taxonomies[record.id]['lineage']))) + '|'
    viral_records = []
    for record in records:
        if 'Virus' in record.description:
            viral_records.append(record)

    SeqIO.write(viral_records, sys.stdout, 'fasta')


records = fasta_seqrecords(sys.argv[1])

print('Classifying contigs', end=' ', file=sys.stderr)
assignments = async_onecodex_classification(records)
print('✓', file=sys.stderr)

assignment_taxids = OrderedDict((k, v['tax_id']) if 'tax_id' in v else (k, '0') for k, v in assignments.items())

print('Retrieving taxonomy information', end=' ', file=sys.stderr)
taxonomies = async_ebi_taxonomy(assignment_taxids)
print('✓', file=sys.stderr)

print('Writing annotated fasta', end=' ', file=sys.stderr)
annotated_fasta(records, assignments, taxonomies)
print('✓', file=sys.stderr)
