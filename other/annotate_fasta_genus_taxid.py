#!/usr/bin/env python3
# Usage: annotate_fasta_genus_taxid.py path/to/file.fasta

import sys
import json
import requests
from Bio import SeqIO

# This is my private API key!
onecodex_api_key = '0f63476dfb8d4b5d95c96bc96af70d7d'

def seqrecords(fasta_path):
    return SeqIO.parse(fasta_path, 'fasta')

def onecodex_lca(seq, onecodex_api_key):
    url = 'https://app.onecodex.com/api/v0/search'
    payload = {'sequence':str(seq)}
    auth = requests.auth.HTTPBasicAuth(onecodex_api_key, '')
    response = requests.post(url, payload, auth=auth, timeout=5)
    result = json.loads(response.text)
    return result

def ebi_taxid_to_lineage(tax_id):
    url = 'http://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{}'
    if tax_id == 0 or tax_id == 1:
        return None, None
    response = requests.get(url.format(tax_id), timeout=5)
    result = json.loads(response.text)
    sciname = result['scientificName']
    if 'lineage' in result:
        taxonomy = [x for x in result['lineage'].split('; ') if x]
    else:
        taxonomy = ('', ['', ''])
    return sciname, taxonomy

for seqrecord in seqrecords(sys.argv[1]):
	lca = onecodex_lca(str(seqrecord.seq), onecodex_api_key)
	taxonomy = ebi_taxid_to_lineage(lca['tax_id'])
	print(seqrecord.id, lca['tax_id'], taxonomy[0], sep='\t')