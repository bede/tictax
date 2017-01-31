#!/usr/bin/env python3

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
        async with session.post(url, data=payload, timeout=60) as response:
            r = await response.json()
    except (asyncio.TimeoutError, aiohttp.errors.ClientOSError):
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
        async with session.get(url.format(taxid), timeout=60) as response:
            r = await response.json()
            template['sciname'] = r.get('scientificName', '')
            template['lineage'] = r['lineage'].strip(' ;').split('; ') if 'lineage' in r else []
            template['rank'] = r.get('rank', '')
            return template
    except (asyncio.TimeoutError, aiohttp.errors.ClientOSError):
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
    conn = aiohttp.TCPConnector(limit=15)
    with aiohttp.ClientSession(auth=oc_auth, connector=conn) as oc_session:
        with aiohttp.ClientSession(connector=conn) as ebi_session:
            tasks = [classify_taxify(oc_session, ebi_session, r.id, str(r.seq)) for r in records]
            # responses = await asyncio.gather(*tasks)
    # return dict(responses)
            return [await f for f in tqdm.tqdm(asyncio.as_completed(tasks), total=len(tasks))]
            # return [await f for f in asyncio.as_completed(tasks)]

if __name__ == '__main__':
    prerequisites()
    records = fasta_seqrecords(sys.argv[1])
    print('Classifying sequences…', file=sys.stderr)
    asyncio.get_event_loop().run_until_complete(classify_taxify_records(records))
    print('✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌', file=sys.stderr)
