#!/usr/bin/env python3

import io
import sys
import time
import datetime
import requests

from Bio import SeqIO

run_endpoint = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
status_endpoint = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/'
results_endpoint = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/'

query = { 'email': 'bede.constantinides@manchester.ac.uk',
          'program': 'blastn',
          'stype': 'dna',
          'database': 'em_rel_vrl',
          'align': 6,
          'sequence': None,
          'match_scores': '1,-3',
          'gapopen': 5, 
          'gapextend': 2,
          'exp': '1e-3',
          'filter': 'T',
          'dropoff': 0,
          'scores': 5,
          'alignments': 5,
          'title': None,
          'sequence': None }

records = []
with open(sys.argv[1], 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        records.append(record.format('fasta'))

hits = {}
for record in records:
    query['sequence'] = record
    query['title'] = record.split('\n')[0][1:]
    hits[query['title']] = []
    job_id = requests.post(run_endpoint, data=query)
    print('# db: {}, jobid: {}, time: {}'
          .format(query['database'], job_id.text, datetime.datetime.now()))
    print(query['title'])
    while True:
        status = requests.get(status_endpoint + job_id.text)
        if status.text == 'FINISHED':
            results = requests.get(results_endpoint + job_id.text + '/out')
            for line in io.StringIO(results.text):
                if ':' in line:
                    fields = line.strip().split('\t')
                    hit_db, hit_id = fields[1].split(':')
                    hit_fasta = requests.get('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/{}/{}/fasta'
                                             .format(hit_db, hit_id))
                    hit_fasta_header= hit_fasta.text.partition('\n')[0][1:100]
                    hit = hit_fasta_header + '\t' + ', '.join(fields[2:])
                    if hit_fasta_header not in hits[query['title']]:
                        hits[query['title']].append(hit_fasta_header)
                        print(hit)
            break
        elif status.text == 'RUNNING':
            time.sleep(5)
        else:
            print(status.text)
            break
