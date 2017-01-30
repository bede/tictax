#!/usr/bin/env python3

import io
import sys
import time
import logging
import requests
import collections
import multiprocessing
import concurrent.futures

import ete3

from Bio import SeqIO

from pprint import pprint

ncbi = ete3.NCBITaxa()

logger = logging.getLogger(__name__)
# logging.basicConfig(level=logging.INFO)

def build_ebi_blast_query(title, sequence):
    return { 'email': 'bede.constantinides@manchester.ac.uk',
             'program': 'blastn',
             'stype': 'dna',
             'database': 'em_rel_vrl',
             'align': 6,
             'match_scores': '1,-3',
             'gapopen': 5, 
             'gapextend': 2,
             'exp': '1e-10',
             'filter': 'T',
             'dropoff': 0,
             'scores': 5,
             'alignments': 5,
             'title': title,
             'sequence': str(sequence) }

def ebi_blast(query):
    run_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
    status_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/'
    results_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/'
    dbfetch_url = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/{}/{}/fasta'
    
    def parse_hits(title, raw_hits):
        hits = []
        for line in io.StringIO(raw_hits):
            if ':' in line:
                fields = line.strip().split('\t')
                hit_db, hit_id = fields[1].split(':')
                hit_fa = requests.get(dbfetch_url.format(hit_db, hit_id))
                hit_fa_header = hit_fa.text.partition('\n')[0][1:].strip()
                hit_fa_accession, hit_fa_name = hit_fa_header.split(' ', 1)
                hit = (hit_fa_accession, hit_fa_name) + tuple(fields[2:])
                hits.append(hit)
        return hits

    results = {}
    call = requests.post(run_url, data=query)
    start_time = time.time()
    logger.info('dispatched blast jobid: ' + call.text)
    while True:
        status = requests.get(status_url + call.text)
        if status.text == 'FINISHED':
            raw_hits = requests.get(results_url + call.text + '/out')
            hits = parse_hits(query['title'], raw_hits.text)
            # results[query['title']] = hits
            logger.info(status.text + ' ' + call.text)
            print(time.time() - start_time)
            break
        elif time.time() - start_time > 120:
            print('blast timeout')
            logger.error('blast timeout')
            break
        elif status.text == 'RUNNING':
            time.sleep(5)
        else:
            print('status: ' + status.text)
            logger.error('status: ' + status.text)
            break
    
    # logger.info(results)
    return hits

def ebi_parallel_blast(fasta):
    records = collections.OrderedDict()
    with open(fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            records[record.id] = record.seq
    
    queries = [build_ebi_blast_query(title, seq) for title, seq in records.items()]
    hits = {}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as x:
        results = {q['title']:r for q, r in zip(queries, x.map(ebi_blast, queries))}

    return results

results = ebi_parallel_blast(sys.argv[1])
# pprint(results)

for query, hits in results.items():
    if hits:
        print(query, hits[0][1], sep=', ')
    else:
        print(query, None)

top_hits = {query:hits[1] for query, hits in results.items() if hits}
# top_hit_names = [{query: hits[1][1]} for query, hits in result.items() for result in results]

# print(top_hit_names)

