#!/usr/bin/env python3

import os
import sys

from Bio import SeqIO

for f in os.listdir():
	if f.endswith('.fasta'):
		print(f)
		records = SeqIO.parse(f, 'fasta')
		for record in records:
			print(len(record.seq))
