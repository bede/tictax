#!/usr/bin/env python3

# import ete3

from Bio import SeqIO

# ncbi = ete3.NCBITaxa()

# with open('/Users/Bede/Research/Notebooks/res/kaiju.out') as kaiju_file:
#     for line in kaiju_file:
#         seq_id = line.strip().split('\t')[1].split(' ')[0]
#         kaiju_lca = int(line.strip().split('\t')[2])
#         lineage = ncbi.get_lineage(kaiju_lca) if kaiju_lca else []
#         lineage_taxids_names = ncbi.get_taxid_translator(lineage) if lineage else []
#         lineage_names = [lineage_taxids_names[taxid] for taxid in lineage]
#         lineage_fmt = ':'.join(lineage_names)
#         name = lineage_names[-1] if lineage_names else ''
#         print(seq_id, kaiju_lca, name, lineage_fmt, sep='|')

records = list(SeqIO.parse('/Users/bede/Research/kaiju_parser/RFF_TrimmedReads.fna', 'fasta'))

for r in records:
	print(r)

# with open('/Users/bede/Research/kaiju_parser/kaiju.out') as kaiju_file:
#     for line in kaiju_file:
#         seq_id = line.strip().split('\t')[1].split(' ')[0]
#         kaiju_lca = int(line.strip().split('\t')[2])
#         lineage = ncbi.get_lineage(kaiju_lca) if kaiju_lca else []
#         lineage_taxids_names = ncbi.get_taxid_translator(lineage) if lineage else []
#         lineage_names = [lineage_taxids_names[taxid] for taxid in lineage]
#         lineage_fmt = ':'.join(lineage_names)
#         name = lineage_names[-1] if lineage_names else ''
#         print(seq_id, kaiju_lca, name, lineage_fmt, sep='|')