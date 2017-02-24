import sys
import argh
import asyncio

import pandas as pd

from Bio import SeqIO

from tictax import tictax


def kmer_lca(fasta_path: 'path to fasta formatted input',
             progress: 'show progress bar (sent to stderr)' = False):
    '''
    Lowest common ancestor sequence assignment using the One Codex API.
    Streams annotated records to stdout in fasta format.
    Taxa assigned using the One Codex 31mer LCA database.
    '''
    conf = tictax.config()
    records = tictax.parse_fasta(fasta_path)
    print('Classifying sequences…', file=sys.stderr)
    asyncio.get_event_loop().run_until_complete(tictax.oc_classify(records,
                                                                   conf['one_codex_api_key'],
                                                                   progress,
                                                                   True))
    print('✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌', file=sys.stderr)


# def filter_taxa(fasta_path: 'path to tictax annotated fasta file',
#            taxids: 'single or multiple comma delimited taxon IDs' = '1',
#            exclude: 'exclude matched records' = False):
#     '''filter a tictax annotated fasta file by taxon ID(s)'''
#     records = SeqIO.parse(fasta_path, 'fasta')
#     filter_taxids = set(map(int, taxids.split(',')) if ',' in taxids else [int(taxids)]) # As set
#     ncbi = ete3.NCBITaxa()
#     for r in records:
#         match = False
#         description = r.description.partition(' ')[2]
#         fields = tuple(description.split('|')[1:-1]) # remove bookend pipes
#         # print(fields)
#         taxid = int(fields[0])
#         if taxid >= 0: # ETE can produce negative taxids
#             lineage = set(ncbi.get_lineage(int(taxid))) if taxid else set([0])
#             if set.intersection(lineage, filter_taxids):
#                 match = True
#             if exclude: # exclude matched records
#                 if match:
#                     continue
#                 else:
#                     print(r.format('fasta'), end='')
#             else: # include matched records
#                 if match:
#                     print(r.format('fasta'), end='')
#                 else:
#                     continue


# def kmer_lca_offline(fasta_path, classifications_path, lineage=False):
#     ncbi = ete3.NCBITaxa()
#     records = SeqIO.parse(fasta_path, 'fasta')
#     lcas = pd.read_csv(classifications_path, sep='\t', usecols=[0,1])
#     lcas['Header'] = lcas['Header'].map(lambda x: x.lstrip('>'))
#     ids_lcas = dict(lcas.to_records(index=False))
#     # amended_records = []
#     for r in records:
#         taxid = int(ids_lcas[r.id])
#         if taxid:
#             sciname = ncbi.translate_to_names([taxid])[0]
#             lineage_taxids = ncbi.get_lineage(taxid)
#             lineage_taxids_names = ncbi.get_taxid_translator(lineage_taxids)
#             lineage_fmt = ':'.join([lineage_taxids_names[t] for t in lineage_taxids])
#         else:
#             sciname, lineage_fmt = '', ''
#         r.description = '|{}|{}|{}|'.format(taxid, sciname, lineage_fmt)
#         # amended_records.append(r)
#         print(r.format('fasta'), end='')
#     # SeqIO.write(amended_records, sys.stdout, 'fasta')


# def annotate_megan_taxids(fasta_path, megan_csv_path):
#     ncbi = ete3.NCBITaxa()
#     records = SeqIO.parse(fasta_path, 'fasta')
#     ids_taxids = dict(pd.read_csv(megan_csv_path, sep='\t').to_records(index=False))
#     for r in records:
#         if r.id in ids_taxids:
#             taxid = ids_taxids[r.id]
#             if taxid > 0:
#                 sciname = ncbi.translate_to_names([taxid])[0]
#                 lineage_taxids = ncbi.get_lineage(taxid)
#                 rank = ncbi.get_rank([taxid]).get(taxid, 'unknown')
#                 lineage_taxids_names = ncbi.get_taxid_translator(lineage_taxids)
#                 lineage_names_fmt = ':'.join([lineage_taxids_names.get(t, 'unknown') for t in lineage_taxids])
#                 r.description = '|{}|{}|{}|{}|'.format(taxid, sciname, rank, lineage_names_fmt)
#             else:
#                 r.description = '|0|unknown|unknown|unknown|'
#             print(r.format('fasta'), end='')


def main():
    parser = argh.ArghParser()
    parser.add_commands([kmer_lca])
    # parser.add_commands([kmer_lca,
    #                      filter_taxa,
    #                      kmer_lca_offline,
    #                      annotate_megan_taxids])
    parser.dispatch()



if __name__ == '__main__':
    main()