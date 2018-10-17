import sys
import argh
import asyncio
import warnings

from Bio import SeqIO
import pandas as pd

from tictax import tictax


def configure_warnings(show_warnings):
    '''Show or suppress warnings, mainly for TreeSwift Tree.mrca() operations'''
    if show_warnings:
        warnings.filterwarnings('always')
    else:
        warnings.filterwarnings('ignore')


def kmer_lca(seqs_path: 'path to (optionally gzipped) fasta/fastq input',
             fastq: 'input is fastq; disable autodetection' = False,
             progress: 'show progress bar (sent to stderr)' = False):
    '''
    Lowest common ancestor sequence assignment using the One Codex API.
    Streams annotated records to stdout in fasta format.
    Taxa assigned using the One Codex 31mer LCA database.
    '''
    conf = tictax.config()
    records = tictax.parse_seqs(seqs_path, fastq)
    print('Classifying sequences…', file=sys.stderr)
    asyncio.get_event_loop().run_until_complete(tictax.oc_classify(records,
                                                                   conf['one_codex_api_key'],
                                                                   progress,
                                                                   True))
    print('✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌 ✓📌', file=sys.stderr)


def annotate_diamond(fasta_path: 'path to fasta input',
                     diamond_path: 'path to Diamond taxonomic classification output'):
    '''
    Annotate fasta headers with taxonomy information from Diamond
    '''
    records = tictax.parse_seqs(fasta_path)
    annotated_records = tictax.annotate_diamond(records, diamond_path)
    SeqIO.write(annotated_records, sys.stdout, 'fasta')


@argh.named('filter')  # Avoids namespace collision in CLI
def filter_taxa(fasta_path: 'path to fasta input',
                taxids: 'comma delimited list of taxon IDs',
                unclassified: 'pass sequences unclassified at superkingdom level >(0)' = False,
                discard: 'discard specified taxa' = False,
                warnings: 'show warnings' = False):
    '''
    Customisable filtering of tictax flavoured fasta files
    '''
    configure_warnings(warnings)
    records = SeqIO.parse(fasta_path, 'fasta')
    filtered_records = tictax.filter_taxa(records,
                                          map(int, taxids.split(',')),
                                          unclassified,
                                          discard)
    SeqIO.write(filtered_records, sys.stdout, 'fasta')


def matrix(fasta_path: 'path to tictax annotated fasta input',
           scafstats_path: 'path to BBMap scaftstats file'):
    '''
    Generate taxonomic count matrix from tictax classified contigs
    '''
    records = SeqIO.parse(fasta_path, 'fasta')
    df = tictax.matrix(records, scafstats_path)
    df.to_csv(sys.stdout)



def main():
    argh.dispatch_commands([kmer_lca,
                            annotate_diamond,
                            filter_taxa,
                            matrix])


if __name__ == '__main__':
    main()