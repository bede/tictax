import sys
import argh
import asyncio

from Bio import SeqIO

from tictax import tictax


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


def main():
    parser = argh.ArghParser()
    parser.add_commands([kmer_lca])
    parser.dispatch()


if __name__ == '__main__':
    main()