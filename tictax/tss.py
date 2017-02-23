import tictax

if __name__ == '__main__':
	records = tictax.parse_fasta('/Users/Bede/Research/Tools/tictax/tictax/tests/test.fa')
	records_classified = tictax.kmer_lca_stream(records, '0f63476dfb8d4b5d95c96bc96af70d7d')
	print(records_classified.format('fasta'))