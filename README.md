# tictax: streaming sequence classification with web services

Rapid lowest common ancestor (LCA) assignment from command line or Python using One
Codex and EBI APIs. Capable of hundreds of requests per second with short queries.  

Scientific names and lineages are retrieved through chained requests to the EBI 
taxonomy or a local `ete3` database if available, and annotated within fasta description fields for
easy viewing and parsing.

## Command line usage
```
$ tictax kmer-lca test.fa
```
![Tictax demo](demo.gif)
```
$ tictax kmer-lca --progress test.fa > test.tt.fa
```
![Tictax demo with progress](demo_progress.gif)  

Tictax uses the One Codex (fast) and/or EBI BLAST web services (accurate) to classify nucleotide sequences.
Taxonomy information is subsequrntly requested from the EBI Taxonomy API and/or a local 
Requests are made concurrently respecting rate limits, and are returned in the order they arrive.

## Python API usage
### `kmer_lca_records(fasta_path, one_codex_api_key, progress=False)`
- Returns Biopython SeqRecords with tictax annotations as the `description` attribute  
- LCAs are assigned using an LCA index of 31mers from the One Codex database
```python
import tictax

records = tictax.parse_fasta('test.fa') # Biopython SeqRecord generator 
records_classified = tictax.kmer_lca_records(records, one_codex_api_key) # List of SeqRecords
print(records_classified.format('fasta')) # Generate multifasta
```

## Installation
Requires Python 3.5+. Installs with pip. No external dependencies.
```
pip3 install tictax
```

## Issues
Feel free to open issues and PRs, else [tweet](https://twitter.com/beconstant) or mail me via `b àt bede dawt im`.

## Todo ✓📌
- [ ] Switch to semicolon delimiters
- [ ] Add ebiblastn LCA functionality
- [ ] Add in onecodexrt / ebiblastn
- [ ] Documentation
- [ ] `kmer-lca` command `--out` option for two column id to LCA info
- [ ] `abundance` command for generating taxonomic abundances matrices at a specified rank
- [ ] Stream input from stdin