# tictax: streaming sequence classification with web services
[![PyPI version](https://badge.fury.io/py/tictax.svg)](https://badge.fury.io/py/tictax)

Fast and convenient lowest common ancestor (LCA) assignment from command line or Python using One
Codex and EBI APIs. Capable of ~100 requests/s with the wind blowing in the right direction. 

Scientific names and lineages are retrieved using the EBI taxonomy API and annotated within fasta
description fields for easy viewing and parsing.  

**Tictax currently supports One Codex classification only. BLAST-based LCA (slow, sensitive) coming soon**

## Command line usage
```
$ tictax kmer-lca test.fa
```
![Tictax demo](demo.gif)
```
$ tictax kmer-lca --progress test.fa > test.tt.fa
```
![Tictax demo with progress](demo_progress.gif)  

### Built-in help
```
$ tictax -h
usage: tictax [-h] {kmer-lca} ...

positional arguments:
  {kmer-lca}
    kmer-lca  Lowest common ancestor sequence assignment using the One Codex
              API. Streams annotated records to stdout in fasta format. Taxa
              assigned using the One Codex 31mer LCA database.

optional arguments:
  -h, --help  show this help message and exit
```
```
$ tictax kmer-lca -h
usage: tictax kmer-lca [-h] [-p] fasta-path

    Lowest common ancestor sequence assignment using the One Codex API.
    Streams annotated records to stdout in fasta format.
    Taxa assigned using the One Codex 31mer LCA database.
    

positional arguments:
  fasta-path      path to fasta formatted input

optional arguments:
  -h, --help      show this help message and exit
  -p, --progress  show progress bar (sent to stderr) (default: False)
```

## Python API
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
Requires Python 3.5+. Installs with pip.
```
pip3 install tictax
```

## Issues
Feel free to open issues and PRs, else [tweet](https://twitter.com/beconstant) or mail me via `b àt bede dawt im`.

## Todo ✓📌
- [ ] `blast-lca` command - add BLAST LCA functionality
- [x] Switch to semicolon delimiters
- [x] Add in onecodexrt / ebiblastn
- [ ] `kmer-lca` command `--out` option for two column id to LCA info
- [ ] `abundance` command for generating taxonomic abundances matrices at a specified rank
- [ ] Stream input from stdin
- [ ] Use`aiohttp.errors` for exception handling
- [ ] Option to return records as async generator 
- [x] Support `fastq` and `gzip` input
