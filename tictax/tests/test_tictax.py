import os
import json

from tictax import tictax, cli


# Fetch config file including API key. Run `$ tictax kmer-lca` to generate
conf_path = os.path.join(os.path.expanduser('~'), '.tictax')
try:
    with open(conf_path, 'r') as conf_fh:
        conf = json.load(conf_fh)
except Exception:
    raise IOError('Error parsing {}'.format(conf_path))


test_fa_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.fa') # 30 records
test_fa_gz_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.fa.gz') # 30 records
test_fq_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.fq') # 10 records
test_fq_gz_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.fq.gz') # 10 records


# ----- CLI ----- 
def test_cli():
    os.system('tictax -h')

def test_kmer_lca_fa(capsys): # Captures and scrutinises stdout
    cli.kmer_lca(test_fa_path)
    stdout, stderr = capsys.readouterr()
    assert stdout.count('>') == 30

def test_kmer_lca_fa_gz(capsys): # Captures and scrutinises stdout
    cli.kmer_lca(test_fa_gz_path)
    stdout, stderr = capsys.readouterr()
    assert stdout.count('>') == 30

def test_kmer_lca_fq(capsys): # Captures and scrutinises stdout
    cli.kmer_lca(test_fq_path)
    stdout, stderr = capsys.readouterr()
    assert stdout.count('>') == 10

def test_kmer_lca_fq_gz(capsys): # Captures and scrutinises stdout
    cli.kmer_lca(test_fq_gz_path)
    stdout, stderr = capsys.readouterr()
    assert stdout.count('>') == 10


# ----- Python API -----
def test_kmer_lca_records():
    classified_records = tictax.kmer_lca_records(test_fa_path) # Assumes config is in place
    assert len(classified_records) == 30
    classified_records = tictax.kmer_lca_records(test_fa_path,
                                                 one_codex_api_key=conf['one_codex_api_key'])
    assert len(classified_records) == 30
