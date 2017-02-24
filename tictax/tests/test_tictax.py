import os

from tictax import tictax, cli


test_fa_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.fa') # 30 records


# ----- CLI ----- 
def test_cli():
    os.system('tictax -h')

def test_kmer_lca(capsys): # Captures and scrutinises stdout
    cli.kmer_lca(test_fa_path)
    stdout, stderr = capsys.readouterr()
    assert stdout.count('>') == 30


# ----- Python API -----
def test_kmer_lca_records():
    classified_records = tictax.kmer_lca_records(test_fa_path) # Assuming config is in place
    assert len(classified_records) == 30
    classified_records = tictax.kmer_lca_records(test_fa_path, # Ignore config
                                                 one_codex_api_key='0f63476dfb8d4b5d95c96bc96af70d7d')
    assert len(classified_records) == 30
