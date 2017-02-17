import os

from tictax import tictax

def test_cli():
    os.system('tictax -h')

def test_kmer_lca(capsys):
    test_fa_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test.fa')
    tictax.kmer_lca(test_fa_path)
    stdout, stderr = capsys.readouterr()
    print(stdout.count('>'))
    assert stdout.count('>') == 30
