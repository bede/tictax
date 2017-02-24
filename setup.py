import re
import sys

from setuptools import setup

__version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        open('tictax/__init__.py').read()).group(1)

if sys.version_info[0] < 3:
      sys.exit('Tictax requires Python 3.5 or greater. Please upgrade')

setup(name='tictax',
      version=__version__,
      description='Streaming nucleotide sequence classification with web services',
      url='https://github.com/bede/tictax',
      author='Bede Constantinides',
      author_email='bedeabc@gmail.com',
      license='LICENSE',
      packages=['tictax'],
      zip_safe=True,
      install_requires=['argh',
                        'tqdm',
                        'pandas',
                        'aiohttp',
                        'biopython'],
      package_data={'tictax': ['tests/test.fa']},
      entry_points = {'console_scripts':['tictax=tictax.cli:main']})
