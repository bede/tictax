#!/usr/bin/env python3

import ete3

ncbi = ete3.NCBITaxa()

print(ncbi.get_taxid_translator([9606]))
# print(ncbi.get_lineage('9606'))
# print(ncbi.get_taxid_translator('9606'))
# print(ncbi.get_name_translator('9606'))
# print(ncbi.translate_to_names('9606'))