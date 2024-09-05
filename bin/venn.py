#!/usr/bin/env python3

import sys, os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

class_per_DB_infile_name=sys.argv[1]


GOLD_taxonIDs=[]
BacMap_taxonIDs=[]
NCBI_taxonIDs=[]

#reads the classification per database file and stores the NcbiTaxonomyIDs in lists
with open(class_per_DB_infile_name) as infile:
  for line in infile:
    (database, NcbiTaxonomyID, classification) = line.strip().split("\t")
    taxonomyid=int(NcbiTaxonomyID)
    if database == "GOLD":
      GOLD_taxonIDs.append(taxonomyid)
    elif database == "BacMap":
      BacMap_taxonIDs.append(taxonomyid)
    elif database == "NCBI":
      NCBI_taxonIDs.append(taxonomyid)


venn3([GOLD_taxonIDs, BacMap_taxonIDs, NCBI_taxonIDs], ('GOLD', 'BacMap', 'NCBI'))
plt.savefig('venn_taxonIDs.png')

sys.stderr.flush()
