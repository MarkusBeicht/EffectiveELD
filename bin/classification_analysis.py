#!/usr/bin/env python3

import sys, os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_unweighted

class_per_DB_infile_name=sys.argv[1]


GOLD_class_by_taxonid={}
BacMap_class_by_taxonid={}
NCBI_class_by_taxonid={}


#reads the classification per database file and stores the NcbiTaxonomyIDs in lists
with open(class_per_DB_infile_name) as infile:
  for line in infile:
    (database, NcbiTaxonomyID, classification) = line.strip().split("\t")
    taxonomyid=int(NcbiTaxonomyID)
    if database == "GOLD":
      GOLD_class_by_taxonid[taxonomyid]=classification
    elif database == "BacMap":
      BacMap_class_by_taxonid[taxonomyid]=classification
    elif database == "NCBI":
      NCBI_class_by_taxonid[taxonomyid]=classification





#venn-diagram of NCBI-taxonomyIDs classified by the 3 databases GOLD, BacMap, NCBI
GOLD_set=set(GOLD_class_by_taxonid.keys())
BacMap_set=set(BacMap_class_by_taxonid.keys())
NCBI_set=set(NCBI_class_by_taxonid.keys())

venn3([GOLD_set, BacMap_set, NCBI_set], ('GOLD', 'BacMap', 'NCBI'))
plt.savefig('venn_taxonIDsperDB_weighted.png')
plt.close()

venn3_unweighted([GOLD_set, BacMap_set, NCBI_set], ('GOLD', 'BacMap', 'NCBI'))
plt.savefig('venn_taxonIDsperDB_unweighted.png')
plt.close()





#comparison of the overlapping NCBItaxonomyIDs between the different databases
#the results are written to stdoutput

sys.stderr.write("\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Databases", "shared_NCBItaxonomyIDs", "shared_classifications", "different_classifications", "overlap (%)", "first_DB_s", "first_DB_n"))

#GOLD_BacMap
shared_NCBItaxonomyIDs = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid}
shared_classifications = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and GOLD_class_by_taxonid[x] == BacMap_class_by_taxonid[x]}
different_classifications = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and GOLD_class_by_taxonid[x] != BacMap_class_by_taxonid[x]}
n_to_s = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and GOLD_class_by_taxonid[x] == "n" and BacMap_class_by_taxonid[x] == "s"}
s_to_n = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and GOLD_class_by_taxonid[x] == "s" and BacMap_class_by_taxonid[x] == "n"}
sys.stderr.write("%s\t%i\t%i\t%i\t%f\t%i\t%i\n" % ("GOLD_BacMap", len(shared_NCBItaxonomyIDs), len(shared_classifications), len(different_classifications), len(shared_classifications)/len(shared_NCBItaxonomyIDs), len(n_to_s) ,len( s_to_n)))



#GOLD_NCBI
shared_NCBItaxonomyIDs = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in NCBI_class_by_taxonid}
shared_classifications = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in NCBI_class_by_taxonid and GOLD_class_by_taxonid[x] == NCBI_class_by_taxonid[x]}
different_classifications = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in NCBI_class_by_taxonid and GOLD_class_by_taxonid[x] != NCBI_class_by_taxonid[x]}
n_to_s = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in NCBI_class_by_taxonid and GOLD_class_by_taxonid[x] == "n" and NCBI_class_by_taxonid[x] == "s"}
s_to_n = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in NCBI_class_by_taxonid and GOLD_class_by_taxonid[x] == "s" and NCBI_class_by_taxonid[x] == "n"}
sys.stderr.write("%s\t%i\t%i\t%i\t%f\t%i\t%i\n" % ("GOLD_NCBI", len(shared_NCBItaxonomyIDs), len(shared_classifications), len(different_classifications), len(shared_classifications)/len(shared_NCBItaxonomyIDs), len(n_to_s) ,len( s_to_n)))

#BacMap_NCBI
shared_NCBItaxonomyIDs = {x: BacMap_class_by_taxonid[x] for x in BacMap_class_by_taxonid if x in NCBI_class_by_taxonid}
shared_classifications = {x: BacMap_class_by_taxonid[x] for x in BacMap_class_by_taxonid if x in NCBI_class_by_taxonid and BacMap_class_by_taxonid[x] == NCBI_class_by_taxonid[x]}
different_classifications = {x: BacMap_class_by_taxonid[x] for x in BacMap_class_by_taxonid if x in NCBI_class_by_taxonid and BacMap_class_by_taxonid[x] != NCBI_class_by_taxonid[x]}
n_to_s = {x: BacMap_class_by_taxonid[x] for x in BacMap_class_by_taxonid if x in NCBI_class_by_taxonid and BacMap_class_by_taxonid[x] == "n" and NCBI_class_by_taxonid[x] == "s"}
s_to_n = {x: BacMap_class_by_taxonid[x] for x in BacMap_class_by_taxonid if x in NCBI_class_by_taxonid and BacMap_class_by_taxonid[x] == "s" and NCBI_class_by_taxonid[x] == "n"}
sys.stderr.write("%s\t%i\t%i\t%i\t%f\t%i\t%i\n" % ("BacMap_NCBI", len(shared_NCBItaxonomyIDs), len(shared_classifications), len(different_classifications), len(shared_classifications)/len(shared_NCBItaxonomyIDs), len(n_to_s) ,len( s_to_n)))



#GOLD_BacMap_NCBI
shared_NCBItaxonomyIDs = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and x in NCBI_class_by_taxonid}
shared_classifications = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and x in NCBI_class_by_taxonid and GOLD_class_by_taxonid[x] == BacMap_class_by_taxonid[x] and GOLD_class_by_taxonid[x] == NCBI_class_by_taxonid[x]}
different_classifications = {x: GOLD_class_by_taxonid[x] for x in GOLD_class_by_taxonid if x in BacMap_class_by_taxonid and x in NCBI_class_by_taxonid and (GOLD_class_by_taxonid[x] != BacMap_class_by_taxonid[x] or GOLD_class_by_taxonid[x] != NCBI_class_by_taxonid[x] or BacMap_class_by_taxonid[x] != NCBI_class_by_taxonid[x])}


sys.stderr.write("%s\t%i\t%i\t%i\t%f\t%s\t%s\n" % ("GOLD_BacMap_NCBI", len(shared_NCBItaxonomyIDs), len(shared_classifications), len(different_classifications), len(shared_classifications)/len(shared_NCBItaxonomyIDs), "", ""))


sys.stderr.write("\n%s\t%s\t%s\t%s\n" % ("NCBItaxonomyID", "GOLD", "BacMap", "NCBI"))
for NCBItaxonomyID in different_classifications.keys():
  sys.stderr.write("%i\t%s\t%s\t%s\n" % (NCBItaxonomyID, GOLD_class_by_taxonid[NCBItaxonomyID], BacMap_class_by_taxonid[NCBItaxonomyID], NCBI_class_by_taxonid[NCBItaxonomyID]))


sys.stderr.flush()
