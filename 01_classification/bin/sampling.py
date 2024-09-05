#!/usr/bin/env python3

import sys, gzip, os
import random

TAXONID_BACTERIA=2
TAXONID_EUKARYOTA=2759

taxnodesfilename=sys.argv[1]

GOLD_infile_name=sys.argv[2]
BacMap_infile_name=sys.argv[3]
NCBI_infile_name=sys.argv[4]

sampling_outfile_name=sys.argv[5]



def open_un_compressed(filename):
  if filename.endswith(".gz"):
    return gzip.open(filename, "rt")
  return open(filename)

def store_in_list(mylist, mykey, myvalue):
  if len(mylist) <= mykey:
    mylist.extend([0]*(1+mykey-len(mylist)))
  mylist[mykey]=myvalue

def isbacterial(parent_by_nodeid, taxonid):
  return isInLineage(parent_by_nodeid, taxonid, TAXONID_BACTERIA)

def isInLineage(parent_by_nodeid, taxonid, lineageid):
  while True:
    if taxonid == lineageid:
      return True
    parentid=parent_by_nodeid[taxonid]
    if parentid==taxonid:
      break
    taxonid=parentid
  return False


#reads taxonomy-nodes file
parent_by_nodeid=[]
with open_un_compressed(taxnodesfilename) as infile:
  for line in infile:
    parts = line.split("|",3)
    nodeid=int(parts[0].strip())
    parentid=int(parts[1].strip())
    store_in_list(parent_by_nodeid, nodeid, parentid)
sys.stderr.flush()



#reading GOLD dataset to extract Bacteria and their ecosystem
GOLD_metadata={}
with open_un_compressed(GOLD_infile_name) as infile:
  for line in infile:
    (GoldOrganismID, OrganismName, Phylum, NcbiTaxonomyID, Ecosystem, EcosystemCategory, HostTaxonomyID) = line[:-1].split("\t")
    taxonomyid=int(NcbiTaxonomyID.split(">")[1].split("<")[0])
    if HostTaxonomyID == "":
      hostID=0
    else:
      hostID=int(HostTaxonomyID)
    if isbacterial(parent_by_nodeid, taxonomyid):
      if Ecosystem != "":
        GOLD_metadata[OrganismName]=[Ecosystem, hostID]
    

#reading BacMap dataset to extract Bacteria and their habitat
BacMap_metadata={}
with open_un_compressed(BacMap_infile_name) as infile:
  next(infile)
  for line in infile:
    (name,NcbiTaxonomyID,habitat,hostID,hostname) = line.strip().split(",")
    taxonomyid=int(NcbiTaxonomyID)
    hostID=int(hostID)
    if isbacterial(parent_by_nodeid, taxonomyid):
      if habitat not in ("Multiple", "Specialized", "NA"):
        BacMap_metadata[name]=[habitat, hostID]


#reading NCBI dataset to extract Bacteria and their environment
NCBI_metadata={}
with open_un_compressed(NCBI_infile_name) as infile:
  next(infile)
  for line in infile:
    (acession,name,NcbiTaxonomyID,environment,hostname,hostID) = line.split("\t")
    taxonomyid=int(NcbiTaxonomyID)
    hostID=int(hostID)
    if isbacterial(parent_by_nodeid, taxonomyid):
      if environment not in ("NA", "Missing"):
        NCBI_metadata[name]=[environment, hostID]



#random sampling of 50 organisms for each dataset
with open(sampling_outfile_name, "w") as sampling_outfile:
  sampling_outfile.write("%s\t%s\t%s\t%s\n" % ("database", "name", "environmental_metadata", "hostID"))
  
  #sampling of the GOLD dataset
  GOLD_sample = random.sample(list(GOLD_metadata.keys()), 50)
  for sample in GOLD_sample:
    sampling_outfile.write("%s\t%s\t%s\t%i\n" % ("GOLD", sample, GOLD_metadata[sample][0], GOLD_metadata[sample][1]))
  #sampling of the BacMap dataset
  BacMap_sample = random.sample(list(BacMap_metadata.keys()), 50)
  for sample in BacMap_sample:
    sampling_outfile.write("%s\t%s\t%s\t%i\n" % ("BacMap", sample, BacMap_metadata[sample][0], BacMap_metadata[sample][1]))
  #sampling of the NCBI dataset
  NCBI_sample = random.sample(list(NCBI_metadata.keys()), 50)
  for sample in NCBI_sample:
    sampling_outfile.write("%s\t%s\t%s\t%i\n" % ("NCBI", sample, NCBI_metadata[sample][0], NCBI_metadata[sample][1]))

sys.stderr.flush()
