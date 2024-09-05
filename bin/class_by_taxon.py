#!/usr/bin/env python3

import sys, gzip, os
import matplotlib.pyplot as plt


TAXONID_BACTERIA=2
TAXONID_EUKARYOTA=2759

taxnodesfilename=sys.argv[1]

GOLD_infile_name=sys.argv[2]
BacMap_infile_name=sys.argv[3]
NCBI_infile_name=sys.argv[4]

BacMap_classification_infile=sys.argv[5]
NCBI_classification_infile=sys.argv[6]
RefeseqAssembliesInput=sys.argv[7]

class_outfile_name=sys.argv[8]
classDB_outfile_name=sys.argv[9]
genome_outfile_name=sys.argv[10]




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

def iseukaryotic(parent_by_nodeid, taxonid):
  return isInLineage(parent_by_nodeid, taxonid, TAXONID_EUKARYOTA)

def isInLineage(parent_by_nodeid, taxonid, lineageid):
  while True:
    if taxonid == lineageid:
      return True
    parentid=parent_by_nodeid[taxonid]
    if parentid==taxonid:
      break
    taxonid=parentid
  return False

parent_by_nodeid=[]
nodecounter=0
sys.stderr.write("Reading taxonomy nodes")
with open_un_compressed(taxnodesfilename) as infile:
  for line in infile:
    if nodecounter % 100000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
    nodecounter += 1
    parts = line.split("|",3)
    nodeid=int(parts[0].strip())
    parentid=int(parts[1].strip())
    store_in_list(parent_by_nodeid, nodeid, parentid)
sys.stderr.write("done (%i nodes).\n" % (nodecounter))
sys.stderr.flush()




GOLD_class_by_taxonid={}
genomescounter=0
sys.stderr.write("Reading GOLD genomes")
with open_un_compressed(GOLD_infile_name) as infile:
  for line in infile:
    if genomescounter % 1000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
    genomescounter += 1
    (GoldOrganismID, OrganismName, Phylum, NcbiTaxonomyID, Ecosystem, EcosystemCategory, HostTaxonomyID) = line[:-1].split("\t")
    taxonomyid=int(NcbiTaxonomyID.split(">")[1].split("<")[0])
    if isbacterial(parent_by_nodeid, taxonomyid):
      genomeclass=""
      if not HostTaxonomyID and Ecosystem == "Environmental":
        genomeclass="n"
      elif (HostTaxonomyID and iseukaryotic(parent_by_nodeid, int(HostTaxonomyID))) or (Ecosystem == "Host-associated" and EcosystemCategory not in ["", "Microbial", "Unclassified"]):
        genomeclass="s"

      if genomeclass:
        if taxonomyid in GOLD_class_by_taxonid:
          if genomeclass != GOLD_class_by_taxonid[taxonomyid]:
            if genomeclass == "n":
              # We already know this taxon is not non-symbiotic, not making it a non-symbiont again
              continue
        GOLD_class_by_taxonid[taxonomyid]=genomeclass

sys.stderr.write("GOLD: done (%i genomes).\n" % (genomescounter))
sys.stderr.flush()

sys.stderr.write("Non-symbionts: %i\n" % (list(GOLD_class_by_taxonid.values()).count("n")))
sys.stderr.write("Symbionts:     %i\n" % (list(GOLD_class_by_taxonid.values()).count("s")))




#BacMap
sys.stderr.write("Reading BacMap habitat classifications. \n")
BacMap_habitats={}
with open_un_compressed(BacMap_classification_infile) as infile:
  next(infile)
  for line in infile:
    (habitat,classification) = line.strip().split("\t")
    BacMap_habitats[habitat] = classification


sys.stderr.write("Reading BacMap genomes. \n")
BacMap_class_by_taxonid={}
with open_un_compressed(BacMap_infile_name) as infile:
  next(infile)
  for line in infile:
    (name,NcbiTaxonomyID,habitat,hostID,hostname) = line.strip().split(",")
    taxonomyid=int(NcbiTaxonomyID)
    if isbacterial(parent_by_nodeid, taxonomyid):
      genomeclass=""
      
      #BacMap organisms with host-associated habitats or eukaryotic-host organisms are classified as symbiotic
      if (habitat in BacMap_habitats and BacMap_habitats[habitat] == 's') or iseukaryotic(parent_by_nodeid, int(hostID)):
         genomeclass="s"
      
      #BacMap with other habitats and no host are classified as non-symbiotic, ambiguous habitats are excluded 
      elif (habitat in BacMap_habitats and BacMap_habitats[habitat] == 'n') and hostname == "NA":
        genomeclass="n"
      else:
        continue

      if genomeclass:
        if taxonomyid in BacMap_class_by_taxonid:
          if genomeclass != BacMap_class_by_taxonid[taxonomyid]:
            if genomeclass == "n":
              # We already know this taxon is not non-symbiotic, not making it a non-symbiont again
              continue
        BacMap_class_by_taxonid[taxonomyid]=genomeclass

sys.stderr.write("Non-symbionts: %i\n" % (list(BacMap_class_by_taxonid.values()).count("n")))
sys.stderr.write("Symbionts:     %i\n" % (list(BacMap_class_by_taxonid.values()).count("s")))




#NCBI
sys.stderr.write("Reading NCBI environment classifications. \n")
NCBI_environments={}
with open_un_compressed(NCBI_classification_infile) as infile:
  next(infile)
  for line in infile:
    (environment,classification) = line.strip().split("\t")
    NCBI_environments[environment] = classification


sys.stderr.write("Reading NCBI genomes.\n")
NCBI_class_by_taxonid={}
with open_un_compressed(NCBI_infile_name) as infile:
  next(infile)
  for line in infile:
    (acession,name,NcbiTaxonomyID,environment,hostname,hostID) = line.split("\t")
    taxonomyid=int(NcbiTaxonomyID)
    if isbacterial(parent_by_nodeid, taxonomyid):
      genomeclass=""
      #host associated environments classified as symbiotic
      if (environment in NCBI_environments and NCBI_environments[environment] == "s") or iseukaryotic(parent_by_nodeid, int(hostID)):
        genomeclass="s"
      #the other environments are classified as non-symbiotic
      elif (environment in NCBI_environments and NCBI_environments[environment] == 'n') and hostname == "NA":
        genomeclass="n"
      else:
        continue

      if genomeclass:
        if taxonomyid in NCBI_class_by_taxonid:
          if genomeclass != NCBI_class_by_taxonid[taxonomyid]:
            if genomeclass == "n":
              # We already know this taxon is not non-symbiotic, not making it a non-symbiont again
              continue
        NCBI_class_by_taxonid[taxonomyid]=genomeclass


sys.stderr.write("Non-symbionts: %i\n" % (list(NCBI_class_by_taxonid.values()).count("n")))
sys.stderr.write("Symbionts:     %i\n" % (list(NCBI_class_by_taxonid.values()).count("s")))




#creates classification by NCBItaxonomyID file per database
with open(classDB_outfile_name, "w") as classDB_outfile:
  for i in GOLD_class_by_taxonid:
    classDB_outfile.write("%s\t%i\t%s\n" % ("GOLD", i, GOLD_class_by_taxonid[i]))
  for i in BacMap_class_by_taxonid:
    classDB_outfile.write("%s\t%i\t%s\n" % ("BacMap", i, BacMap_class_by_taxonid[i]))
  for i in NCBI_class_by_taxonid:
    classDB_outfile.write("%s\t%i\t%s\n" % ("NCBI", i, NCBI_class_by_taxonid[i]))



#create consolidated classification output file from GOLD, BacMap anD NCBI
class_by_taxonid={}
for txid in GOLD_class_by_taxonid:
  class_by_taxonid[txid] = GOLD_class_by_taxonid[txid]

for txid in BacMap_class_by_taxonid:
  if txid not in class_by_taxonid: 
    class_by_taxonid[txid] = BacMap_class_by_taxonid[txid]
  elif txid in class_by_taxonid and class_by_taxonid[txid] == "n":
    class_by_taxonid[txid] = BacMap_class_by_taxonid[txid]

for txid in NCBI_class_by_taxonid:
  if txid not in class_by_taxonid:
    class_by_taxonid[txid] = NCBI_class_by_taxonid[txid]
  elif txid in class_by_taxonid and class_by_taxonid[txid] == "n":
    class_by_taxonid[txid] = NCBI_class_by_taxonid[txid]
  

with open(class_outfile_name, "w") as class_outfile:
  for txid in class_by_taxonid:
    class_outfile.write("%i\t%s\n" % (txid, class_by_taxonid[txid]))




#Reads file containing all Refseq Assemblies and classifies them as symbiotic or non-symbiotic
sys.stderr.write("Reading Refseq Assemblies.\n")
counter = 0
with open_un_compressed(RefeseqAssembliesInput) as infile:
  with open(genome_outfile_name, "w") as genome_outfile:
    for line in infile:
      if len(line.split("\t")) == 5:
        (accession, taxonomyID, species, strain, ftp) = line.split("\t")
        organism="%s %s" % (species, strain)
      elif len(line.split("\t")) == 4:
        (accession, taxonomyID, species, ftp) = line.split("\t")
        organism=species
      refseqtxid=int(taxonomyID)
      if refseqtxid in class_by_taxonid:
        genomeclass=class_by_taxonid[refseqtxid]
        if (parent_by_nodeid[refseqtxid] in class_by_taxonid):
          if class_by_taxonid[parent_by_nodeid[refseqtxid]] == 's' and class_by_taxonid[refseqtxid] == 'n':
            continue
        genome_outfile.write("%s\t%i\t%s\t%s\n" % (accession, refseqtxid, organism, class_by_taxonid[refseqtxid]))
        counter += 1

sys.stderr.write("done (%i genomes).\n" % (counter))
sys.stderr.flush()
