#!/usr/bin/env python3

import sys, gzip
from Bio import SeqIO

taxonomyfilename="nodes.dmp"
mergedfilename="merged.dmp"
speciesfilename="eggnog/eggnog4.species_list.txt"
proteinfilename="eggnog/eggnog4.proteins.all.fa.gz"

def get_class(taxonomyid, parent_by_node):
  gclass = ""

  node = taxonomyid
  while 1:
    parent = parent_by_node[node]

    if parent == 2759:
      gclass = "e"
      break

    if parent == node:
      break

    node = parent
  
  return gclass

parent_by_node={}
infile=open(taxonomyfilename)
for line in infile:
  parts = line.split("|")
  node_id = int(parts[0].strip())
  parent  = int(parts[1].strip())
  parent_by_node[node_id] = parent
infile.close

merged_by_node={}
infile=open(mergedfilename)
for line in infile:
  parts = line.split("|")
  node_id = int(parts[0].strip())
  parent  = int(parts[1].strip())
  parent_by_node[node_id] = parent
infile.close

euk_taxonomyids={}
infile=open(speciesfilename)
for line in infile:
  parts = line[:-1].split("\t")
  if line.startswith("#"):
    continue

  taxonomyid = int(parts[1])
  if taxonomyid in merged_by_node:
    taxonomyid = merged_by_node[taxonomyid]
  gclass = get_class(taxonomyid, parent_by_node)

  if gclass == "e":
    euk_taxonomyids[taxonomyid]=1
  
infile.close

with gzip.open(proteinfilename, "rt") as infile:
  for entry in SeqIO.parse(infile, "fasta"):
    taxonomyid = int(entry.id.split(".")[0])
    if taxonomyid in euk_taxonomyids:
      if entry.seq.endswith("*"):
        entry.seq = entry.seq[:-1]
      SeqIO.write(entry, sys.stdout, "fasta")
