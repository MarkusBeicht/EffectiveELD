#!/usr/bin/env python3

import sys

hmminfilename=sys.argv[1]
namesfilename=sys.argv[2]
hmmoutfilename=sys.argv[3]

accessions={}
infile=open(namesfilename)
for line in infile:
  accessions[line.strip()]=1
infile.close()

infile=open(hmminfilename)
outfile=open(hmmoutfilename,"w+")
HMM=[]
for line in infile:
  HMM.append(line)
  if line.startswith("ACC   "):
    name = line.strip()[6:].split(".")[0]
  if line.startswith("//"):
    if name in accessions:
      for hmmline in HMM:
        outfile.write(hmmline)
    HMM=[]
infile.close()
outfile.close()
