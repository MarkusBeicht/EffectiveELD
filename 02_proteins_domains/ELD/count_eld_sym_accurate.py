#!/usr/bin/env python3

import sys, os.path, numpy, csv
from scipy.stats import wilcoxon

genomesfolder=sys.argv[1]
genomefilename=sys.argv[2]
eukdomainfilename=sys.argv[3]
nonsymdomaincountfilename=sys.argv[4]
precalcfilename=sys.argv[5]

#reads in file containing all refseq accessions and their classification (n...non-symbiotic, s...symbiotc)
refseq_accessions=[]
with open(genomefilename) as infile:
  for line in infile:
    (refseq_accession, gclass) = line.strip().split("\t")
    if gclass == "s":
      refseq_accessions.append(refseq_accession)
sys.stderr.write('Got %i refseq_accessions of symbionts.\n' % len(refseq_accessions))


#reads in the list of eukaryotic domains into dictionary
eukdomains={}
with open(eukdomainfilename) as infile:
  for line in infile:
    (name, interproid, domainid) = line.strip().split("\t")
    eukdomains[domainid]=[]
sys.stderr.write('Got %i eukaryotic domains.\n' % len(eukdomains.keys()))


#reads domaincounts of non-symbiotic bacteria into dictionary
counts=[]
with open(nonsymdomaincountfilename) as infile:
  for line in infile:
    (refseq_accession, domainid, count) = line.strip().split("\t")
    if domainid in eukdomains.keys():
      counts = eukdomains[domainid]
      counts.append(count)
      eukdomains[domainid]=counts
sys.stderr.write('Got %i eukaryotic domains.\n' % len(eukdomains.keys()))


#searches through all pfam files in the genomes-directory
poutfile=open(precalcfilename, "w")
for refseq_accession in refseq_accessions:
  proteins_by_domainid={}
  domaincounts=[refseq_accession]
  with open(os.path.join(genomesfolder, "%s.faa.pfam" % refseq_accession)) as infile:
    for line in infile:
      parts = line.strip().split("\t")
      proteinname = parts[0]
      domainid = parts[4]
      if domainid in eukdomains:
        if domainid not in proteins_by_domainid:
          proteins_by_domainid[domainid]={}
        proteins_by_domainid[domainid][proteinname]=1

  #the number of proteins in each pfam-file which belong to each eukaryotic domain is counted and then stored in the "domaincounts_sym"-csv-file
  for domainid in eukdomains.keys():
    genomecounts = eukdomains[domainid]
    count=0
    if domainid in proteins_by_domainid:
      count = len(proteins_by_domainid[domainid].keys())
      genomecounts.append(count)  
    domaincounts.append(count)

    #does a Wilcoxon signed rank test for each domain and writes the p-value into the output file
    #the differences between the count-value that is being tested and the genomecounts are calculated (ties are given +0.5, such that they function in support of the null-hypothesis)
    nonsym_mean = numpy.mean([int(x) for x in genomecounts])
    differences = numpy.array([int(x) - count  + 0.5 if int(x) == count else int(x) - count for x in genomecounts])
    W1,p = wilcoxon(differences, alternative='less')
    
    #rank-biserial correlation
    #the positive ranks (W0) are calculated using the overall ranks (W) and the negative ranks (W1) from the wilcoxon test
    n = len(differences)
    W = (n * (n + 1)) / 2
    W0 = W - W1
    rb = abs(W1 - W0) / W
    poutfile.write("%s\t%s\t%1.5f\t%1.5f\t%1.5f\n" % (domainid, refseq_accession, p, rb, nonsym_mean))
    sys.stdout.flush()

sys.stderr.write("\n")