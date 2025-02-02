#!/usr/bin/env python3

import sys, os.path, numpy, csv
from scipy.stats import wilcoxon

genomesfolder="genomes"
genomefilename="refseq_accessions_class.txt"
eukdomainfilename="euk_domain_list"
domaincountfilenamelist="domaincounts_nonsym.txt"
domaincountfilenametable="domaincounts_nonsym.csv"
domainthresholdfilename="euk_domain"

#reads in file containing all refseq accessions and their classification (n...non-symbiotic, s...symbiotc)
refseq_accessions=[]
with open(genomefilename) as infile:
  for line in infile:
    (refseq_accession, gclass) = line.strip().split("\t")
    if gclass == "n":
      refseq_accessions.append(refseq_accession)
sys.stderr.write('Got %i refseq_accessions of non-symbionts.\n' % len(refseq_accessions))


#reads in the list of eukaryotic domains and their descriptions
eukdomains={}
with open(eukdomainfilename) as infile:
  for line in infile:
    (domainid, desc) = line.strip().split("\t")
    eukdomains[domainid]=(desc, [])
sys.stderr.write('Got %i eukaryotic domains.\n' % len(eukdomains.keys()))


#writes eukaryotic domains as header row into the "domaincounts_nonsym"-csv-file for storage of all domaincounts in the non-symbiotc training dataset
with open(domaincountfilenametable, "w") as toutfile:
  write = csv.writer(toutfile)
  header = ["refseq_accession"]
  for i in eukdomains.keys():
    header.append(i)
  write.writerow(header)


#searches through all pfam files in the genomes-directory
coutfile=open(domaincountfilenamelist, "w")
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

  #the number of proteins in each pfam-file which belong to each eukaryotic domain is counted and then stored in the "domaincounts_nonsym"-csv-file
  for domainid in eukdomains.keys():
    (description, genomecounts) = eukdomains[domainid]
    count=0
    if domainid in proteins_by_domainid:
      count = len(proteins_by_domainid[domainid].keys())
    genomecounts.append(count)
    domaincounts.append(count)

    #writes domaincount of each domain for each genome
    coutfile.write("%s\t%s\t%i\n" % (refseq_accession, domainid, count))
    sys.stdout.flush()

  
  #writes list of domaincounts for each genome as a table
  with open(domaincountfilenametable, "a") as toutfile:
    write = csv.writer(toutfile)
    write.writerow(domaincounts)

sys.stderr.write("\n")



#computes the Wilcoxon signed-rank-test at the p-value 0.05 and calculates the threshold values at which a certain effect size is reached (0.1, 0.3, 0.5)
#the threshold values are stored the in the "euk_domain_thresholds"-file

doutfile=open(domainthresholdfilename, "w")
for domainid in eukdomains.keys():
  (description, genomecounts) = eukdomains[domainid]

  #no effect size limitation
  threshold, p = 0, 1
  while p > 0.05:
    threshold += 1
    #the differences between the thershold value that is being tested and the genomecounts are calculated (ties are given +0.5 in support of the null-hypothesis)
    differences = numpy.array([x - threshold + 0.5 if x == threshold else x - threshold for x in genomecounts])
    W,p = wilcoxon(differences, alternative='less')
  

  #small effect size: rb > 0.1
  threshold_01, p, rb = 0, 1, 0
  while p > 0.05 or rb < 0.1:
    threshold_01 += 1
    #the differences between the thershold value that is being tested and the genomecounts are calculated (ties are set = +0.5, in support of the null-hypothesis)
    differences = numpy.array([x - threshold_01 + 0.5 if x == threshold_01 else x - threshold_01 for x in genomecounts])
    W1,p = wilcoxon(differences, alternative='less')
    
    #rank-biserial correlation
    #the positive ranks (W0) are calculated using the overall ranks (W) and the negative ranks (W1) from the wilcoxon test
    n = len(differences)
    W = (n * (n + 1)) / 2
    W0 = W - W1
    rb = abs(W1 - W0) / W

  #medium effect size: rb > 0.3
  threshold_03, p, rb = 0, 1, 0
  while p > 0.05 or rb < 0.3:
    threshold_03 += 1
    #the differences between the thershold value that is being tested and the genomecounts are calculated (ties are set = +0.5, in support of the null-hypothesis)
    differences = numpy.array([x - threshold_03 + 0.5 if x == threshold_03 else x - threshold_03 for x in genomecounts])
    W1,p = wilcoxon(differences, alternative='less')
    
    #rank-biserial correlation
    #the positive ranks (W0) are calculated using the overall ranks (W) and the negative ranks (W1) from the wilcoxon test
    n = len(differences)
    W = (n * (n + 1)) / 2
    W0 = W - W1
    rb = abs(W1 - W0) / W


  #large effect size: rb > 0.5
  threshold_05, p, rb = 0, 1, 0
  while p > 0.05 or rb < 0.5:
    threshold_05 += 1
    #the differences between the thershold value that is being tested and the genomecounts are calculated (ties are set = +0.5, in support of the null-hypothesis)
    differences = numpy.array([x - threshold_05 + 0.5 if x == threshold_05 else x - threshold_05 for x in genomecounts])
    W1,p = wilcoxon(differences, alternative='less')
    
    #rank-biserial correlation
    #the positive ranks (W0) are calculated using the overall ranks (W) and the negative ranks (W1) from the wilcoxon test
    n = len(differences)
    W = (n * (n + 1)) / 2
    W0 = W - W1
    rb = abs(W1 - W0) / W


  #very large effect size: rb > 0.7
  threshold_07, p, rb = 0, 1, 0
  while p > 0.05 or rb < 0.7:
    threshold_07 += 1
    #the differences between the thershold value that is being tested and the genomecounts are calculated (ties are set = +0.5, in support of the null-hypothesis)
    differences = numpy.array([x - threshold_07 + 0.5 if x == threshold_07 else x - threshold_07 for x in genomecounts])
    W1,p = wilcoxon(differences, alternative='less')
    
    #rank-biserial correlation
    #the positive ranks (W0) are calculated using the overall ranks (W) and the negative ranks (W1) from the wilcoxon test
    n = len(differences)
    W = (n * (n + 1)) / 2
    W0 = W - W1
    rb = abs(W1 - W0) / W

  doutfile.write("%s\t%s\t%i\t%i\t%i\t%i\t%i\n" % (domainid, description, threshold, threshold_01, threshold_03, threshold_05, threshold_07))
  sys.stdout.flush()
sys.stderr.write("\n")
