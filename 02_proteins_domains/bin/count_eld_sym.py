#!/usr/bin/env python3

import sys, os.path, csv

genomesfolder="genomes"
genomefilename="refseq_accessions_class.txt"
eukdomainfilename="euk_domain"
domaincountfilenamelist="domaincounts_sym.txt"
domaincountfilenametable="domaincounts_sym.csv"
precalcfilename="euk_score"

#reads in file containing all refseq accessions and their classification (n...non-symbiotic, s...symbiotc)
refseq_accessions=[]
with open(genomefilename) as infile:
  for line in infile:
    (refseq_accession, gclass) = line.strip().split("\t")
    if gclass == "s":
      refseq_accessions.append(refseq_accession)
sys.stderr.write('Got %i refseq_accessions of symbionts.\n' % len(refseq_accessions))

#reads in the list of eukaryotic domains, their descriptions and MWU-test-thresholds
eukdomains={}
thresholds={}
with open(eukdomainfilename) as infile:
  for line in infile:
    (domainid, desc, t, t01, t03, t05, t07) = line.strip().split("\t")
    eukdomains[domainid]=(desc, [])
    thresholds[domainid]=(int(t), int(t01), int(t03), int(t05), int(t07))
sys.stderr.write('Got %i eukaryotic domains.\n' % len(eukdomains.keys()))

#writes eukaryotic domains as header row into the "domaincounts_sym"-csv-file for storage of all domaincounts of genomes classified as symbiotc
with open(domaincountfilenametable, "w") as toutfile:
  write = csv.writer(toutfile)
  header = ["refseq_accession"]
  for i in eukdomains.keys():
    header.append(i)
  write.writerow(header)



#searches through all pfam files in the genomes-directory
poutfile=open(precalcfilename, "w")
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

  #the number of proteins in each pfam-file which belong to each eukaryotic domain is counted and then stored in the "domaincounts_sym"-csv-file
  for domainid in eukdomains.keys():
    (description, genomecounts) = eukdomains[domainid]
    count=0
    if domainid in proteins_by_domainid:
      count = len(proteins_by_domainid[domainid].keys())
      genomecounts.append(count)  
    domaincounts.append(count)
    

    #writes each significantly enriched domain with their highest reached effect-size threshold (r >= 0.7, r >= 0.5, r >= 0.3, r >= 0.1, r < 0.1) into "euk_score"-file
    (t, t01, t03, t05, t07) = thresholds[domainid]
    if count >= t07:
      poutfile.write("%s\t%s\t%1.5f\n" % (domainid, refseq_accession, 0.7))
    elif count >= t05:
      poutfile.write("%s\t%s\t%1.5f\n" % (domainid, refseq_accession, 0.5))
    elif count >= t03:
      poutfile.write("%s\t%s\t%1.5f\n" % (domainid, refseq_accession, 0.3))
    elif count >= t01:
      poutfile.write("%s\t%s\t%1.5f\n" % (domainid, refseq_accession, 0.1))
    elif count >= t:
      poutfile.write("%s\t%s\t%1.5f\n" % (domainid, refseq_accession, 0.0))

    #writes domaincount of each domain for each genome
    coutfile.write("%s\t%s\t%i\n" % (refseq_accession, domainid, count))
    sys.stdout.flush()


  #writes list of domaincounts for each genome as a table
  with open(domaincountfilenametable, "a") as toutfile:
    write = csv.writer(toutfile)
    write.writerow(domaincounts)

sys.stderr.write("\n")