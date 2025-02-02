#!/usr/bin/env python3

import sys, os.path
from Bio import SeqIO

genomesfolder="genomes"
genomefilename="refseq_accessions_class.txt"
eukscorefilename="euk_score"

proteinfilename="protein"
proteinfilename_full="protein_full"
eukpredfilename="euk_prediction"

#reads in file containing all refseq accessions and their classification (n...non-symbiotic, s...symbiotc)
refseq_accessions=[]
with open(genomefilename) as infile:
  for line in infile:
    (refseq_accession, gclass) = line.strip().split("\t")
    if gclass == "s":
      refseq_accessions.append(refseq_accession)
sys.stderr.write('Got %i refseq_accessions of symbionts.\n' % len(refseq_accessions))

#reads in the list of eukaryotic domains, their descriptions and Wilcoxon signed-rank test threshold scores
eukdomains_by_refseq_accession={}
scores={}
with open(eukscorefilename) as infile:
  for line in infile:
    (domainid, refseq_accession, score) = line.strip().split("\t")
    scores[(refseq_accession, domainid)]=float(score)
    if refseq_accession not in eukdomains_by_refseq_accession:
      eukdomains_by_refseq_accession[refseq_accession]={}
    eukdomains_by_refseq_accession[refseq_accession][domainid]=1
sys.stderr.write('Got eukaryotic domains for %i refseq_accessions.\n' % len(eukdomains_by_refseq_accession.keys()))


#searches through all pfam-files and faa-files in the genomes-directory
protein_count=0
poutfile=open(proteinfilename, "w")
eoutfile=open(eukpredfilename, "w")
foutfile=open(proteinfilename_full, "w")
for refseq_accession in refseq_accessions:
  if refseq_accession not in eukdomains_by_refseq_accession:
    continue

  domainids_by_protein={}
  with open(os.path.join(genomesfolder, "%s.faa.pfam" % refseq_accession)) as infile:
    for line in infile:
      parts = line.strip().split("\t")
      proteinname = parts[0]
      domainid = parts[4]
      
      if domainid in eukdomains_by_refseq_accession[refseq_accession]:
        if proteinname not in domainids_by_protein:
          domainids_by_protein[proteinname]={}
        domainids_by_protein[proteinname][domainid]=1


  #gets protein-descriptions from faa-files
  annotation_by_protein={}
  with open(os.path.join(genomesfolder, "%s.faa" % refseq_accession)) as infile:
    for entry in SeqIO.parse(infile, "fasta"):
      if entry.id in domainids_by_protein:
        description = entry.description.split(" ",1)[1].split("[")[0].strip()
        if description.startswith("MULTISPECIES: "):
          description=description[14:]
        annotation_by_protein[entry.id]=(description)
  
  #writes refseq_accession, protein-accession, protein-description into "protein"-file
  for proteinname in domainids_by_protein.keys():
    description = ""
    if proteinname in annotation_by_protein:
      description = annotation_by_protein[proteinname]
    else:
      sys.stderr.write(proteinname)
    poutfile.write("%s\t%s\t%s\n" % (refseq_accession, proteinname, description))
    protein_count+=1

    if (refseq_accession,domainid) in scores.keys():
      score=scores[(refseq_accession,domainid)]
      foutfile.write("%s\t%s\t%s\t%s\t%1.5f\n" % (refseq_accession, proteinname, description, domainid, score))
    

    #writes refseq_accession, domainID, protein-accession into "euk_prediction"-file
    for domainid in domainids_by_protein[proteinname].keys():
      eoutfile.write("%s\t%s\t%s\n" % (refseq_accession, proteinname, domainid))
      



sys.stderr.write("\n")
sys.stderr.write('Got %i secreted proteins.\n' % protein_count)
poutfile.close()
eoutfile.close()
foutfile.close()