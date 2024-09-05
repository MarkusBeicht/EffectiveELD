#!/usr/bin/env python3

import sys, os

genome_infile_name=sys.argv[1]
sortorder_infile_name=sys.argv[2]
genome_outfile_name=sys.argv[3]


#reads the refseq_sortorders file and stores the sortorder for each refseq assembly
refseq_sortorders={}
with open(sortorder_infile_name, "rt") as sortorder_infile:
  for line in sortorder_infile:
    (refseq_accession, NcbiTaxonomyID, name, strain, sortorder) = line.strip().split("\t")
    taxonomyid=int(NcbiTaxonomyID)
    refseq_sortorders[refseq_accession]=sortorder



#reads all entries in the genome file 
txid_min_sortorder={}
best_assembly_per_txid={}
with open(genome_infile_name, "rt") as infile:
  for line in infile:
    (refseq_accession, NcbiTaxonomyID, name, classification) = line.strip().split("\t")
    taxonomyid=int(NcbiTaxonomyID)

    #finds the best refseq_assembly per NcbiTaxonomyID (lowest sortorder)
    if classification == "s":
        
      #sets default sortorder at max value
      genome_sortorder="9C6XFFFF9999999999999999999999"
      #updates sortorder of refseq_assembly if it is in refseq_sortorders file 
      if refseq_accession in refseq_sortorders.keys():
        genome_sortorder=refseq_sortorders[refseq_accession]
        
      #checks if the 
      if taxonomyid not in txid_min_sortorder.keys():
        txid_min_sortorder[taxonomyid]=genome_sortorder
      elif genome_sortorder < txid_min_sortorder[taxonomyid]:
        txid_min_sortorder[taxonomyid]=genome_sortorder
      best_assembly_per_txid[taxonomyid]=refseq_accession



counter = 0
with open(genome_infile_name, "rt") as infile:
  with open(genome_outfile_name, "w") as outfile:
    for line in infile:
      (refseq_accession, NcbiTaxonomyID, name, classification) = line.strip().split("\t")
      taxonomyid=int(NcbiTaxonomyID)
      
      #writes all non-symbiotic assemblies to outfile
      if classification == "n":
        outfile.write("%s\t%i\t%s\t%s\n" % (refseq_accession, taxonomyid, name, classification))
        counter += 1
      
      #writes all the best symbiotic assembly per NcbiTaxonomyID to outfile
      elif refseq_accession in best_assembly_per_txid.values():
        outfile.write("%s\t%i\t%s\t%s\n" % (refseq_accession, taxonomyid, name, classification))
        counter += 1

sys.stderr.write("done (%i genomes).\n" % (counter))
sys.stderr.flush()
