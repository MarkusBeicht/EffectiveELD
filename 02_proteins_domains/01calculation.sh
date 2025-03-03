#!/bin/bash
# 
#SBATCH --job-name=calculation
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=calculation.out
#SBATCH --error=calculation.err
#SBATCH --time=00-12:00:00

module load python3/

echo Extracting the RefSeq accessions and their classifications from genome file
cut -f1,4 genome > refseq_accessions_class.txt

echo Determines eukaryotic domains from annotations of eukaryotic EggNOG genomes
cat chunks/chunk.*.pfam | cut -f1,5,6 | sed 's/\.[^\t]*\t/\t/' | uniq | sort | uniq | cut -f2,3 | sort | uniq -c | tr -s " " | sed 's/^ //' | awk '$1>2'| cut -f2- -d " " >euk_domain_list 

echo Counting euk domains in non-symbionts and calculating thresholds for Wilcoxon signed-rank test
bin/count_eld_nonsym.py

echo Counting euk domains in symbionts and scoring them against the Wilcoxon signed-rank test thresholds of the non-symbionts 
bin/count_eld_sym.py

echo Extracting proteins containing an ELD from the faa-files and storing their domain association and annotation
bin/get_eld_proteins.py


echo Removing files that are no longer needed
rm -rf effectors.org
echo calculation completed


echo Extracting ELDs and ED into subsets of Pfam HMM file - ELD needed in protein mode and ED needed in genome mode
rm -f ../data/pfamids_euk?.hmm*
cut -f1 euk_score | sort | uniq >pfamids_euk4.txt
cut -f1 euk_domain >pfamids_euk0.txt
bin/extract_hmm3_models.py /lisc/scratch/mirror/interpro/interproscan-*/data/pfam/*/pfam_a.hmm pfamids_euk4.txt pfamids_euk4.hmm
bin/extract_hmm3_models.py /lisc/scratch/mirror/interpro/interproscan-*/data/pfam/*/pfam_a.hmm pfamids_euk0.txt pfamids_euk0.hmm
/lisc/scratch/mirror/interpro/interproscan-*/bin/hmmer/hmmer3/3.3*/hmmpress -f pfamids_euk4.hmm
/lisc/scratch/mirror/interpro/interproscan-*/bin/hmmer/hmmer3/3.3*/hmmpress -f pfamids_euk0.hmm


echo Removing files that are no longer needed
rm -f pfamids_euk4.txt pfamids_euk0.txt refseq_accessions_class.txt euk_domains.txt log chunks
rm -rf effectors.org
echo calculation completed