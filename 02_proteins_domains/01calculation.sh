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

echo Determining eukaryotic domains from annotations of eukaryotic EggNOG genomes
cat chunks/chunk.*.pfam | cut -f1,5,6 | sed 's/\.[^\t]*\t/\t/' | uniq | sort | uniq | cut -f2,3 | sort | uniq -c | tr -s " " | sed 's/^ //' | awk '$1>2'| cut -f2- -d " " >euk_domain_list 

echo Counting euk domains in non-symbionts and calculating thresholds for Wilcoxon signed-rank test
bin/count_eld_nonsym.py

echo Counting euk domains in symbionts and scoring them against the Wilcoxon signed-rank test thresholds of the non-symbionts 
bin/count_eld_sym.py

echo Extracting proteins containing an ELD from the faa-files and storing their domain association and annotation
bin/get_eld_proteins.py


echo Removing files that are no longer needed
#rm -rf effectors.org
echo calculation completed
