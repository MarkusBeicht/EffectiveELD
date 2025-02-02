#!/bin/bash
# 
#SBATCH --job-name=ELD
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=ELD.out
#SBATCH --error=ELD.err
#SBATCH --time=00-01:00:00

module load python3/

#calculates p-values for required domains of the domains for the validation of the scoring method
python3 ELD/count_eld_sym_accurate.py effectors.org/genomes_ELD ELD/refseq_accessions_class_ELD.txt ELD/euk_domain_ELD domaincounts_nonsym.txt ELD/euk_score_ELD

#calculates ROC and AUC
python3 ELD/ROC.py ELD/ELD_data ELD/euk_score_ELD ELD/ROC.png 0 0

#calculates ROC and AUC, requiring a mean domaincount of at most 0.5
python3 ELD/ROC.py ELD/ELD_data ELD/euk_score_ELD ELD/ROC_mean05.png 0.5 0

#calculates ROC and AUC with varying effect sizes 
python3 ELD/ROC.py ELD/ELD_data ELD/euk_score_ELD ELD/ROC_rb01.png 0 0.1
python3 ELD/ROC.py ELD/ELD_data ELD/euk_score_ELD ELD/ROC_rb03.png 0 0.3
python3 ELD/ROC.py ELD/ELD_data ELD/euk_score_ELD ELD/ROC_rb05.png 0 0.5
python3 ELD/ROC.py ELD/ELD_data ELD/euk_score_ELD ELD/ROC_rb07.png 0 0.7

