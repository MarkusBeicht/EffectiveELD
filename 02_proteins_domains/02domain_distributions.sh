#!/bin/bash
# 
#SBATCH --job-name=domain_distributions
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=calculation.out
#SBATCH --error=calculation.err
#SBATCH --time=00-04:00:00

module load python3/

python3 domain_distributions/domain_analysis.py domaincounts_nonsym.txt domain_distributions/domain_analysis.txt domain_distributions/plots/
python3 domain_distributions/goodnessoffit_control.py domain_distributions/goodnessoffit_control.txt
