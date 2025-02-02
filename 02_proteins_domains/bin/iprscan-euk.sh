#!/bin/bash
#
#SBATCH --job-name=iprscan
#SBATCH --cpus-per-task=2
#SBATCH --mem=200
#SBATCH --output=log/iprscan-%A_%a.out
#SBATCH --error=log/iprscan-%A_%a.err
#SBATCH --partition=basic
#SBATCH --nice=10000
#SBATCH --mail-type=ALL

/localmirror/monthly/interpro/interproscan-*/interproscan.sh -cpu $SLURM_CPUS_PER_TASK -appl Pfam -f tsv -t p -i chunks/chunk.${SLURM_ARRAY_TASK_ID} -o $TMPDIR/chunk.${SLURM_ARRAY_TASK_ID}.pfam -T $TMPDIR && mv $TMPDIR/chunk.${SLURM_ARRAY_TASK_ID}.pfam chunks/