#!/bin/bash
#SBATCH --output=submit.out
#SBATCH --error=submit.err
#SBATCH --license=interpro

mkdir -p effectors.org/chunks effectors.org/log effectors.org/genomes
ln -s -f effectors.org/chunks
ln -s -f effectors.org/log
ln -s -f effectors.org/genomes


CHUNKSIZE_PROK=10
cut -f1 genome > effectors.org/submit.txt
genomes=$(wc -l < effectors.org/submit.txt)
echo $((${genomes} / ${CHUNKSIZE_PROK})) array jobs
sbatch -a 0-$((${genomes} / ${CHUNKSIZE_PROK})) bin/iprscan-prok.sh ${CHUNKSIZE_PROK}


CHUNKSIZE_EUK=10000
NCHUNKS=$(bin/dump_eggnog_euk_proteins.py | tr '*' 'X' | bin/split_seqfile.py chunks/chunk ${CHUNKSIZE_EUK})
sbatch -a 1-$NCHUNKS bin/iprscan-euk.sh

