#!/bin/bash
#
#SBATCH --job-name=iprscan
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --output=log/iprscan-%A_%a.out
#SBATCH --error=log/iprscan-%A_%a.err
#SBATCH --partition=basic
#SBATCH --nice=1000
#SBATCH --mail-type=ALL
#SBATCH --license=interpro
#SBATCH --time=00-01:00:00


CHUNKSIZE=$1
REFSEQACCESSIONFILE=/lisc/scratch/cube/beicht/ELD/02_proteins_domains/effectors.org/submit.txt
STARTINDEX=$(( (${SLURM_ARRAY_TASK_ID} * $CHUNKSIZE) +1))
ENDINDEX=$(( $STARTINDEX+$CHUNKSIZE ))

#iterates through all refseq-accessions in the "submit.txt"-file 
for refseq_accession in $(sed -n "${STARTINDEX},${ENDINDEX}p" $REFSEQACCESSIONFILE); do
  
  #downloads the faa-protein-assembly, if it is not already in the directory
  if [ ! -s genomes/$refseq_accession.faa.pfam ] ; then
    g1=$(grep -P "$(echo $refseq_accession | sed 's/[^ftp]*\(ftp.*\)/\1/')" /lisc/scratch/cube/beicht/ELD/01_classification/results/all_refseq_assemblies_per_txid.csv | sed 's/^.*ftp:/ftp:/')
    g2=$(grep -P "$(echo $refseq_accession | sed 's/[^ftp]*\(ftp.*\)/\1/')" /lisc/scratch/cube/beicht/ELD/01_classification/results/all_refseq_assemblies_per_txid.csv | sed 's/^.*ftp:/ftp:/' | awk -F"/" '{print $NF}')
    wget $g1'/'$g2'_protein.faa.gz' -O genomes/$refseq_accession.faa.gz

    #unzips all ".gz"-zipped files in the directory
    if [ $(ls genomes/ | grep -c .faa.gz) -gt 0 ] ; then 
      gunzip genomes/$refseq_accession.faa.gz
      
      #excludes genomes which have less/equal than 100 proteins in their assembly (check for incompleteness)
      if [ $(grep -c '^>' genomes/$refseq_accession.faa) -le 100 ] ; then
        #echo genomes/$refseq_accession.faa >> failed.txt
        rm -f $refseq_accession.faa
      
      #performs an interpro-scan with the faa-file for the Pfam-application with precalculations enabled
      else
        echo $g1'/'$g2'_protein.faa.gz'
	/lisc/scratch/mirror/interpro/interproscan-*/interproscan.sh -appl Pfam -f TSV -i genomes/$refseq_accession.faa -o $TMPDIR/$refseq_accession.faa.pfam -T $TMPDIR && mv $TMPDIR/$refseq_accession.faa.pfam genomes/
      fi
    fi
  fi
done

rm -f $REFSEQACCESSIONFILE
rm -f genomes/*.faa.gz

