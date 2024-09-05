#!/bin/bash
#
#SBATCH --job-name=genome_classification
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=genome_classification.out
#SBATCH --error=genome_classification.err
#SBATCH --time=00-03:00:00


module load python3
module load edirect

#getting NCBI taxdump nodes from the lisc ncbi mirror
echo get NCBI taxdump nodes
tar -xvf /scratch/mirror/ncbi/current/taxonomy/taxdump.tar.gz files/nodes.dmp

#searching the ncbi taxonomy database for bacterial refseq assemblies and writing them into a csv-file
echo search NCBI for all bacteria refseq assemblies
esearch -db assembly -query '"latest refseq" AND "ok"[Taxonomy Check Status] AND "Bacteria[Organismgroup]"' | efetch -format docsum | xtract -pattern DocumentSummary -first element AssemblyAccession, Taxid, SpeciesName, Sub_value, FtpPath_RefSeq  > results/all_refseq_assemblies_per_txid.csv

#the bacterial assemblies from NCBI are classified using environmental metadata from GOLD, BacMap and NCBI-nuccore (environment, hosts)
echo classifying assemblies
python3 bin/class_by_taxon.py files/nodes.dmp results/goldorganisms.txt.gz results/BacMapHabitats.txt results/NCBI_metadata.txt results/BacMap_classification.txt results/NCBI_classification.txt results/all_refseq_assemblies_per_txid.csv results/classification results/classification_per_db results/genome_all


#rm -f files/nodes.dmp
