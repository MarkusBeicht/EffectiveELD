#!/bin/bash
#
#SBATCH --job-name=1_assembly_per_txid
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=1_assembly_per_txid.out
#SBATCH --error=1_assembly_per_txid.err
#SBATCH --time=00-03:00:00

module load python3
module load edirect

#retrieves all bacterial refseq assemblies with their sort order
esearch -db assembly -query '"latest refseq" AND "ok"[Taxonomy Check Status] AND "Bacteria[Organismgroup]"' | efetch -format docsum | xtract -pattern DocumentSummary -def " " -first element AssemblyAccession, Taxid, SpeciesName, Sub_value, SortOrder > results/refseq_assemblies_SortOrder.txt


python3 bin/best_assembly_per_txid.py results/genome_all results/refseq_assemblies_SortOrder.txt results/genome

