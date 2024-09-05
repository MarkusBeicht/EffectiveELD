#!/bin/bash
#
#SBATCH --job-name=bacmap
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=bacmap.out
#SBATCH --error=bacmap.err
#SBATCH --time=00-03:00:00

module load edirect
mkdir -p results


TMPFILE=$(mktemp)
ncbilink=$(mktemp)
NcbiTaxonomyID=$(mktemp)
name=$(mktemp)
habitat=$(mktemp)
hostname=$(mktemp)
hostid=$(mktemp)

name,NcbiTaxonomyID,habitat,hostid,hostname > results/BacMapHabitats.txt

#opens the website of each organism in the bacmap database, download the html code, and searches through the html file
#organisms which do not have a habitat entry on their BacMap page are filtered out

echo searches through all organism-pages in the bacmap-database
for i in {1..1768} 
do
	#retrieves the organims bacmap page html-code
	curl -s "http://bacmap.wishartlab.com/organisms/$i" > $TMPFILE

	#extracts the name, habitat, hostname
	name="$(grep -i "<th>Names</th>" $TMPFILE -A1 | tail -n1 | cut -f2 -d ">" | cut -f1 -d "<" | sed -e $'s/,/\\\n/g')"
	habitat="$(grep -i "<th>Habitat</th>" $TMPFILE -A1 | tail -n1 | cut -f2 -d ">" | cut -f1 -d "<" | sed -e $'s/,/\\\n/g')"
	hostname="$(grep -i "<th>Host name</th>" $TMPFILE -A1 | tail -n1 | cut -f2 -d ">" | cut -f1 -d "<" | sed 's/,/ /g')"	

	#the host names are then searched against the NCBI nucleotide database to assign an NCBI taxonomy ID to the host, if possible
	if [[ "$hostname" != "NA" ]]; then
		hostid="$(esearch -db taxonomy -query "$hostname" | efetch -stop 1 -format docsum | xtract -pattern DocumentSummary -element Id)"
	fi

	if [[ "$hostid" == "" ]]; then
		hostid=0
	fi


	#the ncbilink is used to retrieve the NCBItaxonomyID from the NCBI nuccore database
	ncbilink="$(grep -m 1 "nuccore" $TMPFILE | cut -f2 -d"=" | cut -f1 -d ">" | sed -e 's/^.//' -e 's/.$//')"

	if [[ "$ncbilink" = *"http:"* ]]; then
		ncbilink=${ncbilink/http/https}
	fi

	curl -s $ncbilink > ncbi
	NcbiTaxonomyID="$(grep -m 1 "ORGANISM=" ncbi | sed 's/^.*ORGANISM=//' | cut -f1 -d "&" )"
	
	
	if [[ "$NcbiTaxonomyID" == "" ]]; then
		NcbiTaxonomyID=0
	fi
	

	#if [[ "$habitat" != "NA" ]]; then
	echo $name,$NcbiTaxonomyID,$habitat,$hostid,$hostname >> results/BacMapHabitats.txt
	#fi
done

echo bacmap organism-environmental-metadata written into "results/BacMapHabitats.txt"
rm -f $TMPFILE
rm -f $ncbilink
rm -f ncbi
rm -f $name
rm -f $NcbiTaxonomyID
rm -f $habitat