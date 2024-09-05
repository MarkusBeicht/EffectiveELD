#!/bin/bash
#
#SBATCH --job-name=ncbi_nuccore
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=ncbi_nuccore.out
#SBATCH --error=ncbi_nuccore.err
#SBATCH --time=00-03:00:00


#The NCBI-nuccore datasbase is searched for environmental metadata

module load edirect
mkdir -p files
mkdir -p results


#The accession numbers of the master records with environment metadata (MIGS-6) are retrieved from the NCBI-nuccore datasbase. 
echo search NCBI-nuccore datasbase for environmental metadata
esearch -db nuccore -query '(((environment)) AND ("Bacteria"[Organism]) AND ( "tsa master"[Properties] OR "wgs master"[Properties] ))' | efetch -format docsum | xtract -pattern DocumentSummary -element Caption > files/NCBI_environment_records.txt
echo List of NCBI nuccore master records with environment metadata MIGS-6 written into "files/NCBI_environment_records.txt"


name=$(mktemp)
NCBITaxonomyID=$(mktemp)
environment=$(mktemp)
hostname=$(mktemp)
hostID=$(mktemp)

echo accession$'\t'name$'\t'NCBITaxonomyID$'\t'environment$'\t'hostname$'\t'hostID > files/NCBI_environment_data.txt

while read -u 9 p; do

#retrieving organism's name from DocumentSummary file
echo  "$p"$'\t'| tr -d '\n' >> files/NCBI_environment_data.txt
	name="$(efetch -format docsum -db nuccore -id "$p" | xtract -pattern DocumentSummary -element Organism, Strain | sed 's/\t/ /g')"
	

	#retrieving ASN.1 file of each NCBI nucleotide master entry and searching it for environmental metadata (environment & host) 
	efetch -db nuccore -format asn1 -id "$p" > MIGS.asn1	

	if grep -q "db \"taxon\","  MIGS.asn1; then
		NCBITaxonomyID="$(cat MIGS.asn1 | grep "db \"taxon\"," -A1 -m 1 | tail -n 1 | cut -f2 -d "d" | cut -f2 -d " " | sed -e 's/^/\t/')"
	else
		NCBITaxonomyID="NA"
	fi


	if grep -q "label str \"environment\"," MIGS.asn1; then
		environment="$(cat MIGS.asn1 | grep "label str \"environment\"," -A1 -m 1 | tail -n 1 | cut -f2 -d "\"" | cut -f1 -d "\"" | sed -r 's/[,]+/;/g' | sed -e 's/^/\t/' | tr -d '\n')"
	else
		environment="NA"
	fi


	#retrieving hostname from asn1 file
	if grep -q "subtype nat-host," MIGS.asn1; then
		hostname=""
		hostname="$(cat MIGS.asn1 | grep "subtype nat-host," -A1 -m 1 | tail -n 1 | cut -f2 -d "\"" | cut -f1 -d "\"" | sed -r 's/[,]+/;/g' | sed -e 's/^/\t/')" >> files/NCBI_environment_data.txt


		#search for hostID in NCBI taxonomy database using hostname
		if [[ "$hostname" != "" ]]; then
			hostID="$(esearch -db taxonomy -query "$hostname" | efetch -stop 1 -format docsum | xtract -pattern DocumentSummary -element Id)"
		fi

		if [[ "$hostID" == "" ]]; then
			hostID=0
		fi

	else 
		hostname="NA"
		hostID=0
	fi

	echo "$p"$'\t'$name$'\t'$NCBITaxonomyID$'\t'$environment$'\t'$hostname$'\t'$hostID >> files/NCBI_environment_data.txt

done 9<files/NCBI_environment_records.txt
rm -f MIGS.asn1




echo NCBI-nuccore-environmental-metadata written into "files/NCBI_environment_data.txt"
rm -f files/NCBI_environment_records.txt


#searching NCBI nuccore database for host entries
echo search NCBI nuccore database for host entries
esearch -db nuccore -query '(((host)) AND ( "tsa master"[Properties] OR "wgs master"[Properties] ))' | efetch -format docsum | xtract -pattern DocumentSummary -def " " -element Caption, Organism, Strain, TaxId, SubType, SubName > files/NCBI_host_records.txt
NCBI_host.py extracts the hostnames from NCBI_host_records.txt and searches them against the NCBI-taxonomy database to retrieve the NCBItaxonomyID of the host
python3 bin/NCBI_host.py files/NCBI_host_records.txt files/NCBI_host_data.txt

#rm -f files/NCBI_host_records.txt


#the environment data is written into NCBI_metadata.txt, and the host_data entries that are not yet included in the environment_data file are added
echo consolidating the environment and host metadata from NCBI in NCBI_metadata.txt
cat files/NCBI_environment_data.txt > results/NCBI_metadata.txt && awk 'NR==FNR{a[$1]; next} {if (!($1 in a)) print}' files/NCBI_environment_data.txt files/NCBI_host_data.txt >> results/NCBI_metadata.txt 
rm -f files/NCBI_environment_data.txt
rm -f files/NCBI_host_data.txt


