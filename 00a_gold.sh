#!/bin/bash
#
#SBATCH --job-name=gold
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=gold.out
#SBATCH --error=gold.err
#SBATCH --time=00-03:00:00

module load python3

#maximal allowed search size per page on GOLD is 500
CHUNKSIZE=500
FINISHED=0
PAGE=1

mkdir -p results
rm -f results/goldorganisms.txt.gz


#loop searches through GOLD organsism-database 500 entries each time, and downloads the page as a html file 
#the html file is then searched by the pyhton3 file bin/goldhtml2txt.py each time, which stores the environmental-metadata as a structured txt file


echo search GOLD organsism-database for environmental metadata
while [ $FINISHED -eq 0 ] ; do
  TMPFILE=$(mktemp)
  curl -s "https://gold.jgi.doe.gov/organisms?page=$PAGE&count=$CHUNKSIZE&Organism.NCBI+Taxonomy+ID_options=equals&count=25&setColumns=yes&Organism.NCBI+Taxonomy+ID=&Organism.Host+Taxonomy+ID=&Organism.Ecosystem=&Organism.Ecosystem+Category=&Organism.Habitat=" | bin/goldhtml2txt.py >$TMPFILE
  NLINES=$(cat $TMPFILE | wc -l)
  cat $TMPFILE >>results/goldorganisms.txt
  rm -f $TMPFILE
  if [ $NLINES -lt $CHUNKSIZE ] ; then
    FINISHED=1
  else
    echo -n '.'
    PAGE=$(( $PAGE + 1 ))
  fi
done

gzip results/goldorganisms.txt
echo bacmap gold-organism-environmental-metadata written into "results/goldorganisms.txt.gz"
