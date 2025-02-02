#!/usr/bin/env python3

"""split_seqfile.py Split larger sequence files into smaller files
usage: split_seqfile.py infilename outfileprefix new_size
"""

import string, sys
from Bio import SeqIO

outfileprefix = sys.argv[1]
new_size = int(sys.argv[2])
zerofill100=0
if len(sys.argv)>3:
  zerofill100=int(sys.argv[3])

n_seq = 0
if not zerofill100:
  newfileidx = 1
  outfile = open("%s.%i" % (outfileprefix, newfileidx), "w+")
else:
  newfileidx = 0
  outfile = open("%s.%02i" % (outfileprefix, newfileidx), "w+")

for entry in SeqIO.parse(sys.stdin, "fasta"):
  n_seq += 1
  if n_seq > new_size:
    newfileidx += 1
    n_seq -= new_size
    outfile.close()
    if not zerofill100:
      outfile = open("%s.%i" % (outfileprefix, newfileidx), "w+")
    else:
      outfile = open("%s.%02i" % (outfileprefix, newfileidx), "w+")
   
  SeqIO.write(entry, outfile, "fasta")

outfile.close()

if not zerofill100:
  print(newfileidx)
