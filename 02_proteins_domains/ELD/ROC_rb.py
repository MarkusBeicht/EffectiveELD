#!/usr/bin/env python3

import sys, os.path, numpy
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, confusion_matrix
import copy

ELDfilename=sys.argv[1]
predictionfilename=sys.argv[2]
outfilename=sys.argv[3]


#loading ELD data from the literature 
ELDcounts={}
InterproPfamMapping={}
min_pvalue={}
with open(ELDfilename) as infile:
  for line_number, line in enumerate(infile, start=1):
    (InterproName, InterproID, PfamIDs, *data) = line.strip().split("\t")
    if line_number == 1:
      organisms = [organism.strip() for organism in data]
    elif line_number == 2:
      refseq_accessions = [refseq_accession for refseq_accession in data]
    else:
      ELDs = [int(ELD) for ELD in data]
      for i, ELD in enumerate(ELDs):
        ELDcounts[InterproID, refseq_accessions[i]]=min(ELD, 1)
        min_pvalue[InterproID, refseq_accessions[i]]=1
      InterproPfamMapping[InterproID]=[PfamID.strip() for PfamID in PfamIDs.split(',')]

saved_state = copy.deepcopy(min_pvalue)

#loading Wilcoxon signed-rank test p-values from euk_score_ELD file
effect_sizes = [0.7, 0.5, 0.3, 0.1]
for effect_size in effect_sizes:
  min_pvalue = saved_state
  with open(predictionfilename) as infile:
    for line in infile:
      (PfamID, refseq_accession, pvalue, rb, nonsym_mean) = line.strip().split("\t")
      p=float(pvalue)
      nonsym_mean = float(nonsym_mean)
      rb = float(rb)
      if rb < effect_size:
        p = 1
      InterproID = next((InterproID for InterproID, PfamIDs in InterproPfamMapping.items() if PfamID in PfamIDs))
      if p < min_pvalue[(InterproID, refseq_accession)]:
        min_pvalue[(InterproID, refseq_accession)]=p
  val = numpy.array(list(ELDcounts.values()))
  prob = 1 - numpy.array(list(min_pvalue.values()))
  fpr, tpr, thresholds = roc_curve(val, prob, pos_label=1)
  roc_auc = auc(fpr, tpr)
  print('AUC: ', roc_auc)
  plt.plot(fpr, tpr, label=f'ROC, r = {effect_size:.1f} (AUC = {roc_auc:.4f})')

plt.plot([0, 1], [0, 1], 'r--', label='Random Guess')
plt.title('ROC Curve of the EffectiveELD scoring method\n(requiring minimal effect sizes)')
plt.xlabel('False Positive Rate (1 - Specificity)')
plt.ylabel('True Positive Rate (Sensitivity)')
plt.text(2, 7, f"AUC: {roc_auc:.2f}", fontsize=12, color='red')
plt.legend()
plt.savefig(outfilename)
