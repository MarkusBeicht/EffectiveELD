#!/usr/bin/env python3

import sys, os.path, numpy
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, confusion_matrix


ELDfilename=sys.argv[1]
predictionfilename=sys.argv[2]
outfilename=sys.argv[3]
mean_filter=float(sys.argv[4])
min_effectsize=float(sys.argv[5])


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


#loading Wilcoxon signed-rank test p-values from euk_score_ELD file
with open(predictionfilename) as infile:
  for line in infile:
    (PfamID, refseq_accession, pvalue, r, nonsym_mean) = line.strip().split("\t")
    p=float(pvalue)
    nonsym_mean = float(nonsym_mean)
    r = float(r)
    if mean_filter != 0 and nonsym_mean > mean_filter:
      p = 1
    if min_effectsize != 0 and r < min_effectsize:
      p = 1
    InterproID = next((InterproID for InterproID, PfamIDs in InterproPfamMapping.items() if PfamID in PfamIDs))
    if p < min_pvalue[(InterproID, refseq_accession)]:
      min_pvalue[(InterproID, refseq_accession)]=p



val = numpy.array(list(ELDcounts.values()))
prob = 1 - numpy.array(list(min_pvalue.values()))
fpr, tpr, thresholds = roc_curve(val, prob, pos_label=1)
roc_auc = auc(fpr, tpr)
print('AUC: ', roc_auc)

plt.plot(fpr, tpr, label=f'ROC (AUC = {roc_auc:.4f})')
plt.plot([0, 1], [0, 1], 'r--', label='Random Guess')

# Set labels and title
if mean_filter == 0 and min_effectsize == 0:
  plt.title('ROC Curve of the EffectiveELD scoring method')
elif mean_filter != 0 and min_effectsize == 0:
  plt.title('ROC Curve of the EffectiveELD scoring method\n(maximal mean domain count in non-symbionts of %s)' % mean_filter)
elif mean_filter == 0 and min_effectsize != 0:
  plt.title('ROC Curve of the EffectiveELD scoring method\n(minimal effect size of r = %s)' % min_effectsize)
elif mean_filter != 0 and min_effectsize != 0:
  plt.title('ROC Curve of the EffectiveELD scoring method\n(maximal mean domain count in non-symbionts of %s)\n(minimal effect size of r = %s)' % mean_filter, min_effectsize)


plt.xlabel('False Positive Rate (1 - Specificity)')
plt.ylabel('True Positive Rate (Sensitivity)')
plt.text(2, 7, f"AUC: {roc_auc:.2f}", fontsize=12, color='red')
plt.legend()
plt.savefig(outfilename)


#calculate statistical key figures for threshold (p = 0.05)
def confusionmatrix_values(val, prob, p):
  threshold = 1 - p
  prob_p = numpy.where(prob >= threshold, 1, 0)
  tn, fp, fn, tp = confusion_matrix(val, prob_p).ravel()
  sensitivity = tp / (tp + fn)
  specificity = tn / (tn + fp)
  PPV = tp / (tp + fp)
  NPV = tn / (tn + fn)
  accuracy = (tp + tn) / (tp + tn + fp + fn)

  print(['P-value', 'Sensitivity', 'Specificity', 'PPV', 'NPV', 'Accuracy', 'tp', 'tn', 'fp', 'fn'])
  print([p, float(sensitivity), float(specificity), float(PPV), float(NPV), float(accuracy), int(tp), int(tn), int(fp), int(fn)])

confusionmatrix_values(val, prob, 0.05)