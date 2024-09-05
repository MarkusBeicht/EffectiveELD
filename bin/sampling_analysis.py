#!/usr/bin/env python3

import sys, gzip, os, pandas
from statsmodels.stats.contingency_tables import mcnemar 
import matplotlib.pyplot as plt


sampling_infile_name=sys.argv[1]
sampling_outfile_name=sys.argv[2]

class MatchAll:
    def __eq__(self, other):
        return True


#compares data from 1 column with the study_classifications and writes them in new column in the Dataframe
def new_results1(df, column, name, s_vals, n_vals):
  results=[]
  for index in df.index:
    #symbiotic accoring to literature
    if df.loc[index, "study_classification"] == "s":
      if df.loc[index, column] in s_vals:
        results.append("TP")
      elif df.loc[index, column] in n_vals:
        results.append("FN")
      else:
        results.append("A")

    #non-symbiotic accoring to literature    
    elif df.loc[index, "study_classification"] == "n":
      if df.loc[index, column] in s_vals:
        results.append("FP")
      elif df.loc[index, column] in n_vals:
        results.append("TN")
      else:
        results.append("A")

  df[name] = results
  return df



#compares data from 3 columns, if the values is provided in at least 1 of them, with the study_classifications and writes them in new column in the Dataframe
def new_results3(df, column1, column2, column3, name, s_vals, n_vals):
  results=[]
  for index in df.index:
    
    #symbiotic accoring to literature
    if df.loc[index, "study_classification"] == "s":
      if df.loc[index, column1] in s_vals or df.loc[index, column2] in s_vals or df.loc[index, column3] in s_vals:
        results.append("TP")
      elif df.loc[index, column1] in n_vals or df.loc[index, column2] in n_vals or df.loc[index, column3] not in n_vals:
        results.append("FN")
      else:
        results.append("A")

    #non-symbiotic accoring to literature    
    elif df.loc[index, "study_classification"] == "n":
      if df.loc[index, column1] in s_vals or df.loc[index, column2] in s_vals or df.loc[index, column3] in s_vals:
        results.append("FP")
      elif df.loc[index, column1] in n_vals or df.loc[index, column2] in n_vals or df.loc[index, column3] not in n_vals:
        results.append("TN")
      else:
        results.append("A")

  df[name] = results
  return df



#Adds the results of a PhenDB prediction to the database results, compares them with the study_classifications and writes them in new column in the Dataframe
def new_results_db1(df, columndb, column1, name, db_s_vals, db_n_vals, phen_s_vals):
  results=[]
  for index in df.index:
    
    #symbiotic accoring to literature
    if df.loc[index, "study_classification"] == "s":
      if df.loc[index, columndb] in db_s_vals or df.loc[index, column1] in phen_s_vals:
        results.append("TP")
      elif df.loc[index, columndb] in db_n_vals:
        results.append("FN")
      else:
        results.append("A")

    #non-symbiotic accoring to literature    
    elif df.loc[index, "study_classification"] == "n":
      if df.loc[index, columndb] in db_s_vals or df.loc[index, column1] in phen_s_vals:
        results.append("FP")
      elif df.loc[index, columndb] in db_n_vals:
        results.append("TN")
      else:
        results.append("A")

  df[name] = results
  return df



#Adds the results of 3 PhenDB predictions (at least 1 secretion system) to the database results, compares them with the study_classifications and writes them in new column in the Dataframe
def new_results_db3(df, columndb, column1, column2, column3, name, db_s_vals, db_n_vals, phen_s_vals):
  results=[]
  for index in df.index:
    
    #symbiotic accoring to literature
    if df.loc[index, "study_classification"] == "s":
      if df.loc[index, columndb] in db_s_vals or df.loc[index, column1] in phen_s_vals or df.loc[index, column2] in phen_s_vals or df.loc[index, column3] in phen_s_vals:
        results.append("TP")
      elif df.loc[index, columndb] in db_n_vals:

        results.append("FN")
      else:
        results.append("A")

    #non-symbiotic accoring to literature    
    elif df.loc[index, "study_classification"] == "n":
      if df.loc[index, columndb] in db_s_vals or df.loc[index, column1] in phen_s_vals or df.loc[index, column2] in phen_s_vals or df.loc[index, column3] in phen_s_vals:
        results.append("FP")
      elif df.loc[index, columndb] in db_n_vals:
        results.append("TN")
      else:
        results.append("A")

  df[name] = results
  return df



#calculates the count, sensitivity, specificity, accuracy of the results_columns
def calc_values(df, column, database):
  filtered_df = df[df["database"] == database]
  TP = (filtered_df[column] == "TP").sum()
  TN = (filtered_df[column] == "TN").sum()
  FP = (filtered_df[column] == "FP").sum()
  FN = (filtered_df[column] == "FN").sum()
  A = (filtered_df[column] == "A").sum()

  count = TP+TN+FP+FN+A
  count_classified = TP+TN+FP+FN
  classified_percent = 100*count_classified/count
  sensitivity = TP / (TP + FN)
  specificity = TN / (TN + FP)
  accuracy = (TP + TN) / (TP + TN + FP + FN)
  if database == "<__main__.MatchAll object at 0x7fdd756e9160>":
    database = "all_databases"
  return column, database, count_classified.item(), round(classified_percent.item(),2), round(sensitivity.item(),2), round(specificity.item(),2), round(accuracy.item(),2)


#Performs exact mcnemar test on the results of two classifications to test whether the second one is significantly better. 
def mcnemar_func(df, before, after):
  true_true = 0
  true_false = 0
  false_true = 0
  false_false = 0
  
  for index in df.index:
    if df.loc[index, before] == "TP" or df.loc[index, before] == "TN":
      if df.loc[index, after] == "TP" or df.loc[index, after] == "TN":
        true_true += 1
      elif df.loc[index, after] == "FP" or df.loc[index, after] == "FN":
        true_false += 1
    elif df.loc[index, before] == "FP" or df.loc[index, before] == "FN":
      if df.loc[index, after] == "TP" or df.loc[index, after] == "TN":
        false_true += 1
      elif df.loc[index, after] == "FP" or df.loc[index, after] == "FN":
        false_false += 1
        
  #Checks if the results were improved (more false_to_true changes than true_to_false chnages in the classifications)
    data = [[true_true,   true_false], [false_true, false_false]]
    result = mcnemar(data, exact=True)
    #Adapts the p-value to fit the one-sided version of the test
    one_sided_p_value = result.pvalue.item() / 2
  else:
    one_sided_p_value = 1.00
  return after, round(one_sided_p_value,2), true_true, true_false, false_true, false_false




#reading GOLD dataset to extract Bacteria and their ecosystem
sampling_data=[]
with open(sampling_infile_name, encoding='iso-8859-1') as infile:
  next(infile)
  for line in infile:
    (database, name, environmental_metadata, hostID, metadata_classification, study_citation, study_classification, refseq_accession, symbiont, T3SS, T4SS, T6SS) = line.strip().split("\t")
    sampling_data.append([database, name, environmental_metadata, hostID, metadata_classification, study_citation, study_classification, refseq_accession, symbiont, T3SS, T4SS, T6SS])


header=["database", "name", "environmental_metadata", "hostID", "metadata_classification", "study_citation", "study_classification", "refseq_accession", "symbiont", "T3SS", "T4SS", "T6SS"]
df = pandas.DataFrame(sampling_data, columns=header)


df=new_results1(df, "metadata_classification", "database_result", ["s"], ["n"])
df=new_results1(df, "T3SS", "T3SS_result", ["+"], ["-"])
df=new_results1(df, "T4SS", "T4SS_result", ["+"], ["-"])
df=new_results1(df, "T6SS", "T6SS_result", ["+"], ["-"])
df=new_results1(df, "symbiont", "symbiont_result", ["+"], ["-"])
df=new_results3(df, "T3SS", "T4SS", "T6SS", "at_least_1_SS_result", ["+"], ["-"])
df=new_results_db1(df, "metadata_classification", "T3SS", "db_T3SS_result", ["s"], ["n"], ["+"])
df=new_results_db1(df, "metadata_classification", "T4SS", "db_T4SS_result", ["s"], ["n"], ["+"])
df=new_results_db1(df, "metadata_classification", "T6SS", "db_T6SS_result", ["s"], ["n"], ["+"])
df=new_results_db1(df, "metadata_classification", "symbiont", "db_symbiont_result", ["s"], ["n"], ["+"])
df=new_results_db3(df, "metadata_classification", "T3SS", "T4SS", "T6SS", "db_at_least_1_SS_result", ["s"], ["n"], ["+"])

df.to_csv(sampling_outfile_name, sep="\t") 


GOLD = (calc_values(df, "database_result", "GOLD"))
NCBI = (calc_values(df, "database_result", "NCBI"))
BacMap = (calc_values(df, "database_result", "BacMap"))
all_databases = (calc_values(df, "database_result", MatchAll()))

T3SS = (calc_values(df, "T3SS_result", MatchAll()))
T4SS = (calc_values(df, "T4SS_result", MatchAll()))
T6SS = (calc_values(df, "T6SS_result", MatchAll()))
SYMBIONT = (calc_values(df, "symbiont_result", MatchAll()))
one_SS = (calc_values(df, "at_least_1_SS_result", MatchAll()))

db_T3SS = (calc_values(df, "db_T3SS_result", MatchAll()))
db_T4SS = (calc_values(df, "db_T4SS_result", MatchAll()))
db_T6SS = (calc_values(df, "db_T6SS_result", MatchAll()))
db_SYMBIONT = (calc_values(df, "db_symbiont_result", MatchAll()))
db_one_SS = (calc_values(df, "db_at_least_1_SS_result", MatchAll()))


print(("result_type", "database", "classified_instances", "classified_%", "sensitivity", "specificity", "accuracy"))
print(GOLD)
print(NCBI)
print(BacMap)
print(all_databases)
print(T3SS)
print(T4SS)
print(T6SS)
print(SYMBIONT)
print(one_SS)
print(db_T3SS)
print(db_T4SS)
print(db_T4SS)
print(db_T6SS)
print(db_SYMBIONT)
print(db_one_SS)

print("\nmcnemar tests")
print(mcnemar_func(df, "database_result", "db_T3SS_result"))
print(mcnemar_func(df, "database_result", "db_T4SS_result"))
print(mcnemar_func(df, "database_result", "db_T6SS_result"))
print(mcnemar_func(df, "database_result", "db_symbiont_result"))
print(mcnemar_func(df, "database_result", "db_at_least_1_SS_result"))


results1 = [GOLD, NCBI, BacMap, all_databases]
results2 = [T3SS, T4SS ,T6SS, SYMBIONT, one_SS]
results3 = [db_T3SS, db_T4SS, db_T6SS, db_SYMBIONT, db_one_SS]
result_names1 = ['GOLD', 'NCBI', 'BacMap', 'all 3 databases']
result_names2 = ['T3SS','T4SS' ,'T6SS', 'SYMBIONT', 'one SS']
result_names3 = ['EM + T3SS', 'EM + T4SS', 'EM + T6SS', 'EM + SYMBIONT', 'EM + one SS']
colors1 = ['cyan', 'olive', 'brown', 'pink']
colors2 = ['red', 'blue', 'green', 'orange', 'purple']



for i, result in enumerate(results1):
  plt.scatter(result[5], result[4], label=result_names1[i], marker='.', s=200, color=colors1[i])
for i, result in enumerate(results2):
  plt.scatter(result[5], result[4], label=result_names2[i], marker='^', s=100, color=colors2[i])
for i, result in enumerate(results3):
  plt.scatter(result[5], result[4], label=result_names3[i], marker='s', s=100, color=colors2[i])


plt.xlabel('Specificity')
plt.ylabel('Sensitivity')
#plt.title('Comparison of the different methods for the classification of bacteria')
plt.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=1)
plt.tight_layout()
plt.savefig('sampling_results.png')
plt.show()
