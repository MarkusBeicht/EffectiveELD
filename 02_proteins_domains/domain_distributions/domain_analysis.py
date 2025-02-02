#!/usr/bin/env python3

import sys, os.path, numpy, scipy.stats, math, csv, random
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import statsmodels.api as sm

domaincountfilename=sys.argv[1]
outfilename=sys.argv[2]
outdirname=sys.argv[3]

#loading nonsymbiotic domaincounts from file
domaincounts={}
with open(domaincountfilename) as infile:
  for line in infile:
    (refseq_accession, domainid, count) = line.strip().split("\t")
    domaincount=int(count)
    if domaincounts.get(domainid) == None:
      domaincounts[domainid]=[domaincount]
    else:
      domaincounts[domainid].append(domaincount)
sys.stderr.write('Got the counts of %i eukaryotic domains for %i genomes.\n' % (len(domaincounts.keys()), len(next(iter(domaincounts.values())))))



#calculating statistical values
def statistical_key_figures(domainid, counts):
  count = len(counts)
  non_zero_counts = len(counts) - counts.count(0)
  max_count = max(counts)
  min_count = min(counts)
  mean = numpy.mean(counts)
  stddev = numpy.std(counts)
  kurtosis = scipy.stats.kurtosis(counts)
  skew = scipy.stats.skew(counts)

  if mean == 0:
  #print(domainid, count, 0, 0, 0.0, 0.0, 0.0, 0.0, "nan", "nan", "nan", 0.0001, "nan", 0.0001, "nan", 0.0001, "nan", 0.0001, "nan", 0.0001, 1.0)
    outfile.write("%s\t%i\t%i\t%i\t%i\t%f\t%f\t%s\t%s\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t%f\n" % (domainid, 0, 0, 0, 0, 0.0, 0.0, "nan", "nan", "nan", 0.0001, "nan", 0.0001, "nan", 0.0001, "nan", 0.0001, "nan", 0.0001, 1.0))
    sys.stdout.flush()
  else:
    #ks_test for normality 
    rng = numpy.random.default_rng()
    norm = scipy.stats.goodness_of_fit(scipy.stats.norm, counts, statistic='ks', random_state=rng)	

    #Log transformation and ks-test
    counts_log = [math.log10(x+1) for x in counts]
    norm_log = scipy.stats.goodness_of_fit(scipy.stats.norm, counts_log, statistic='ks', random_state=rng)	
	
    #Sqrt transformation and ks-test
    counts_sqrt = [math.sqrt(x) for x in counts]
    norm_sqrt = scipy.stats.goodness_of_fit(scipy.stats.norm, counts_sqrt, statistic='ks', random_state=rng)	

    #Reciprocal transformation and ks-test
    counts_cbrt = [math.cbrt(x) for x in counts]
    norm_cbrt = scipy.stats.goodness_of_fit(scipy.stats.norm, counts_cbrt, statistic='ks', random_state=rng)	

    #Reciprocal transformation and ks-test
    counts_yeojohnson, lmbda = scipy.stats.yeojohnson(counts)
    norm_yeojohnson = scipy.stats.goodness_of_fit(scipy.stats.norm, counts_yeojohnson, statistic='ks', random_state=rng)	


    #writing row with statistical values to result file
    #print(domainid, non_zero_counts, min_count, max_count, mean, stddev, median, iqr, kurtosis, skew, norm.statistic, norm.pvalue, norm_log.statistic, norm_log.pvalue, norm_sqrt.statistic, norm_sqrt.pvalue, norm_cbrt.statistic, norm_cbrt.pvalue, norm_yeojohnson.statistic, norm_yeojohnson.pvalue, lmbda)
    outfile.write("%s\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (domainid, count, non_zero_counts, min_count, max_count, mean, stddev, kurtosis, skew, norm.statistic, norm.pvalue, norm_log.statistic, norm_log.pvalue, norm_sqrt.statistic, norm_sqrt.pvalue, norm_cbrt.statistic, norm_cbrt.pvalue, norm_yeojohnson.statistic, norm_yeojohnson.pvalue, lmbda))
    sys.stdout.flush()

sys.stderr.write("\n") 




#visualize domaincount distributions in plots
def domaincount_distribution_plot(domain):
  data = domaincounts[domain]
  mean = numpy.mean(data)
  stddev = numpy.std(data)


  #defining plot
  fig = plt.figure(figsize=(30,20), constrained_layout=True)
  spec = fig.add_gridspec(ncols=2, nrows=1)
  fig.suptitle("%s" % (domain),fontsize=48)


  #domaincount distribuition subplot
  bin_number = len(set(data))
  bins= [x -0.5 for x in range(bin_number+1)]
  ticks = numpy.arange(bin_number)
  x = numpy.linspace(mean - 3 * stddev, mean + 3 * stddev, 100)
  skew = round(scipy.stats.skew(data), 2)
  kurtosis = round(scipy.stats.kurtosis(data), 2)
  textstr = '\n'.join(('skew = ' + str(skew), 'kurtosis = ' + str(kurtosis)))

  ax1 = fig.add_subplot(spec[0, 0])
  ax1.hist(data, bins=bins, density=True)
  ax1.set_xlabel("Domaincount", fontsize = 48)
  ax1.set_ylabel("Probability Density", fontsize = 48)
  ax1.set_title("Domaincount Distribution", fontsize = 48)
  ax1.tick_params(axis='both', which='major', labelsize=32)
  ax1.plot(x, scipy.stats.norm.pdf(x, mean, stddev), label="normal", color="red", linewidth=10)
  ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax1.legend(loc='upper right', shadow=True, fontsize = 48)
  ax1.text(0.95, 0.75, textstr, transform=ax1.transAxes, fontsize=48,va='center',ha='right')


  #qq-subplot
  data_ms = [(x-mean)/stddev for x in data]
  data_np = numpy.asarray(data_ms)
  ax2 = fig.add_subplot(spec[0, 1])
  sm.qqplot(data_np, line ='45', ax=ax2, markersize=10)

  ax2.set_title("Q-Q-Plot", fontsize = 48)
  ax2.set_xlabel("Theoretical Quantiles", fontsize = 48)
  ax2.set_ylabel("Sample Quantiles", fontsize = 48)
  ax2.get_lines()[1].set_linewidth("8")
  ax2.tick_params(axis='both', which='major', labelsize=32)
	
  plt.savefig('%s%s.png' % (outdirname, domain))
  plt.close()


#loop through the domaincounts of every domain and write to csv
outfile=open(outfilename, 'w')
outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("domainid", "n", "non_zero_counts", "min_count", "max_count", "mean", "stddev", "kurtosis", "skew", "ks_test_statistic", "ks_test_pvalue", "norm_log_statistic", "norm_log_pvalue", "norm_sqrt_statistic", "norm_sqrt_pvalue", "norm_cbrt_statistic", "norm_cbrt_pvalue", "norm_yeojohnson_statistic", "norm_yeojohnson_pvalue", "yeojohnson_lmbda"))

#for domainid, counts in domaincounts.items():
#  statistical_key_figures(domainid, counts)


#loops through all domains for plots
domains = list(domaincounts.keys())
for domain in domains:
  domaincount_distribution_plot(domain)

