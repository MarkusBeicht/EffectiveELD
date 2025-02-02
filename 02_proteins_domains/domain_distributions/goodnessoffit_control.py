#!/usr/bin/env python3

import numpy, scipy.stats, sys
outfilename=sys.argv[1]

outfile=open(outfilename, 'w')
outfile.write("%s\t%s\t%s\t%s\n" % ("mean", "stddev", "ks_statistic", "ks_pvalue"))


#tests goodness of fit test x times with normally distributed data with the sample size n
x = 100
n = 7059

for i in range(x):
	#create random normal distribution with varying means and stddev
	mu = numpy.random.randint(-100, 100)
	sigma = numpy.random.randint(0, 100)
	data = numpy.random.normal(mu, sigma, n)

	mean = numpy.mean(data)
	stddev = numpy.std(data)
	
	#perform goodness of fit test (Kolmogornov Smirnov)
	rng = numpy.random.default_rng()
	norm = scipy.stats.goodness_of_fit(scipy.stats.norm, data, statistic='ks', random_state=rng)			
	outfile.write("%f\t%f\t%f\t%f\n" % (mean, stddev, norm.statistic, norm.pvalue))
	
