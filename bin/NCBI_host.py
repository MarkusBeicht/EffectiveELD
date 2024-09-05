#!/usr/bin/env python3

import csv, os, sys

infilename=sys.argv[1]
outfilename=sys.argv[2]


def extract_subtype_data(full_string, data_string, seperator, field):
	if seperator+field+seperator in full_string:
		field_substrings=full_string.split(seperator)
		index=field_substrings.index(field)
		data_substrings=data_string.split(seperator)
		data_field=data_substrings[index]
	else:
		data_field="NA"
	return(data_field)


def get_hostid_from_hostname(hostname):
	if hostname =="NA":
		hostid=0
	else:
		command = "esearch -db taxonomy -query '" + hostname + "' | efetch -format docsum -stop 1 | xtract -pattern DocumentSummary -element TaxId"
		hostid = os.popen(command).read()[:-1]
		if hostid == "":
			hostid=0
		searched_hostnames[hostname] = hostid
	return(hostid)




#reads host-data input file
data = []
with open(infilename) as infile:
	for line in infile:
		data.append(line.split("\t"))




#loops through input data, calls functions and writes output to file
searched_hostnames = {}
with open(outfilename, "w") as outfile:
	for line in data:
		accession=line[0]
		name=line[1]+line[2].rstrip()
		NCBITaxonomyID=int(line[3])
		environment="NA"
		hostname=extract_subtype_data(line[4], line[5], "|", "host")
		hostid=0
		
		#checks if hostname has already been searched against the NCBI taxonomy database 
		if hostname in searched_hostnames:
			hostid=searched_hostnames[hostname]

		#searches new hostnames against the NCBI taxonomy database by calling get_hostid_from_hostname()
		else:
			hostid=get_hostid_from_hostname(hostname)

		outfile.write("%s\t%s\t%i\t%s\t%s\t%s\n" % (accession, name, NCBITaxonomyID, environment, hostname, hostid))
		sys.stderr.flush()
sys.stderr.flush()