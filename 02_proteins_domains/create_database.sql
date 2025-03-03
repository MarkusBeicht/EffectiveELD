#Creating tables
#CREATE table genome (refseq_accession varchar(255) NOT NULL, taxonomyid int(20) NOT NULL, description varchar(255), class varchar(1), PRIMARY KEY(refseq_accession)); 
#CREATE table euk_domain (domainid VARCHAR(7) NOT NULL, description varchar(255), median INT(5) NOT NULL, iqr INT(5) NOT NULL, threshold_0 INT(5) NOT NULL, threshold_1 INT(5) NOT NULL, threshold_3 INT(5) NOT NULL, threshold_5 INT(5) NOT NULL, threshold_7 INT(5) NOT NULL, PRIMARY KEY(domainid));
#CREATE table protein (refseq_accession varchar(255) NOT NULL, protein_accession VARCHAR(50) NOT NULL, description varchar(255), PRIMARY KEY(refseq_accession, protein_accession)); 
#CREATE TABLE euk_score (domainid VARCHAR(7) NOT NULL, refseq_accession varchar(255) NOT NULL, max_threshold DECIMAL(2,1) NOT NULL, PRIMARY KEY(refseq_accession, domainid));
#CREATE TABLE euk_prediction (refseq_accession VARCHAR(255) NOT NULL, protein_accession VARCHAR(50) NOT NULL, domainid VARCHAR(7) NOT NULL, PRIMARY KEY(refseq_accession, protein_accession, domainid)); 

#CREATE TABLE genome_NF3 (refseq_accession varchar(255) NOT NULL, taxonomyid int(20) NOT NULL, description varchar(255) DEFAULT NULL, PRIMARY KEY (refseq_accession), KEY taxonomyid (taxonomyid)); 
#CREATE TABLE classification_NF3 (taxonomyid int(20) NOT NULL, class varchar(1), PRIMARY KEY(taxonomyid), KEY(class)); 


#Loading data
#LOAD DATA local INFILE 'C:/Users/marku/Desktop/genome' INTO TABLE genome;
#LOAD DATA local INFILE 'C:/Users/marku/Desktop/euk_domain' INTO TABLE euk_domain;
#LOAD DATA local INFILE 'C:/Users/marku/Desktop/protein' INTO TABLE protein;
#LOAD DATA local INFILE 'C:/Users/marku/Desktop/euk_score' INTO TABLE euk_score;
#LOAD DATA local INFILE 'C:/Users/marku/Desktop/euk_prediction' INTO TABLE euk_prediction;

#Creating indices (in addition to primary-key indices)
#CREATE INDEX taxonomyid ON genome(taxonomyid);
#CREATE INDEX class ON genome(class);
#CREATE INDEX description ON genome(description);

#CREATE INDEX description ON euk_domain(description);
#CREATE INDEX threshold_0 ON euk_domain(threshold_0);
#CREATE INDEX threshold_1 ON euk_domain(threshold_1);
#CREATE INDEX threshold_3 ON euk_domain(threshold_3);
#CREATE INDEX threshold_5 ON euk_domain(threshold_5);
#CREATE INDEX threshold_7 ON euk_domain(threshold_7);

#CREATE INDEX domainid ON euk_score(domainid);
#CREATE INDEX refseq_accession ON euk_score(refseq_accession);
#CREATE INDEX max_threshold ON euk_score(max_threshold);

#CREATE INDEX description ON protein(DESCRIPTION);
#CREATE INDEX refseq_accession ON protein(refseq_accession);
#CREATE INDEX protein_accession ON protein(protein_accession);

#CREATE INDEX refseq_accession_domainid ON euk_prediction(refseq_accession, domainid);
#CREATE INDEX refseq_accession_protein_accession ON euk_prediction(refseq_accession, protein_accession);


