
#ELDs and effectors in Legionella pneumophila Paris
WITH temp AS (
SELECT genome.refseq_accession, genome.`description` AS genome_description, euk_score.domainid, euk_score.max_threshold, euk_domain.`description` AS domain_description, euk_domain.`median`, euk_domain.iqr, protein.protein_accession, protein.`description` AS protein_description

FROM genome 

LEFT JOIN euk_score
	ON euk_score.refseq_accession = genome.refseq_accession

LEFT JOIN euk_domain
	ON euk_score.domainid = euk_domain.domainid

LEFT JOIN euk_prediction
	ON euk_score.refseq_accession = euk_prediction.refseq_accession
	AND euk_score.domainid = euk_prediction.domainid

LEFT JOIN protein
	ON euk_prediction.refseq_accession = protein.refseq_accession
	AND euk_prediction.protein_accession = protein.protein_accession

WHERE genome.refseq_accession = 'GCF_000048645.1'
#AND euk_score.max_threshold IN (0.7)
#AND euk_domain.`median` IN (0)
)

SELECT distinct temp.refseq_accession, temp.domainid, temp.domain_description, temp.median, temp.iqr, temp.max_threshold, temp.protein_accession, temp.protein_description
FROM temp;



/*
#number of ELDs predited with varying effect sizes
WITH temp AS (
SELECT euk_score.domainid, euk_score.max_threshold
FROM euk_score
#WHERE euk_score.refseq_accession = 'GCF_000048645.1'
)

SELECT temp.max_threshold, COUNT(*) ANZ
FROM temp
GROUP BY temp.max_threshold;
*/


/*
#number of effectors predited with varying effect sizes
WITH temp AS (
#SELECT distinct euk_score.max_threshold, protein.protein_accession
SELECT distinct protein.protein_accession

FROM euk_score
LEFT JOIN euk_domain
	ON euk_score.domainid = euk_domain.domainid
LEFT JOIN euk_prediction
	ON euk_score.refseq_accession = euk_prediction.refseq_accession
	AND euk_score.domainid = euk_prediction.domainid
LEFT JOIN protein
	ON euk_prediction.refseq_accession = protein.refseq_accession
	AND euk_prediction.protein_accession = protein.protein_accession

WHERE euk_score.refseq_accession = 'GCF_000048645.1'
)

#SELECT temp.max_threshold, COUNT(*) ANZ
SELECT *
FROM temp
#GROUP BY temp.max_threshold;
*/



/*
#Medians of known ELDs
SELECT * 
FROM euk_domain
WHERE euk_domain.domainid IN (
'PF00856' #SET
,'PF08123' #DOT1
,'PF02201' #SWIB/MDM2
,'PF04564' #U-box
,'PF00646', 'PF12937', 'PF13013', 'PF15966' #F-box
,'PF02338' #OTU deubiquitinase 
,'PF02902' #ULP1/Peptidase_C48_C
,'PF00651' #BTB/POZ domain TYPE
,'PF02985' #HEAT repeat
,'PF00415' #RCC1
,'PF00560', 'PF13516', 'PF13855' #LRR
,'PF00400' #WD40 repeat
,'PF01535', 'PF12854', 'PF13041', 'PF13812' #PPR
,'PF02493' #MORN
,'PF08238' #Sel1-like repeat
,'PF13174', 'PF13176', 'PF13181' #TPR
,'PF00023', 'PF12796', 'PF13606', 'PF13637', 'PF13857' #ANK
)
ORDER BY euk_domain.domainid;
*/