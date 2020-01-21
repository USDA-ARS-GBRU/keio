# Welcome to Keio tutorial !!!
This tutorial walks you through the process of analyzing illumina read for developing Keio collection type results in our model organism: Phaeobacter_inhibens_DSM_17395. 


# Background

# Background
The Keio collection in Escherichia coli K-12 represents a collection of single-gene deleted mutants. This collection was created manually one by one by replacing the predicted ORF with a kanamycin cassette to inactivate chromosomal genes. Then primers were designed to create in-frame deletions upon excision of the resistance cassette. Of 4288 genes targeted, 3985 mutants were obtained. Majority of these mutants represents mutation of non-essential genes. Summary table of mutant from the Keio collection.

(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1681482/pdf/msb4100050.pdf).

![alt text](https://github.com/ravinpoudel/keio/blob/master/KEIO_mutant_summary.png)


Primarily, the Keio collection provides a molecular tool/resource to understand the functional and physiological aspects of gene at the system levels. Although, creating such molecular collection / tools takes lot of resources and time. Thus, here we explored RB-TnSeq (Randomly Barcode Transposons) method to create a single gene mutant type collection in Phaeobacter_inhibens_DSM_17395. The molecular construct of Randomly Barcode Transposon is similar as follow:

<img src="https://github.com/ravinpoudel/keio/blob/master/RbTransposon.png" align="center" height="550" width="350"/>
 

Once the Randomly Barcode Transposons are created, the constructs are randomly inserted into bacterial genome to create mutated clone. Each clone theoretically should represent a single gene mutation. 
 

<img src="https://github.com/ravinpoudel/keio/blob/master/RB_Clone.png" align="center" height="550" width="350" />

More information about the RB-TnSeq and the methods is described by Wetmore et. al.(2015) (https://mbio.asm.org/content/6/3/e00306-15).


# What we need to do in our project?
Our Illumina reads represent the clone library as generated with the above described methods. Now, we need to map the location of random barcode sequence (about 20 base pair in lengths). Following diagram shows the construct for each read:

![alt text](https://github.com/ravinpoudel/keio/blob/master/keio.png)


# Flowchart- Framework of analyses
Script involved in each step is described in the following flowchart, and also available in `main.py` script file.  

![alt text](https://github.com/ravinpoudel/keio/blob/master/Keio_Flowchart.png)



# Final output file

`df_mapping_with_geneinfo.csv` has the final output. Following table shows first 11 rows of the output table. 

|    | cutseq               | ref_barcode          | distance | cutseq_size | Plate_Number | Clone_Number | scaffold  | strand | pos     | reads | gene_info-[protein_name,locus_tag start, stop, protein_product, rb20_position]                                                                                               |
| -- | -------------------- | -------------------- | -------- | ----------- | ------------ | ------------ | --------- | ------ | ------- | ----- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 0  | CATGAGGCAGCGTGGATAGT | CATGAGGCAGCGTGGATAGT | 0        | 20          | 1            | Clone-7      | PGA1_c    | +      | 627691  | 18    |                                                                                                                                                                              |
| 1  | GGTCGACGTAGAGGGCGGGG | GGTACGCGTAGAGGGCGGTG | 3        | 20          | 1            | Clone-11     | PGA1_c    | -      | 3215303 | 36    | [('WP_014881201.1', 'PGA1_RS15290', 3214836, 3215348, 'hypothetical protein', 3215303)]                                                                                      |
| 2  | TGTGCCAAATGCGGGGCACC | TGTGCCAAATGCGGGGCACC | 0        | 20          | 1            | Clone-13     | PGA1_c    | +      | 486076  | 55    |                                                                                                                                                                              |
| 3  | GCGACAGGCAGCACTAGGTT | GCGACAGGCAGCACTAGGTT | 0        | 20          | 1            | Clone-15     | PGA1_c    | +      | 2307030 | 35    |                                                                                                                                                                              |
| 4  | TTGTTCTGAGTCTAAGCAGT | TTTGTTCTGGTTCAAGCAGT | 4        | 20          | 1            | Clone-25     | PGA1_c    | -      | 747956  | 38    | [('WP_014879339.1', 'PGA1_RS03600', 747173, 748378, 'Gfo/Idh/MocA family oxidoreductase', 747956)]                                                                           |
| 5  | ATGATGTTAGGGATCATAGA | ATGATGTTAGGGATCATAGA | 0        | 20          | 1            | Clone-27     | PGA1_262p | +      | 215693  | 78    | [('WP_014878956.1', 'PGA1_RS01155', 215183, 216454, 'MFS transporter permease', 215693), ('WP_014881826.1', 'PGA1_RS19085', 215123, 215875, 'hypothetical protein', 215693)] |
| 6  | AAGGATGTAAGAGCTATCGA | AAGGATGTAAGAGCTATCGA | 0        | 20          | 1            | Clone-29     | PGA1_c    | -      | 2726842 | 125   | [('WP_014875518.1', 'PGA1_RS13000', 2726604, 2726861, 'hypothetical protein', 2726842)]                                                                                      |
| 7  | GCCTTACAACCGGCTGTGCT | GCCTTACAACCGGCTGTGCT | 0        | 20          | 1            | Clone-30     | PGA1_c    | +      | 706052  | 95    | [('WP_014879308.1', 'PGA1_RS03405', 705621, 706442, 'hypothetical protein', 706052)]                                                                                         |
| 8  | AGGGTGCTTTGTGCCGGGGG | AGGGTGCTTTGTGCCGGGGG | 0        | 20          | 1            | Clone-31     | PGA1_c    | +      | 1410032 | 202   | [('WP_014879849.1', 'PGA1_RS06795', 1409400, 1410080, 'hypothetical protein', 1410032)]                                                                                      |
| 9  | CCGGTCGCGGCCCGCGGGGA | CCGGTCGCGGCCCGCGGGGA | 0        | 20          | 1            | Clone-32     | PGA1_c    | -      | 2298573 | 156   | [('WP_014880515.1', 'PGA1_RS10980', 2297921, 2298847, 'alpha/beta hydrolase', 2298573)]                                                                                      |
| 10 | CTTGTTGCGATGGGTGGAGG | GCATGGTGGATGGGTGGAGG | 4        | 20          | 1            | Clone-39     | PGA1_c    | +      | 3215295 | 15    | [('WP_014881201.1', 'PGA1_RS15290', 3214836, 3215348, 'hypothetical protein', 3215295)]                                                                                      |


### Legend for the output

| Column Header                                                               | Information                                                                    |
| --------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| cutseq                                                                      | random barcode from this analyses                                              |
| ref_barcode                                                                 | random barcode present in "Phaeo_ML1.loconf.pool.txt"                          |
| distance                                                                    | leven distance; based on KNN (NMSLIB)                                          |
| cutseq_size                                                                 | size of random barcode                                                         |
| Plate_Number                                                                | Plate from where the clone was selected                                        |
| Clone_Number                                                                | Specific clone number in a provided Plate_Number                               |
| scaffold                                                                    | information from "Phaeo_ML1.loconf.pool.txt"                                   |
| strand                                                                      | + / -                                                                     |
| pos                                                                         | start position of gene                                                         |
| reads                                                                       | Total number of reads with the given cutseq observed in the current experiment |
| gene_info-[protein_name;locus_tag;start;stop;protein_product;rb20_position] | additional information                                                         |
