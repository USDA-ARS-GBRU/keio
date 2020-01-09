# Welcome to keio tutorial !!!
These tutorials walk you through the process of analyzing illumina read for developing keio collection type results in our model organism. 


# Background

The Keio collection in  Escherichia coli K-12 represents  a collection of single-gene deleted mutants. This collection was created by manually one by one by replacing predicted ORF with a kanamycin cassette to inactivate chromosomal genes. Then primers were designed to create in-frame deletions upon excision of the resistance cassette. Of 4288 genes targeted, 3985 mutants were obtained and major of these mutants represents the mutation of non-essential genes. Summary table of mutant from the main Keio collection paper(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1681482/pdf/msb4100050.pdf).

![alt text](https://github.com/ravinpoudel/keio/blob/master/KEIO_mutant_summary.png)


Primarily,the KEIO collection provides a new molecular tool/resource to understand the functional and physiological aspects of gene at the system levels. 
Although, creating such molecular collection / tools takes lot of resorce and daunting. Thus, here we explore and RB-TnSeq (Randomly Barcode Transposons) method to crate a single gene mutant type collection in Phaeobacter_inhibens_DSM_17395. The molecular construct of Randomly Barcode Transposon is similar as follow:

<img src="https://github.com/ravinpoudel/keio/blob/master/RbTransposon.png" align="center" height="550" width="350"/>
 

Then, the Randomly Barcode Transposon are randomly inserted into bacterial genome to create mutated clone. Each clone therotically should represent a single gene mutation. 
 

<img src="https://github.com/ravinpoudel/keio/blob/master/RB_Clone.png" align="center" height="550" width="350" />


# What we need to do in our project?
Our Illumina reads represent the clone library generated similarly as above describe methods. Now, we need to map the location of random barcode sequence (about 20 base pair) in lengths. Following diagram repsent the construct for each reads:

![alt text](https://github.com/ravinpoudel/keio/blob/master/keio.png)








