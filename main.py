#!/usr/bin/env python
"""
KEIO: A python module to process illumina reads for keio-collection type project.

"""
# Author: Ravin Poudel


#!/usr/bin/env python
# Author: Ravin Poudel

import os
import sys
import textdistance
import hashlib
import subprocess
import Bio
import pickle
import nmslib
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.SeqRecord import SeqRecord


def run_vsearch(maping_fasta, reads_fasta, cluster_id=0.75):
    """ Returns mapping information
    Args:
        input1(str): barcodefile: fasta
            input2(str): reads: fastafile
    input3 (int): cluster_id
    Returns:
        output: file: vsearch output containing alignment positon and quality
    """
    try:
        out_info = 'query+target+ql+tl+id+tcov+qcov+ids+gaps+qrow+trow+id4+qilo+qihi+qstrand+tstrand'
        outputfile = maping_fasta.split(".")[0] + "__output.txt"
        parameters = ["vsearch", "--usearch_global", reads_fasta, "--db", maping_fasta, "--id", str(cluster_id), "--userfield", out_info, "--strand", "plus", "--userout", outputfile]
        p0 = subprocess.run(parameters, stderr=subprocess.PIPE)
        print(p0.stderr.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(str(e))


def parse_vsearch(file, strand="+"):
    """Parse vsearch file returning a dictionary of top hits for each primer and seq.

    Args:
        input (str): file: the input file to be parsed

    Returns:
        (dict): A dictionary, e.g. {seq1: {primer1: {spos,epos,pmatch,tcov, gaps},...},...}
    """
    sdict = {}
    with open(file, 'r') as f:
        for line in f:
            ll = line.strip().split()
            qname = ll[0]
            tname = ll[1]
            pmatch = float(ll[4])
            tcov = float(ll[5])
            gaps = int(ll[8])
            spos = int(ll[12])
            epos = int(ll[13])
            qstrand = ll[14]
            tstrand = ll[15]
            if qstrand == strand and tcov >= 50:
                if gaps < 5:
                    td = {'spos': spos, 'epos': epos, 'pmatch': pmatch, 'tcov': tcov,
                          'gaps': gaps, 'qstrand': qstrand, 'tstrand': tstrand}
                    if qname in sdict:
                        if tname in sdict[qname]:
                            if sdict[qname][tname]['pmatch'] < pmatch:
                                sdict[qname][tname] = td
                        else:
                            sdict[qname][tname] = td
                    else:
                        sdict[qname] = {}
                        sdict[qname][tname] = td
        return sdict

# filter the parsed vserach file based on the number of matching barcode type
def filter_vsearch(sdict, nhits):
    outdict = {}
    for items in sdict.items():
        if len(items[1]) >= nhits:
            outdict[items[0]] = items[1]
    return outdict

# Now create a dictionary with start and end position for each barcode type
def get_randombarcode(keio_fasta, filter_vsearch_dict):
    out_dict={}
    for seq_record in SeqIO.parse(keio_fasta, "fasta"):
        if seq_record.id in filter_vsearch_dict.keys():
            sequence = str(seq_record.seq)
            listkey = list(filter_vsearch_dict[seq_record.id].keys())
            a = filter_vsearch_dict[seq_record.id][listkey[0]]['spos'] # start position for Ffasta
            b = filter_vsearch_dict[seq_record.id][listkey[0]]['epos'] # start position for Ffasta
            c = filter_vsearch_dict[seq_record.id][listkey[1]]['spos']
            d = filter_vsearch_dict[seq_record.id][listkey[1]]['epos']
            e = filter_vsearch_dict[seq_record.id][listkey[2]]['spos']
            f = filter_vsearch_dict[seq_record.id][listkey[2]]['epos']
            F_Index = sorted(listkey)[1]
            R_Index = sorted(listkey)[2]
            ll = [a, b, c, d,e,f]
            ll=sorted(ll)
            sp = ll[3]
            ep = ll[4]
            seq = sequence[sp:ep-1]
            if len(seq) >=18 and len(seq) <=22: # cap the size of random barcode
                out_dict[seq_record.id]= {"cutseq":seq, "F_Index":F_Index,"R_Index":R_Index}
    return out_dict


# Just get a fasta information for random barcode.
def randombarcode_fasta(get_randombarcode_dict):
    barhash = []
    for k in get_randombarcode_dict.keys():
        rb = get_randombarcode_dict[k]['cutseq']
        record = SeqRecord(Seq(rb), id=k ,description="", name="")
        barhash.append(record)
    SeqIO.write(barhash, "rb.fasta", "fasta")


def cluster_db(rbfasta, threads=1, cluster_id=0.9, min_seqlength=10):
    """Runs Vsearch clustering to create a FASTA file of non-redundant sequences. Selects the most abundant sequence as the centroid
    Args:
        threads (int or str):the number of processor threads to use

    Returns:
            (file): uc file with cluster information
            (file): a centroid fasta file
    """
    try:
        centroid_fasta = rbfasta.split(".")[0] + "_centroid_representative_fasta"
        uc_file = rbfasta.split(".")[0] + ".uc"
        parameters0 = ["vsearch",
                      "--cluster_size", rbfasta,
                      "--id", str(cluster_id),
                      "--sizeout", "--sizeorder","--relabel",
                      "Cluster_",
                      "--centroids", centroid_fasta,
                      "--uc", uc_file,
                      "--strand", "both",
                      "--minseqlength", str(min_seqlength),
                      "--threads", str(threads)]
        p0 = subprocess.run(parameters0, stderr=subprocess.PIPE)
        print(p0.stderr.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(str(e))
    except FileNotFoundError as f:
        print(str(f))



def mapR2clusterdb(fastafile, centroid_representative_fasta,cluster_id=0.9, min_seqlength=20):
    """Map reads and the cluster centroid information"""
    try:
        uc_map = fastafile.split(".")[0] + ".cluster_table_mapping.uc"
        readfile = fastafile
        centroidfile = centroid_representative_fasta
        parameters1 = ["vsearch",
                      "--usearch_global", readfile,
                      "--db", centroidfile,
                      "--strand", "plus",
                      "--id", str(cluster_id),
                      "--uc", uc_map,
                      "--strand", "both",
                      "--minseqlength", str(min_seqlength)]
        p1 = subprocess.run(parameters1, stderr=subprocess.PIPE)
        print(p1.stderr.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(str(e))
    except FileNotFoundError as f:
        print(str(f))

# use if NMSLIB
# index the reference list
def create_index(strings):
    index = nmslib.init(space='leven',
                            dtype=nmslib.DistType.INT,
                            data_type=nmslib.DataType.OBJECT_AS_STRING,
                            method='small_world_rand')
    index.addDataPointBatch(strings)
    index.createIndex(print_progress=True)
    return index

# get knn in bactch mode for all query
def get_knns(index, vecs):
    return zip(*index.knnQueryBatch(vecs, k=1, num_threads=4)) # zip creates a tupple with first element as knn location and second element as distance

# Display the actual strings_KNN
def display_knn(query, knn_ids_array):
    print("Query string:\n", query,"\n")
    print ("Ten nearest neighbours:\n")
    for v in ids.tolist():
        print(reflist[v])


### filter based on distance
def filter_knn_dist(qdict, mindist):
    outdict = {}
    for items in qdict.items():
         if items[1]['distance'] <= mindist:
            outdict[items[0]] = items[1]
    return outdict

def head_dict(dict, n):
    """Print first 'n ~ size' number of entries in a dictionary.

    Args:
        input (str): dict: the input dictionary to be parsed

        input (int): n: interger specifying the size of return dictionary

    Returns:
        print n entries for dictionary
    """
    k = 0
    for items in dict.items():
        if k < n:
            print("\n", items)
            k += 1


################# Run steps #################

## run mappingwith each mapping fasta information
run_vsearch("Ffasta.fasta", "merged.fasta", cluster_id=0.65) # merged fasta contain illumina reads from the experiment
run_vsearch("RCfasta.fasta", "merged.fasta", cluster_id=0.65)
run_vsearch("conserved.fasta", "merged.fasta", cluster_id=0.65)

# conctenate  all the output form vearch to make a large file
f_list = []
for filename in os.listdir():
    if filename.endswith('__output.txt'):
        f_list.append(filename)

with open('combinedfile.txt', 'w') as outfile:
    for fname in f_list:
        file = os.path.join(os.getcwd(), fname)
        with open(file) as infile:
            for line in infile:
                outfile.write(line)


# parse vearch mapping file
datplus = parse_vsearch("combinedfile.txt", strand="+")
len(datplus)
head_dict(datplus, 5)

# filter, nhits greater than three hits-- number might change dependig on the number of file your are mapping to each read.
datplus_filter = filter_vsearch(datplus, nhits=3)
len(datplus_filter)
head_dict(datplus_filter, 5)

### get mapping information for random barcodes
rbdict = get_randombarcode("merged.fasta", datplus_filter)

## get random barcode as fasta file
randombarcode_fasta(rbdict)

#Runs Vsearch clustering to create a FASTA file of non-redundant sequences.
# Selects the most abundant sequence as the centroid
cluster_db("rb.fasta", threads=8, cluster_id=0.9, min_seqlength=20)

### Combine information from cluster_db and fasta file
mapR2clusterdb("rb.fasta","rb_centroid_representative_fasta",cluster_id=0.9, min_seqlength=20)


# now filter the output from mapR2clusterdb -- output file is "rb.cluster_table_mapping.uc"
filter_uc=[]
with open("rb.cluster_table_mapping.uc", "r")as f:
    for line in f:
        ll = line.strip().split("\t")
        if ll[0]=='H':
            rb_id = ll[8]
            cluster_id = ll[9].split(";")[0]
            size = int(ll[9].split("=")[1])
            if size >= 1:
                filter_uc.append([rb_id,cluster_id,size])

### Convert list to dataframe
df_filter_uc = pd.DataFrame(filter_uc, columns=["SeqID","ClusterID","ClusterSize"])

# Convert randombacode dict to dataframe
df_rbdict = pd.DataFrame.from_dict(rbdict,orient='index')
# reindex such that SeqID is first column, which we will be using a hook to merge two data frames.
df_rbdict_ri = df_rbdict.rename_axis('SeqID').reset_index()

# Merge df_filter_uc and random barcode dataframes
df_merge = pd.merge(df_rbdict_ri, df_filter_uc, on="SeqID", how='inner')
df_merge.to_csv("df_merge.csv")

################# Processing of "df_merge.csv" is R, where we have applied some fltering criteria
# #
# library(tidyverse)
# ## input file
# data = read.csv("df_merge.csv", header=T, row.names = 1, stringsAsFactors = FALSE)
# data['rb_size'] <- nchar(data$cutseq)
#
#
# library(rio)
# data_list <- import_list("Clone_Information_byPlates.xlsx", setclass = "tbl", rbind = TRUE)
# clone <- data_list %>%
#   select(Reverse, Forward ,Clone_Number,`_file`)
#
# colnames(clone) <- c("R_Index", a"F_Index","Clone_Number","Plate_Number")
# clone$R_Index <- gsub("BarSeq_R","BarSeq_RC",clone$R_Index)
#
# library(dplyr)
# clone_data <- inner_join(data, clone, by = c("R_Index" = "R_Index", "F_Index" = "F_Index"))
# length(unique(clone_data$Clone_Number))
# sum(is.na(clone_data))
# clone_data[is.na(clone_data),]
#
#
# library(tidyverse)
# library(gtools)
# library(stringr)
#
# filtered_clone_data <- clone_data %>%
#   group_by(cutseq, Clone_Number,Plate_Number,rb_size) %>%
#   count(sort=TRUE) %>%
#   group_by(Clone_Number) %>%
#   top_n(1) %>%   # for each clone select top 1
#   filter(n>10) %>% # Remove clone if coutn < 10
#   group_by(cutseq) %>% # if same random barcode occurs at multiple clone, select random barcode with top count
#   top_n(1)
#
# filtered_clone_data["order"] <- as.numeric(gsub("Clone-","",filtered_clone_data$Clone_Number))
# filtered_clone_data_arr <- filtered_clone_data %>% arrange(order)
#
# write.csv(filtered_clone_data_arr,"filtered_clone_data_arr.csv")
#
# # Note: Eventhough we use criteria to select clone with highest count `top_n(1)`, following clone are repeated because of
# # exact count. For example, Clone-246
# filtered_clone_data_arr[duplicated(filtered_clone_data_arr$cutseq), ]
# filtered_clone_data_arr %>% filter(Clone_Number =="Clone-246")
#
# # Likewise for random barcode
# filtered_clone_data_arr[duplicated(filtered_clone_data_arr$Clone_Number), ]
# filtered_clone_data_arr %>% filter(cutseq =="CTATTCACCGACCGCGTGAT")
#
#
# ## Missing Clone
# missing_clone <- clone[!(clone$Clone_Number %in% filtered_clone_data_arr$Clone_Number), ]
#
# ## Get available information if present
# missing_clone_data <- inner_join(clone_data,missing_clone ,by = c("R_Index" = "R_Index", "F_Index" = "F_Index","Clone_Number"="Clone_Number", "Plate_Number"="Plate_Number"))
# length(unique(missing_clone_data$Clone_Number))
# sum(is.na(missing_clone_data))
# clone_data[is.na(missing_clone_data),]
#
# missing_clone_data ["order"] <- as.numeric(gsub("Clone-","",missing_clone_data$Clone_Number))
# missing_clone_file <- missing_clone_data %>%
#   group_by(F_Index,R_Index,Clone_Number, Plate_Number, order) %>%
#   count() %>%
#   arrange(order)
#
# ## Some of the clone has no sequence.
# no_info_clone <- missing_clone[!(missing_clone$Clone_Number %in% missing_clone_file$Clone_Number), ]
# no_info_clone ["order"] <- as.numeric(gsub("Clone-","",no_info_clone$Clone_Number))
# no_info_clone["n"] <- 0
# no_info_clone <- no_info_clone[, c(2,1,3,4,5,6)]
#
#
# ### Combined the clone that were removed based on filtering criteria or has no reads.
# combine_missing_clone <- bind_rows(missing_clone_file,no_info_clone) %>%
#   arrange(order)
#
# write.csv(combine_missing_clone,"missing_clone.csv")
#
#
# #########3 all vs all
# all <- read.delim("blast.out", header = FALSE)
#
#
# ## plot by plate
# combine_missing_clone %>%
#   group_by(Plate_Number) %>%
#   count()
#
# xx<- combine_missing_clone %>%
#   group_by(Clone_Number, Plate_Number) %>%
#   count()
#
#
#
# filter_df <- filtered_clone_data_arr[, c(2,3,5,6)]
#
# missing_df <- combine_missing_clone[, c(3,4,6,5)]
#
#
# all_df <- bind_rows(filter_df,missing_df)
#
# all_df %>%
#   group_by(Plate_Number) %>%
#   count()
#
#
# #
# # Remove duplicate rows of the dataframe
# all_df_distinct <- distinct(all_df)
# all_df_distinct %>%
#   group_by(Plate_Number) %>%
#   count()
#
#
# #######
# all_df_distinct_arr <- all_df_distinct %>%
#   arrange(order)
#
# ### wide table
#
# Clone_ID <- all_df_distinct_arr$Clone_Number[1:384]
# Plate_ID <- paste0("Plate-",1:12)
# readcounts_all <- all_df_distinct_arr$n
#
# mat <- matrix(readcounts_all, nrow=384, ncol=12, byrow=FALSE)
#
# mat_df <- as.data.frame(mat, row.names = Clone_ID)
# colnames(mat_df)<- Plate_ID
#
#
# install.packages("ztable")
# library(ztable)
#
# mycolor=gradientColor(low="gray",mid="orange",high="red",n=50,plot=TRUE)
#
# ztable(mat_df) %>%
#   makeHeatmap(mycolor=mycolor) %>%
#   print(caption="Table 6. Heatmap table with user-defined palette")
#
#
# ztable(mat_df) %>%
#   makeHeatmap(palette="YlOrRd") %>%
#   print(caption="Table 6. Heatmap table with user-defined palette")
#
# ztable(mat_df) %>%
#   makeHeatmap() %>%
#   print(caption="Table 1. Heatmap table showing reads per clone at each plate")
#
#
# ztable(mat_df, size=4, zebra=1, colnames.bold=TRUE, digits = 0) %>%
#   makeHeatmap() %>%
#   print(caption="Table 1. Heatmap table showing reads per clone at each plate")
#
#
#
# plate_matrix <- read.csv("Plate-matrix_count.csv", header = TRUE, row.names = 1)
#
# ztable(plate_matrix, size=4, zebra=1, colnames.bold=TRUE, digits = 0) %>%
#   makeHeatmap() %>%
#   print(caption="Table 1. Heatmap table showing reads per clone at each plate")
#

##########################################################################################3333

## Now using NMSLIB library, search for the matching entry in the referece random barcode file (provided file/information)
# Create an index with random barcodes from reference pool data
df_refinfo = pd.read_table("Phaeo_ML1.loconf.pool.txt", sep='\t', header='infer')
reflist = df_refinfo['barcode']
# initialize a new index
ref_index = create_index(list(df_refinfo['barcode']))

# Query list- use random bacode from the KEIO experiment as query.
qlist = df_merge["cutseq"].tolist()
qdict = df_merge.to_dict('index')

# query for the nearest neighbours of the first datapoint
ids, distances = ref_index.knnQuery(qlist[0], k=10)

# Display the actual strings_KNN
display_knn(qlist[0], ids)

# apply for all querys in batch
idxs,dists = get_knns(ref_index, qlist)

## load the knn and mapping location of reference back to qdict
for record in qdict.keys():
    qdict[record]['refbarloc'] =  np.asscalar(idxs[record])
    qdict[record]['distance'] =  np.asscalar(dists[record])
    qdict[record]['ref_barcode'] = reflist[np.asscalar(idxs[record])]

# filter based on KNN distance
fqdict = filter_knn_dist(qdict, 4)  ## 20% distance -- by chance will be the same == 520,000/520,000*(2^4): 520,000 is the number of reads

# convert qdict to pandas df
df_fqdict = pd.DataFrame.from_dict(fqdict,orient='index')

# inner merge information from reference and query.
df_fqdict_merge = pd.merge(df_fqdict,
                           df_refinfo,
                           left_on="ref_barcode",
                           right_on="barcode",
                           how='inner')

# Remove some columns that have same information
df_fqdict_merge_sel = df_fqdict_merge[['SeqID', 'F_Index', 'R_Index',
                                      'ClusterID', 'ClusterSize','cutseq','ref_barcode',
                                      'distance', 'rcbarcode', 'nTot', 'n', 'scaffold', 'strand',
                                       'pos', 'n2','scaffold2', 'strand2', 'pos2', 'nPastEnd']]


# save file as csv.
df_fqdict_merge_sel.to_csv("df_fqdict_merge.csv")

########### parsing genbank file
annon_df = pd.read_table("Phaeobacter_inhibens_DSM_17395_ProteinTable13044_174218.txt", sep='\t', header='infer')
df_fqdict_merge_sel_annot = pd.merge(df_fqdict_merge_sel,
                           annon_df,
                           left_on="pos",
                           right_on="Start",
                           how='inner')
                           
# save file as csv.
df_fqdict_merge_sel.to_csv("df_fqdict_merge.csv")


# Check if the position of 20bp random barcode lies with the range of start:stop position of gene.
# "filtered_clone_data_arr.csv" is output from R.
filtered_clone_data = pd.read_csv("filtered_clone_data_arr.csv",header='infer',index_col=0)
filtered_clone_data_dict= filtered_clone_data.to_dict('index')


# for filtered_clone_data_dict = get protein information from protein table file
make_list=[]
with open("Phaeobacter_inhibens_DSM_17395_ProteinTable13044_174218.txt", "r") as f:
    lines = f.readline()[1:]
    for line in f:
        ll = line.strip().split("\t")
        start = int(ll[2])
        stop = int(ll[3])
        protein_product=ll[8]
        protein_name = ll[10]
        for record in filtered_clone_data_dict.keys():
            rb_position = filtered_clone_data_dict[record]['pos']
            if start <= rb_position <=stop:
                make_list.append([protein_product,start, stop, protein_name,rb_position])


# check if the key is in the dict, if yes append else add.
def add_or_append(dictionary, key, value):
    if key not in dictionary:
        dictionary[key] = []
    dictionary[key].append(value)

# check if the key is in the dict, if yes append else add.
sdict={}
with open("Phaeobacter_inhibens_DSM_17395_ProteinTable13044_174218.txt", "r") as f:
    lines = f.readline()[1:]
    for line in f:
        ll = line.strip().split("\t")
        start = int(ll[2])
        stop = int(ll[3])
        locus_tag = ll[7]
        protein_product=ll[8]
        protein_name = ll[10]
        for record in filtered_clone_data_dict.keys():
            rb_position = filtered_clone_data_dict[record]['pos']
            to_add = (protein_product,locus_tag,start,stop,protein_name,rb_position)
            if start <= rb_position <=stop:
                add_or_append(sdict, record, to_add)

##################
# combined two dict based on common keys
d1 = filtered_clone_data_dict
d2 = sdict
combined_dict = defaultdict(list)
for d in (d1, d2): # you can list as many input dicts as you want here
    for key, value in d.items():
        combined_dict[key].append(value)

####
# From the combined dictionary, select only the information that we want in the final output file
final_list=[]
for key in combined_dict.keys():
        cutseq = combined_dict[key][0]['cutseq']
        ref_barcode = combined_dict[key][0]['ref_barcode']
        distance = combined_dict[key][0]['distance']
        cutseq_size = combined_dict[key][0]['cutseq_size']
        Plate_Number = combined_dict[key][0]['Plate_Number']
        Clone_Number = combined_dict[key][0]['Clone_Number']
        scaffold = combined_dict[key][0]['scaffold']
        strand = combined_dict[key][0]['strand']
        pos = combined_dict[key][0]['pos']
        reads = combined_dict[key][0]['n']
        if len(combined_dict[key]) == 1:
            final_list.append([cutseq,ref_barcode,distance,cutseq_size,Plate_Number,
                               Clone_Number,scaffold,strand,pos,reads])
        elif len(combined_dict[key]) == 2:
            gene_info = combined_dict[key][1]
            final_list.append([cutseq,ref_barcode,distance,cutseq_size,Plate_Number,
                               Clone_Number,scaffold,strand,pos,reads,gene_info])

## Convert list to dataframe
df_final_list = pd.DataFrame(final_list,
                            columns=["cutseq","ref_barcode","distance","cutseq_size",
                                     "Plate_Number","Clone_Number","scaffold",
                                     "strand","pos","reads","gene_info-[protein_name,locus_tag start, stop, protein_product, rb20_position]"])

# save dataframe as csv
df_final_list.to_csv("df_mapping_with_geneinfo.csv")
