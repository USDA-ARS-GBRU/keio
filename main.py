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


def run_vsearch(maping_fasta, reads.fasta, cluster_id=0.75):
    """ Returns mapping information.
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
        parameters = ["vsearch", "--usearch_global", reads.fasta, "--db", maping_fasta,
                      "--id", str(cluster_id), "--userfield", out_info, "--strand", "plus", "--userout", outputfile]
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

def filter_vsearch(sdict, nhits):
    outdict = {}
    for items in sdict.items():
        if len(items[1]) >= nhits:
            outdict[items[0]] = items[1]
    return outdict


def get_randombarcode(keio.fasta, filter_vsearch_dict):
    out_dict={}
    for seq_record in SeqIO.parse(keio.fasta, "fasta"):
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


def randombarcode_fasta(get_randombarcode_dict):
    barhash = []
    for k in out_dict.keys():
        rb = out_dict[k]['cutseq']
        record = SeqRecord(Seq(rb), id=k ,description="", name="")
        barhash.append(record)
    SeqIO.write(barhash, "rb.fasta", "fasta")


### apply to actual data
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
for filename in os.listdir('queryBarcode'):
    if filename.endswith('__output.txt'):
        f_list.append(filename)

with open('combinedfile.txt', 'w') as outfile:
    for fname in f_list:
        file = os.path.join(os.getcwd(), "queryBarcode", fname)
        with open(file) as infile:
            for line in infile:
                outfile.write(line)


datplus = parse_vsearch("combinedfile.txt", strand="+")
len(datplus)
head_dict(datplus, 5)

# filter, nhits greater than
datplus_filter = filter_vsearch(datplus, nhits=3)
len(datplus_filter)
head_dict(datplus_filter, 5)

### get mapping information for random barcodes
rbdict = get_randombarcode("merged.fasta", datplus_filter)

## get random barcode as fasta file
randombarcode_fasta(rbdict)

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
ids, distances = ref_index.knnQuery(qlist[0], k=1)

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
fqdict = filter_knn_dist(qdict, 4)  ## 20% distance -- by chance will be the same == 520,000/520,000*(2^4)

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
# Check if the position of 20bp random barcode lies with the range of start:stop position of gene.
filtered_clone_data = pd.read_csv("filtered_clone_data_arr.csv",header='infer',index_col=0)
filtered_clone_data_dict= filtered_clone_data.to_dict('index')

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
                print(protein_product,start, stop, protein_name,rb_position)
                make_list.append([protein_product,start, stop, protein_name,rb_position])

def add_or_append(dictionary, key, value):
    if key not in dictionary:
        dictionary[key] = []
    dictionary[key].append(value)


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


d1 = filtered_clone_data_dict
d2 = sdict

combined_dict = defaultdict(list)

for d in (d1, d2): # you can list as many input dicts as you want here
    for key, value in d.items():
        combined_dict[key].append(value)


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


df_final_list = pd.DataFrame(final_list,
                            columns=["cutseq","ref_barcode","distance","cutseq_size",
                                     "Plate_Number","Clone_Number","scaffold",
                                     "strand","pos","reads","gene_info-[protein_name,locus_tag start, stop, protein_product, rb20_position]"])

df_final_list.to_csv("df_mapping_with_geneinfo.csv")
