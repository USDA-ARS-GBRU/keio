#!/usr/bin/env python
"""
keio: A python module to process illumina reads for keio-collection type project.

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
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.SeqRecord import SeqRecord



def run_vsearch(barcodefile, reads, cluster_id=0.75):
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
        outputfile = barcodefile.split(".")[0] + "__output.txt"
        parameters = ["vsearch", "--usearch_global", reads, "--db", barcodefile,
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





#########step to run the program ############
