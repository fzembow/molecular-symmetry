#!/usr/bin/env python

#
# a Hadoop streaming mapper that downloads a data shard 
# from NCBI and runs a molecular calculation on it
# 

TWO_D_FOLDER = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
THREE_D_FOLDER = "ftp://ftp.ncbi.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/"

import sys, subprocess

import extract

for line in sys.stdin:
    
    two_d_url = TWO_D_FOLDER + line.strip()
    three_d_url = THREE_D_FOLDER + line.strip()
    
    #download the archives
    subprocess.call(["wget","-O","test.sdf.gz", "-q", two_d_url])
    
    #uncompress the archives
    subprocess.call(["gunzip","test.sdf.gz"])
    
    #process the symmetries
    extract.process_multiple("test.sdf")
    
    #remove the archive to free space
    subprocess.call(["rm","test.sdf"])