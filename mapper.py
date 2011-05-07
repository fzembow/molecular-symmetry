#!/usr/bin/env python

#
# a Hadoop streaming mapper that downloads a data shard 
# from NCBI and runs a molecular calculation on it
# 

TWO_D_FOLDER = "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
THREE_D_FOLDER = "ftp://ftp.ncbi.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/"

import sys, subprocess
import extract

# f = open("urls-output","w")
# seen = False

for line in sys.stdin:
    
    # compound = line.strip()
    # 
    # compound2 = compound[10:].replace("_0","_")
    # 
    # if not seen:
    #     f.write(TWO_D_FOLDER + compound + "\t" + THREE_D_FOLDER + compound2 + "\n")
    # else:
    #     f.write(TWO_D_FOLDER + compound + "\n")
    # 
    # if compound2.startswith("49825001"):
    #     seen = True
    # 
    # continue
    
    urls = line.strip().split("\t")
    if len(urls) == 2:
        has3D = True
    else:
        has3D = False
    
    #download the archives
    subprocess.call(["wget","-O","data.sdf.gz", "-q", urls[0]])
    if has3D:
        subprocess.call(["wget","-O","data-3d.sdf.gz", "-q", urls[1]])
    
    #uncompress the archives
    subprocess.call(["gunzip","data.sdf.gz"])
    if has3D:
        subprocess.call(["gunzip","data-3d.sdf.gz"])
    
    #process the symmetries
    if has3D:
        extract.process_multiple("data.sdf", "data-3d.sdf")
    else:
        extract.process_multiple("data.sdf")
    
    #remove the archive to free space
    subprocess.call(["rm","data.sdf"])
    if has3D:
        subprocess.call(["rm","data-3d.sdf"])