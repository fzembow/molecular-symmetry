#!/usr/bin/python2.6

#
# a Hadoop streaming mapper that downloads a data shard 
# from NCBI and runs a molecular calculation on it
# 

import sys, subprocess, json
import extract

for line in sys.stdin:
    
    #get the URLS for this task
    data = line.strip().split('\t')
    urls = json.loads(data[1])
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
    subprocess.call(["rm","-f","data.sdf"])
    if has3D:
        subprocess.call(["rm","-f","data-3d.sdf"])
