#!/usr/bin/python2.6

#
# a Hadoop streaming mapper that downloads a data shard 
# from NCBI and runs a molecular calculation on it
# 

import sys, json
import extract

count = 0

for line in sys.stdin:

    count += 1
    
    urls = line.strip().split("\t")
    if len(urls) == 2:
        has3D = True
    else:
        has3D = False

    if has3D:
        sys.stdout.write(str(count) + "\t" + json.dumps([urls[0], urls[1]]) + "\n")
    else:
        sys.stdout.write(str(count) + "\t" + json.dumps([urls[0]]) + "\n")