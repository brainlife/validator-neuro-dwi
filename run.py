#!/usr/bin/env python

import os
import json
import re
import subprocess

# Things that this script checks
# 
# * make sure mrinfo runs successfully on specified dwi file
# * make sure dwi is 4d
# * raise warning if dwi transformation matrix isn't unit matrix (identity matrix)
# * make sure bvecs and bvals can be read
# * make sure bvecs's cols count matches dwi's 4th dimension number
# * make sure bvecs has 3 rows
# * make sure bvals's cols count matches dwi's 4th dimension number
# * make sure bvals has 1 row

#display where this is running
import socket
print(socket.gethostname())

with open('config.json') as config_json:
    config = json.load(config_json)

results = {"errors": [], "warnings": []}

directions = None

#check dwi
if config['dwi'] is None: 
    results['errors'].append("dwi not set")
else:
    try: 
        print("running dwi mrinfo")
        info = subprocess.check_output(["mrinfo", config['dwi']], shell=False)
        results['detail'] = info
        info_lines = info.split('\n')

        #check dimentions
        dim=info_lines[4]
        dims=dim.split("x")
        if len(dims) != 4:
            results['errors'].append("DWI file specified doesn't have 4 dimentions")
        else:
            directions = int(dims[3])
            if directions < 10: 
                results['errors'].append("DWI's 4D seems too small",directions)

        #check transform
        tl = info_lines[-5:-1] #grab last 4 lines (minus very last which is newline)
        tl[0] = tl[0][12:] #remove "Transform:"
        m = []
        for line in tl:
            line = line.strip()
            strs = re.split('\s+', line)
            ml = []
            for s in strs:
                ml.append(float(s))   
            m.append(ml)

        if m[0][0] == 1 and m[0][1] == 0 and m[0][1] == 0 and m[0][0] == 1 and m[0][1] == 0 and m[0][1] == 0 and m[0][0] == 1 and m[0][1] == 0 and m[0][1] == 0:  
            None #good
        else:
            results['warnings'].append("DWI has non-optimal transformation matrix. It should be 1 0 0 / 0 1 0 / 0 0 1")
        
        #TODO - normalize (for now, let's just symlink)
        os.symlink(config['dwi'], "dwi.nii.gz")
        
    except subprocess.CalledProcessError as err:
        results['errors'].append("mrinfo failed on dwi. error code: "+str(err.returncode))

def readb(name):
    #read bvecs
    b = open(name)
    rows = b.readlines()
    normalized = [];
    for row in rows:
        row = row.strip().replace(",", " ").replace("\t", " ") #make it space delimited
        row = re.sub(' +', ' ', row) #remove double spaces
        normalized.append(row.split(' '))
    b.close()
    return normalized

#load bvecs (and check size)
if config['bvecs'] is None: 
    results['errors'].append("bvecs not set")
else:
    try: 
        bvecs = readb(config['bvecs'])
        bvecs_cols = bvecs[0]

        if directions:
            if directions != len(bvecs_cols):
                results['errors'].append("bvecs column count which is "+str(len(bvecs_cols))+" doesn't match dwi's 4d number:"+str(directions))

        if len(bvecs) != 3:
            results['errors'].append("bvecs should have 3 rows but it has "+str(len(bvecs)))

        #write out normalized bvecs (write in FSL format - space delimited)
        f = open('dwi.bvecs', 'w')
        for row in bvecs:
            f.write(" ".join(row))
            f.write("\n")

    except IOError:
        print("failed to load bvecs:"+config['bvecs'])
        results['errors'].append("Couldn't read bvecs")

#load bvals (and check cols)
if config['bvals'] is None: 
    results['errors'].append("bvals not set")
else:
    try: 
        bvals = readb(config['bvals'])
        bvals_cols = bvals[0]

        if directions:
            if directions != len(bvals_cols):
                results['errors'].append("bvals column count which is "+str(len(bvals_cols))+" doesn't match dwi's 4d number:"+str(directions))

        if len(bvals) != 1:
            results['errors'].append("bvals should have 1 row but it has "+str(len(bvals)))

        results['datatype_tags'] = []

        #analyze single / multi shell (0, 2000 is single. 0,2000,3000 is multi shell)
        unique_bvals = list(set(bvals_cols))
        results['tags'] = ['b'+bval for bval in unique_bvals]

        bvalues = list(set(bvals_cols))

        #is normalized?
        normalized = True
        for bvalue in bvalues:
            if float(bvalue) != round(float(bvalue), -2):
                normalized = False
        if normalized:
            results['datatype_tags'].append("normalized")

        #is single shell?
        if len(bvalues) <= 2:
            results['datatype_tags'].append("single_shell")
            results['tags'] = ["b"+str(bvalues[1])];
        
        #write out normalized bvals (write in FSL format - space delimited)
        f = open('dwi.bvals', 'w')
        for row in bvals:
            f.write(" ".join(row))
            f.write("\n")

    except IOError:
        print("failed to load bvals:"+config['bvals'])
        results['errors'].append("Couldn't read bvals")

#deprecated
with open("products.json", "w") as fp:
    json.dump([results], fp)

with open("product.json", "w") as fp:
    json.dump(results, fp)

