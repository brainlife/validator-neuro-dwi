#!/bin/env python

import os
import json
import re
import subprocess

# Things that this script checks
# 
# (for t1)
# * make sure mrinfo runs successfully on specified t1 file
# * make sure t1 is 3d
# * raise warning if t1 transformation matrix isn't unit matrix (identity matrix)
#
# (for dwi)
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

#TODO - how should I keep up wit this path?
mrinfo="/N/soft/rhel6/mrtrix/0.3.15/mrtrix3/release/bin/mrinfo"

directions = None

if 'dwi' in config:

    #check dwi
    if config['dwi'] is None: 
        results['errors'].append("dwi not set")
    else:
        try: 
            print "running dwi mrinfo"
            info = subprocess.check_output([mrinfo, config['dwi']], shell=False)
            results['detail'] = info
            info_lines = info.split("\n")

            #check dimentions
            dim=info_lines[3]
            dims=dim.split("x")
            if len(dims) != 4:
                results['errors'].append("DWI file specified doesn't have 4 dimentions")
            else:
                directions = int(dims[3])
                if directions < 10: 
                    results['errors'].append("DWI's 4D seems too small",directions)

            #check transform
            tl = info_lines[9:12]
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
            
        except subprocess.CalledProcessError as err:
            results['errors'].append("mrinfo failed on dwi. error code: "+str(err.returncode))

    #load bvecs (and check size)
    if config['bvecs'] is None: 
        results['errors'].append("bvecs not set")
    else:
        try: 
            bvecs = open(config['bvecs'])
            bvecs_rows = bvecs.readlines()
            bvecs_cols = bvecs_rows[0].strip().replace(",", " ")

            #remove double spaces
            bvecs_cols_clean = re.sub(' +', ' ', bvecs_cols) 
            bvecs_cols = bvecs_cols_clean.split(' ')		

            if directions:
                if directions != len(bvecs_cols):
                    results['errors'].append("bvecs column count which is "+str(len(bvecs_cols))+" doesn't match dwi's 4d number:"+str(directions))

            if len(bvecs_rows) != 3:
                results['errors'].append("bvecs should have 3 rows but it has "+str(len(bvecs_rows)))

        except IOError:
            print "failed to load bvecs:"+config['bvecs']
            results['errors'].append("Couldn't read bvecs")

    #load bvals (and check cols)
    if config['bvals'] is None: 
        results['errors'].append("bvals not set")
    else:
        try: 
            bvals = open(config['bvals'])
            bvals_rows = bvals.readlines()
            bvals_cols = bvals_rows[0].strip().replace(",", " ")

            #remove double spaces
            bvals_cols_clean = re.sub(' +', ' ', bvals_cols) 
            bvals_cols = bvals_cols_clean.split(" ") 

            if directions:
                if directions != len(bvals_cols):
                    results['errors'].append("bvals column count which is "+str(len(bvals_cols))+" doesn't match dwi's 4d number:"+str(directions))

            if len(bvals_rows) != 1:
                results['errors'].append("bvals should have 1 row but it has "+str(len(bvals_rows)))

            #analyze single / multi shell (0, 2000 is single. 0,2000,3000 is multi shell)
            #TODO - I think I should round bvalue?
            results['bvalues'] = list(set(bvals_cols))
            if len(results['bvalues']) <= 2:
                results['tags'] = ["single_shell"]

        except IOError:
            print "failed to load bvals:"+config['bvals']
            results['errors'].append("Couldn't read bvals")

with open("products.json", "w") as fp:
    json.dump([results], fp)

