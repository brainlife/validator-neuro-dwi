#!/usr/bin/env python

import os
import json
import re
import subprocess
import nibabel
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table

# Things that this script checks
# 
# * make sure nibabel runs successfully on specified dwi file
# * make sure dwi is 4d
# * raise warning if dwi transformation matrix isn't unit matrix (identity matrix)
# * make sure bvecs and bvals can be read
# * make sure bvecs's cols count matches dwi's 4th dimension number
# * make sure bvecs has 3 rows
# * make sure bvals's cols count matches dwi's 4th dimension number
# * make sure bvals has 1 row

# display where this is running
import socket

print(socket.gethostname())

with open('config.json') as config_json:
    config = json.load(config_json)

results = {"errors": [], "warnings": []}
directions = None
gtab = []

def check_affine(affine):
    if affine[0][0] != 1: results['warnings'].append("transform matrix 0.1 is not 1")
    if affine[0][1] != 0: results['warnings'].append("transform matrix 0.2 is not 0")
    if affine[0][2] != 0: results['warnings'].append("transform matrix 0.2 is not 0")
    if affine[1][0] != 0: results['warnings'].append("transform matrix 1.0 is not 0")
    if affine[1][1] != 1: results['warnings'].append("transform matrix 1.1 is not 1")
    if affine[1][2] != 0: results['warnings'].append("transform matrix 1.2 is non 0")
    if affine[2][0] != 0: results['warnings'].append("transform matrix 2.0 is not 0")
    if affine[2][1] != 0: results['warnings'].append("transform matrix 2.1 is not 0")
    if affine[2][2] != 1: results['warnings'].append("transform  matrix 2.2 is not 1")

#def isInt(v):
#    v = v.strip()
#    return v=='0' or (v if v.find('..') > -1 else v.lstrip('-+').rstrip('0').rstrip('.')).isdigit()

def isFloat(v):
    try:     i = float(v)
    except:  return False
    return True

def isInt(v):
    try:     i = int(v)
    except:  return False
    return True

print("checking input paramerters")
if config['dwi'] is None:
    results['errors'].append("dwi not set")

if config['bvecs'] is None:
    results['errors'].append("bvecs not set")

if config['bvals'] is None:
    results['errors'].append("bvals not set")

if len(results['errors']) == 0:
    try:
        print("validating bvecs")
        bvecs = open(config['bvecs'])
        bvecs_rows = bvecs.readlines()
        bvecs_cols = bvecs_rows[0].strip().replace(",", " ")

        # remove double spaces
        bvecs_cols_clean = re.sub(' +', ' ', bvecs_cols)
        bvecs_cols = bvecs_cols_clean.split(' ')
        directions = len(bvecs_cols)

        # check 4d size
        if directions < 10:
            results['warnings'].append("direction seems too small.." + str(directions))

        # bvecs should have 3 row
        if len(bvecs_rows) != 3:
            results['errors'].append("bvecs should have 3 rows but it has " + str(len(bvecs_rows)))

        # write out bvecs (write in FSL format - space delimited)
        f = open('dwi.bvecs', 'w')
        for row in bvecs_rows:
            row = row.strip().replace(",", " ")
            row_clean = re.sub(' +', ' ', row)

            f.write(row_clean)
            f.write("\n")

            #check to make sure all values are float
            for v in row_clean.split(' '):
                if not isFloat(v):
                    results['errors'].append("bvecs contains non float:" +v)

        f.close();

    except IOError:
        print("failed to load bvecs:" + config['bvecs'])
        results['errors'].append("Couldn't read bvecs")

    try:
        print("validating bvals")
        bvals = open(config['bvals'])
        bvals_rows = bvals.readlines()
        bvals_cols = bvals_rows[0].strip().replace(",", " ")

        # remove double spaces
        bvals_cols_clean = re.sub(' +', ' ', bvals_cols)
        bvals_cols = bvals_cols_clean.split(" ")

        if directions != len(bvals_cols):
            results['errors'].append(
                "bvals column count which is " + str(len(bvals_cols)) + " doesn't match dwi's 4d number:" + str(
                    directions))

        if len(bvals_rows) != 1:
            results['errors'].append("bvals should have 1 row but it has " + str(len(bvals_rows)))

        results['datatype_tags'] = []

        # analyze single / multi shell (0, 2000 is single. 0,2000,3000 is multi shell)
        unique_bvals = list(set(bvals_cols))
        # results['tags'] = ['b'+bval for bval in unique_bvals]

        bvalues = list(set(bvals_cols))

        # is normalized?
        normalized = True
        for bvalue in bvalues:
            if not isInt(bvalue):
                results['errors'].append("bvals contains non int:" +v)
            else:
                if float(bvalue) != round(float(bvalue), -2):
                    normalized = False
        if normalized:
            results['datatype_tags'].append("normalized")

        # is single shell?
        if len(bvalues) <= 2:
            results['datatype_tags'].append("single_shell")
            results['tags'] = ["b" + str(bvalues[1])];

        # write out normalized bvals (write in FSL format - space delimited)
        f = open('dwi.bvals', 'w')
        f.write(bvals_cols_clean)
        f.write("\n")
        f.close()

    except IOError:
        print("failed to load bvals:" + config['bvals'])
        results['errors'].append("Couldn't read bvals")

    print("analyzing bvecs/bvals")
    bvals, bvecs = read_bvals_bvecs('dwi.bvals', 'dwi.bvecs')
    gtab = gradient_table(bvals, bvecs)
    results['gtab_info'] = gtab.info

    #sort into shells (1000th)
    shells = {}
    for i in range(len(gtab.bvals)):
        bval = gtab.bvals[i]
        bvec = gtab.bvecs[i]
        shell = str(round(bval, -3))
        if shell not in shells:
            shells[shell] = []
        shells[shell].append((i, bval, bvec*bval))

    #output into plotly format
    data = []
    for shell in shells:
        xs = []
        ys = []
        zs = []
        texts = []
        for v in shells[shell]:
            texts.append(v[0])
            xs.append(v[2][0])
            ys.append(v[2][1])
            zs.append(v[2][2])

        if shell == "0.0":
            color = "black"
        elif shell == "1000.0":
            color = "blue"
        elif shell == "2000.0":
            color = "green"
        elif shell == "3000.0":
            color = "purple"
        elif shell == "4000.0":
            color = "cyan"
        else:
            color = "red"

        data.append({
            'type': 'scatter3d',
            'mode': 'text', 
            'name': shell, 
            'x': xs,
            'y': ys,
            'z': zs,
            'text': texts,
            'textfont': {
                'color': color,
                'size': 8
            }
        })
    results['brainlife'] = [{
        'type': 'plotly',
        'name': 'Gradient Table(bvecs/bvals)',
        'layout': {},
        'data': data,
    }]

    try:
        print("validating dwi")
        img = nibabel.load(config['dwi'])
        results['dwi_headers'] = str(img.header)
        results['dwi_base_affine'] = str(img.header.get_base_affine())

        # check dimensions
        dims = img.header['dim'][0]
        if dims != 4:
            results['errors'].append("DWI should be 4D but has " + str(dims))

        # 4d should have same col num as bvecs/bval directions
        if directions != img.header['dim'][4]:
            results['errors'].append(
                "bvecs column count which is " + str(img.header['dim'][4]) + " doesn't match bvecs col counts:" + str(
                    directions))

        # check affine
        check_affine(img.header.get_base_affine())

        #create symlink
        os.symlink(config['dwi'], "dwi.nii.gz")

    except Exception as e:
        results['errors'].append("nibabel failed on dwi. error code: " + str(e))

with open("product.json", "w") as fp:
    json.dump(results, fp)

print(results)
