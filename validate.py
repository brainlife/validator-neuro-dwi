#!/usr/bin/env python3

import os
import json
import re
import subprocess
import nibabel
import shutil
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from PIL import Image, ImageDraw
import binascii

import math
import numpy as np

from scipy.ndimage import zoom

from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

# display where this is running (debug?)
import socket
print(socket.gethostname())

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
# * check for bvecs flipping 

#Returns the unit vector of the vector. 
def unit_vector(vector):
    return vector / np.linalg.norm(vector)

#Returns the angle in radians between vectors 'v1' and 'v2' 
def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

#flip angle that's >90 to face the same direction
def flip_angle(a):
    if a > math.pi/2:
        return math.pi - a
    return a

#find the most common bvals used
def most_common(bvals):
    round_bvals = []
    for bval in bvals:
        round_bvals.append(round(bval, -2))
    return max(set(round_bvals), key=bvals.count)

#the heart of bvecs/bvals detection
def sum_diag(img, shift):
    sum=img[0]
    for i in range(1, img.shape[0]):
        sum = np.roll(sum, shift)
        sum = np.add(sum, img[i])
    return sum

with open('config.json') as config_json:
    config = json.load(config_json)

results = {"errors": [], "warnings": [], "brainlife": [], "datatype_tags": [], "tags": []}
#directions = None
gtab = None

def warning(msg):
    global results
    results['warnings'].append(msg) 
    #results['brainlife'].append({"type": "warning", "msg": msg}) 
    print(msg)

def error(msg):
    global results
    results['errors'].append(msg) 
    #results['brainlife'].append({"type": "error", "msg": msg}) 
    print(msg)

def check_affine(affine):
    if abs(affine[0][1]) > 0.01: warning("transform matrix 0.2 is not near 0")
    if abs(affine[0][2]) > 0.01: warning("transform matrix 0.2 is not near 0")
    if abs(affine[1][0]) > 0.01: warning("transform matrix 1.0 is not near 0")
    if abs(affine[1][2]) > 0.01: warning("transform matrix 1.2 is non near 0")
    if abs(affine[2][0]) > 0.01: warning("transform matrix 2.0 is not near 0")
    if abs(affine[2][1]) > 0.01: warning("transform matrix 2.1 is not near 0")

def isFloat(v):
    try:     i = float(v)
    except:  return False
    return True

def get_change(current, previous):
    if current == previous:
        return 100.0
    try:
        return (abs(current - previous) / previous) * 100.0
    except ZeroDivisionError:
        return 0

def fix_level(image):
    image = image - image.mean(axis=0)
    image_max = np.max(image)
    return (image / image_max)*300

def generate_sample(slice1, slice2):
    slice1 = fix_level(slice1).T
    slice2 = fix_level(slice2).T

    pos = np.subtract(slice1, slice2).clip(min=0)*2
    neg = np.subtract(slice2, slice1).clip(min=0)*2

    image1_and_pos = Image.fromarray(slice1 + pos*2).convert('L')
    image2_and_neg = Image.fromarray(slice2 + neg*2).convert('L')

    image1 = Image.fromarray(slice1).convert('L')
    image2 = Image.fromarray(slice2).convert('L')

    image1_annotated = Image.merge('RGB', (image1_and_pos, image1, image1))
    image2_annotated = Image.merge('RGB', (image2_and_neg, image2, image2))

    draw1 = ImageDraw.Draw(image1_annotated)
    draw1.text((10, 10), "Should see features in /") 
    draw2 = ImageDraw.Draw(image2_annotated)
    draw2.text((10, 10), "Should see features in \\") 

    return (image1_annotated.convert('RGB'), image2_annotated.convert('RGB'))

print("checking input paramerters")
if config['dwi'] is None:
    error("dwi not set")

if config['bvecs'] is None:
    error("bvecs not set")

if config['bvals'] is None:
    error("bvals not set")

with open(config['dwi'], 'rb') as test_f:
    if binascii.hexlify(test_f.read(2)) != b'1f8b':
        error("dwi doesn't look like gzipped");

if len(results['errors']) == 0:
    print("loading bvecs/bvals with dipy.io")
    bvals, bvecs = read_bvals_bvecs(config["bvals"], config["bvecs"])
    try: 
        gtab = gradient_table(bvals, bvecs)
    except ValueError as err:
        warning(str(err))

        #re-try with rediculous atol to bypass the check (some data has [0,0,0] vector!
        print("... failed to parse bvals/bvecs.. trying with atol=1")
        gtab = gradient_table(bvals, bvecs, atol=1)

    # write out bvecs (write in FSL format - space delimited)
    x = ""
    y = ""
    z = ""
    for bvec in bvecs:
        if x != "":
            x += " "
        if y != "":
            y += " "
        if z != "":
            z += " "
        x += str(bvec[0])
        y += str(bvec[1])
        z += str(bvec[2])

    if not os.path.exists("output"):
        os.mkdir("output")

    if not os.path.exists("secondary"):
        os.mkdir("secondary")

    f = open('output/dwi.bvecs', 'w')
    f.write(x+"\n")
    f.write(y+"\n")
    f.write(z+"\n")
    f.close()

    #round bvals
    bvals = [int(round(x)) for x in bvals]

    f = open('output/dwi.bvals', 'w')
    f.write(" ".join(str(x) for x in bvals))
    f.close()

    #reload bvals/bvecs just written to make sure it's still good (mainly for debugging)
    bvals_check, bvecs_check = read_bvals_bvecs("output/dwi.bvals", "output/dwi.bvecs")

    #copy bvecs/bvals to secondary
    #shutil.copyfile("output/dwi.bvals", "secondary/dwi.bvals")
    #shutil.copyfile("output/dwi.bvecs", "secondary/dwi.bvecs")

    # analyze single / multi shell (0, 2000 is single. 0,2000,3000 is multi shell)
    bvalues = list(set(bvals))
    bvalues.sort()

    # is normalized?
    normalized = True
    for bvalue in bvalues:
        if float(bvalue) != round(float(bvalue), -2):
            normalized = False
    if normalized:
        results['datatype_tags'].append("normalized")

    # is single shell?
    if len(bvalues) == 2:
        results['datatype_tags'].append("single_shell")
        results['tags'].append("b" + str(bvalues[1]))

    if len(bvalues) == 1 and bvalues[0] == "0":
        results['tags'].append("b0")

    try:
        print("loading dwi to check bvecs flipping")
        img = nibabel.load(config['dwi'])

        #deprecated.. does bids export still use this?
        #results['dwi_headers'] = str(img.header) #need to str() so that we can save it to product.json
        #results['dwi_affine'] = str(img.affine) #need to str() as array is not serializable

        results['meta'] = {'nifti_headers': {}}
        for key in img.header:
            value = img.header[key]
            results['meta']['nifti_headers'][key] = value
        results['meta']['nifti_headers']['base_affine'] = img.header.get_base_affine()

        results["bvals"] = bvals
        results["bvecs"] = bvecs

        dimX = img.header["pixdim"][1]
        dimY = img.header["pixdim"][2]
        dimZ = img.header["pixdim"][3]
        if abs(dimX - dimY) > dimX*0.1 or abs(dimX - dimZ) > dimX*0.1 or abs(dimY - dimZ) > dimX*0.1:
            warning("pixdim is not close to isomorphic.. some dwi processing might fail")

        #results['meta'] = {
        #    "dim":img.header['dim'].tolist(),
        #    "pixdim":img.header['pixdim'].tolist()
        #}
        
        #determine storage orientation
        #http://community.mrtrix.org/t/mrconvert-flips-gradients/581/6
        det = np.linalg.det(img.affine)
        results["meta"]['dwi_affine_determinant'] = det
        radiological=False
        print("affine determinant", det)
        if det < 0:
            radiological=True
            results["meta"]['storage_orientation'] = 'radiological'
            print('storage_orientation: radiological - matches bvecs orientation')
        else:
            results["meta"]['storage_orientation'] = 'neurological'
            print('storage_orientation: neurological - flipping x')
            warning("storage orientation is neurologial (det>0). Watch out!")
            results['tags'].append("neurological")

        # check dimensions
        dims = img.header['dim'][0]
        if dims != 4:
            error("DWI should be 4D but has " + str(dims))

        # check affine
        check_affine(img.header.get_base_affine())

        #create symlink
        if(os.path.islink('output/dwi.nii.gz')):
            os.remove('output/dwi.nii.gz')
        os.symlink("../"+config['dwi'], "output/dwi.nii.gz")

    except Exception as err:
        error("nibabel failed on dwi. error code: " + str(err))

    ###############################################################################################
    #
    # check bvecs flipping
    #
    #find the most common bvals (most likely to find the right directions)
    #TODO if if there are near identical number of bvalues, should I use higher bvalue?
    b=most_common(bvals)
    print("bvalue for flip check: %d" % b)
    
    #calculate bvecs angle from various reference angles
    angs = []
    for idx in range(len(bvecs)):
        bvec = bvecs[idx]
        bval = bvals[idx]

        #ignore bvecs with low bval
        if bval < 500:
            continue

        #ignore bvecs that's too off
        if abs(bval - b) > 300:
            continue

        #ignore vec like [0,0,0] with non-0 bval? maybe it means b0?
        if np.linalg.norm(bvec) == 0:
            continue

        #calculate angle from x/y/z directions
        x1_ang = flip_angle(angle_between(bvec, (1,1,0)))
        x2_ang = flip_angle(angle_between(bvec, (-1,1,0)))
        y1_ang = flip_angle(angle_between(bvec, (0,1,1)))
        y2_ang = flip_angle(angle_between(bvec, (0,-1,1)))
        z1_ang = flip_angle(angle_between(bvec, (1,0,1)))
        z2_ang = flip_angle(angle_between(bvec, (1,0,-1)))
        angs.append((x1_ang, x2_ang, y1_ang, y2_ang, z1_ang, z2_ang, bvec, bval, idx));

    if len(angs) < 10: #TODO - not sure how small is too small
        warning("we don't have enough directions to run flip test.. skipping")
    else:
        print("loading data for flip test")
        #img_data = img.get_fdata()
        #print(nibabel.is_proxy(img.dataobj))
        #img_data = np.asanyarray(img.dataobj)
        
        #pull slices for each tests (pick 3 closest to the axis and sum them)
        angs.sort(key=lambda tup: tup[0])
        print("pulling slices for x1: %d(%f) %d(%f) %d(%f)" % (angs[0][8], angs[0][0], angs[1][8], angs[1][0], angs[2][8], angs[2][0]))
        vol_x1 = img.dataobj[..., angs[0][8]] + img.dataobj[..., angs[1][8]] + img.dataobj[..., angs[2][8]]

        angs.sort(key=lambda tup: tup[1])
        print("pulling slices for x2: %d(%f) %d(%f) %d(%f)" % (angs[0][8], angs[0][1], angs[1][8], angs[1][1], angs[2][8], angs[2][1]))
        vol_x2 = img.dataobj[..., angs[0][8]] + img.dataobj[..., angs[1][8]] + img.dataobj[..., angs[2][8]]

        angs.sort(key=lambda tup: tup[2])
        print("pulling slices for y1: %d(%f) %d(%f) %d(%f)" % (angs[0][8], angs[0][2], angs[1][8], angs[1][2], angs[2][8], angs[2][2]))
        vol_y1 = img.dataobj[..., angs[0][8]] + img.dataobj[..., angs[1][8]] + img.dataobj[..., angs[2][8]]

        angs.sort(key=lambda tup: tup[3])
        print("pulling slices for y2: %d(%f) %d(%f) %d(%f)" % (angs[0][8], angs[0][3], angs[1][8], angs[1][3], angs[2][8], angs[2][3]))
        vol_y2 = img.dataobj[..., angs[0][8]] + img.dataobj[..., angs[1][8]] + img.dataobj[..., angs[2][8]]

        angs.sort(key=lambda tup: tup[4])
        print("pulling slices for z1: %d(%f) %d(%f) %d(%f)" % (angs[0][8], angs[0][4], angs[1][8], angs[1][4], angs[2][8], angs[2][4]))
        vol_z1 = img.dataobj[..., angs[0][8]] + img.dataobj[..., angs[1][8]] + img.dataobj[..., angs[2][8]]

        angs.sort(key=lambda tup: tup[5])
        print("pulling slices for z2: %d(%f) %d(%f) %d(%f)" % (angs[0][8], angs[0][5], angs[1][8], angs[1][5], angs[2][8], angs[2][5]))
        vol_z2 = img.dataobj[..., angs[0][8]] + img.dataobj[..., angs[1][8]] + img.dataobj[..., angs[2][8]]

        noflip_v = []
        flip_v = []

        xy_scores_f = []
        xy_scores_nf = []
        yz_scores_f = []
        yz_scores_nf = []
        xz_scores_f = []
        xz_scores_nf = []

        ###############################################################################################
        # x/y test
        print("running x/y flip test")

        p=0
        m=0
        mid=int(vol_x1.shape[2]/2)
        for i in range(vol_x1.shape[2]):
            slice1 = vol_x1[:, :, i].astype('float32')
            slice2 = vol_x2[:, :, i].astype('float32')

            slice1 = zoom(slice1, [dimX, dimY])
            slice2 = zoom(slice2, [dimX, dimY])

            #generate sample
            if i == mid:
                (image1, image2) = generate_sample(slice1, slice2)
                image1.save('secondary/xy1.png')
                image2.save('secondary/xy2.png')
          
            pos = np.subtract(slice1, slice2).clip(min=0)
            pos=np.pad(pos, ((0,0),(0, pos.shape[0])), 'constant')
            neg = np.subtract(slice2, slice1).clip(min=0)
            neg=np.pad(neg, ((0,0),(0, neg.shape[0])), 'constant')

            l=np.max(sum_diag(pos, 1))
            r=np.max(sum_diag(pos, -1))
            l+=np.max(sum_diag(neg, -1))
            r+=np.max(sum_diag(neg, 1))

            if l<=r:
                p+=1.0
                xy_scores_f.append(None)
                xy_scores_nf.append(float(r-l))
            else:
                m+=1.0
                xy_scores_f.append(float(r-l))
                xy_scores_nf.append(None)

        xy_flipped=False
        print ("  noflip", p)
        print ("  flip", m)
        noflip_v.append(p)
        flip_v.append(m)
        if p < m:
            print("x/y-flipped!")
            xy_flipped=True

        ###############################################################################################

        print("running y/z flip test")

        p=0
        m=0
        mid=int(vol_y1.shape[0]/2)
        for i in range(vol_y1.shape[0]):
            slice1 = vol_y1[i, :, :].astype('float32')
            slice2 = vol_y2[i, :, :].astype('float32')

            slice1 = zoom(slice1, [dimY, dimZ])
            slice2 = zoom(slice2, [dimY, dimZ])

            #generate sample
            if i == mid:
                (image1, image2) = generate_sample(np.fliplr(slice1), np.fliplr(slice2))
                image1.save('secondary/yz1.png')
                image2.save('secondary/yz2.png')
          
            pos = np.subtract(slice1, slice2).clip(min=0)
            pos=np.pad(pos, ((0,0),(0, pos.shape[0])), 'constant')
            neg = np.subtract(slice2, slice1).clip(min=0)
            neg=np.pad(neg, ((0,0),(0, neg.shape[0])), 'constant')

            l=np.max(sum_diag(pos, 1))
            r=np.max(sum_diag(pos, -1))
            l+=np.max(sum_diag(neg, -1))
            r+=np.max(sum_diag(neg, 1))

            if l<=r:
                p+=1.0
                yz_scores_f.append(None)
                yz_scores_nf.append(float(r-l))
            else:
                m+=1.0
                yz_scores_f.append(float(r-l))
                yz_scores_nf.append(None)

        yz_flipped=False
        print ("  noflip", p)
        print ("  flip", m)
        noflip_v.append(p)
        flip_v.append(m)
        if p < m:
            print("y/z-flipped!")
            yz_flipped=True

        ###############################################################################################
        print("running x/z flip test")

        p=0
        m=0
        mid=int(vol_z1.shape[1]/2)
        for i in range(vol_z1.shape[1]):
            slice1 = vol_z1[:, i, :].astype('float32')
            slice2 = vol_z2[:, i, :].astype('float32')

            slice1 = zoom(slice1, [dimX, dimZ])
            slice2 = zoom(slice2, [dimX, dimZ])

            #store sample image
            if i == mid:
                (image1, image2) = generate_sample(np.fliplr(slice1), np.fliplr(slice2))
                image1.save('secondary/xz1.png')
                image2.save('secondary/xz2.png')
         
            pos = np.subtract(slice1, slice2).clip(min=0)
            pos=np.pad(pos, ((0,0),(0, pos.shape[0])), 'constant')
            neg = np.subtract(slice2, slice1).clip(min=0)
            neg=np.pad(neg, ((0,0),(0, neg.shape[0])), 'constant')

            l=np.max(sum_diag(pos, 1))
            r=np.max(sum_diag(pos, -1))
            l+=np.max(sum_diag(neg, -1))
            r+=np.max(sum_diag(neg, 1))

            if l<=r:
                p+=1.0
                xz_scores_f.append(None)
                xz_scores_nf.append(float(r-l))
            else:
                m+=1.0
                xz_scores_f.append(float(r-l))
                xz_scores_nf.append(None)
            
        xz_flipped=False
        print ("  noflip", p)
        print ("  flip", m)
        noflip_v.append(p)
        flip_v.append(m)
        if p < m:
            print("x/z-flipped!")
            xz_flipped=True

        if not xy_flipped and not yz_flipped and not xz_flipped:
            print("no flip!")
        elif xy_flipped and xz_flipped:
            print("x is flipped !")
            warning("bvecs-x seems to be flipped.")
        elif xy_flipped and yz_flipped:
            print("y is flipped !")
            warning("bvecs-y seems to be flipped.")
        elif yz_flipped and xz_flipped:
            print("z is flipped !")
            warning("bvecs-z seems to be flipped.")
        else:
            print("inconclusive");
            warning("bvecs-z flipping could not be determined. Please check the data quality.")

        #output results in plotly graph
        x_labels = ['x/y', 'y/z',  'x/z']
        noflip = {
            'type': 'bar',
            'name': 'No Flip',
            'x': x_labels,
            'y': noflip_v,
        }
        flip = {
            'type': 'bar',
            'name': 'Flip',
            'x': x_labels,
            'y': flip_v,
        }
        results['brainlife'].append({
            'type': 'plotly',
            'name': 'Flip Evidence',
            'desc': 'All axies should show more no-flip evidence. The more apart the flip v.s. no-flip evidence is, the better image quality.',
            'layout': {},
            'data': [noflip, flip],
        })

        #add xy scores
        results['brainlife'].append({
            'type': 'plotly',
            'name': 'Feature stddev (x/y)',
            'desc': 'stddev computed for each slices in Z axis',
            'layout': {
                'barmode': 'stack',
                'xaxis': {
                    'title': 'z-voxel index'
                },
                'yaxis': {
                    'title': 'orientation correctness'
                }
            },
            'data': [
                {
                    'type': 'bar',
                    'name': 'no-flip',
                    #'x': x_labels, #0:i
                    'y': xy_scores_nf,
                },
                {
                    'type': 'bar',
                    'name': 'flip',
                    #'x': x_labels, #0:i
                    'y': xy_scores_f,
                }
            ],
        })

        #add yz scores
        results['brainlife'].append({
            'type': 'plotly',
            'name': 'Feature stddev (yz)',
            'desc': 'stddev computed for each slices in X axis',
            'layout': {
                'barmode': 'stack',
                'xaxis': {
                    'title': 'x-voxel index'
                },
                'yaxis': {
                    'title': 'orientation correctness'
                }
            },
            'data': [
                {
                    'type': 'bar',
                    'name': 'no-flip',
                    #'x': x_labels, #0:i
                    'y': yz_scores_nf,
                },
                {
                    'type': 'bar',
                    'name': 'flip',
                    #'x': x_labels, #0:i
                    'y': yz_scores_f,
                }
            ],
        })

        #add xz scores
        results['brainlife'].append({
            'type': 'plotly',
            'name': 'Feature stddev (xz)',
            'desc': 'stddev computed for each slices in Y axis',
            'layout': {
                'barmode': 'stack',
                'xaxis': {
                    'title': 'y-voxel index'
                },
                'yaxis': {
                    'title': 'orientation correctness'
                }
            },
            'data': [
                {
                    'type': 'bar',
                    'name': 'no-flip',
                    #'x': x_labels, #0:i
                    'y': xz_scores_nf,
                },
                {
                    'type': 'bar',
                    'name': 'flip',
                    #'x': x_labels, #0:i
                    'y': xz_scores_f,
                }
            ],
        })

    #create bvecs vis .. sort into shells (100th)
    shells = {}
    for i in range(len(gtab.bvals)):
        bval = gtab.bvals[i]
        bvec = gtab.bvecs[i]
        shell = str(round(bval, -2))
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
            'name': str(shell), 
            'x': xs,
            'y': ys,
            'z': zs,
            'text': texts,
            'textfont': {
                'color': color,
                'size': 8
            }
        })

    results['brainlife'].append({
        'type': 'plotly',
        'name': 'Gradients (bvecs/bvals)',
        'layout': {},
        'data': data,
    })

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
            np.int16, np.int32, np.int64, np.uint8,
            np.uint16, np.uint32, np.uint64)):
            ret = int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            ret = float(obj)
        elif isinstance(obj, (np.ndarray,)): 
            ret = obj.tolist()
        else:
            ret = json.JSONEncoder.default(self, obj)

        if isinstance(ret, (float)):
            if math.isnan(ret):
                ret = None

        if isinstance(ret, (bytes, bytearray)):
            ret = ret.decode("utf-8")

        return ret


with open("product.json", "w") as fp:
    json.dump(results, fp, cls=NumpyEncoder)

if len(results["errors"]) > 0:
    print("test failed")

