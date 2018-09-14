# neuro/dwi datatype validator

Brainlife uses this code to validate and normalize input datasets uploaded by users to make sure that the content and format of the data matches the Brainlife's neuro/dwi datatype.

Currently this App performs following checks

* make sure nibabel runs successfully on specified dwi file
* make sure dwi is 4d
* raise warning if dwi transformation matrix isn't unit matrix (identity matrix)
* make sure bvecs and bvals can be read
* make sure bvecs's cols count matches dwi's 4th dimension number
* make sure bvecs has 3 rows
* make sure bvals's cols count matches dwi's 4th dimension number
* make sure bvals has 1 row
* check for bvecs flipping (x/z and z) - using a similar mechanism for brain-life/app-testflip

This service is not meant to be executed outside Brainlife.

In the future, this repo might be incooporated into to [brain-life/datatypes repo](https://github.com/brain-life/datatypes)

### Authors
- Soichi Hayashi (hayashis@iu.edu)

### Project directors
- Franco Pestilli (franpest@indiana.edu)

### Funding 
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
