#!/bin/bash
../run.py

module load nodejs #for jsonlint
cat products.json | jsonlint
