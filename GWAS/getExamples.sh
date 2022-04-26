#!/usr/bin/env bash
# getExamples.sh

# Retrieves the example data in GWAS directory and unzips the file

wget \
  -O "hapmap1.zip" \
  "https://zzz.bwh.harvard.edu/plink/hapmap1.zip"

unzip "hapmap1.zip"
