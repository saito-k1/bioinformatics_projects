#!/usr/bin/env bash
# cancerGenomics.sh

R -e "rmarkdown::render('cancerGenomics.Rmd', output_format='all')"
