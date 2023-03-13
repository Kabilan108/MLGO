#!/bin/bash
#
# This script will download over 11000 datasets from the DEE2 repository using the
# `loaddatasets.R` script.
# Datasets will be loaded in batches of ~3500 datasets. After each batch, the
# `cleandata.R` script will be run to clean the data format the results as dataframes.
# Once all 3 batches are complete, the `prepdata.py` script will be run to prepare
# the feature matrix and target vector for the machine learning models.


for file in data/raw/batch-*-GSE.txt; do
    batch=${file#batch-}
    batch=${batch%-GSE.txt}

    echo $batch

    #Rscript loaddatasets.R "$file"

    if [[ -e data/processed/DEG_B${batch}.RData ]]; then
        echo "Batch $batch already processed"
    else
        echo ""
        #Rscript cleandata.R
    fi


    #echo "Processing $file"
    #Rscript loaddatasets.R $file
    #Rscript cleandata.R
done
