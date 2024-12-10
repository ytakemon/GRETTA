#!/bin/bash -v

# Step 1: Set up Conda environment
conda env create -f environment.yaml

# Step 2: Activate the environment
conda activate gretta_env

# Step 3: Run R setup script - installs gretta from github
Rscript setup.R

# Step 4: Download depmap files
folder_name="GRETTA_project"
if [ ! -f data/dataset.csv ]; then
    mkdir -p ${folder_name}
    wget -O ${folder_name}/GRETTA_DepMap_23Q4_data.tar.gz  https://www.bcgsc.ca/downloads/ytakemon/GRETTA/23Q4/GRETTA_DepMap_23Q4_data.tar.gz
    cd ${folder_name}
    tar -xzvf GRETTA_DepMap_23Q4_data.tar.gz
fi

