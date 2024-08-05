#!/bin/bash
cd /Users/mdarifulislam/Library/CloudStorage/OneDrive-IndianaUniversity/Research/Github/FM2STRESS/FM2STRESS_project/code/InteractiveFM/SKHASH2/SKHASH

# Activate conda pythoon environment named 'obspy'
# source ~/miniconda3/etc/profile.d/conda.sh
conda init
conda activate obspy

python3 SKHASH.py None
