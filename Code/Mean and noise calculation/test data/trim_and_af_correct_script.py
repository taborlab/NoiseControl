# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 10:13:49 2017

@author: Evan Olson
"""

# Add path to the trim_and_af_correct.py file
import sys
sys.path.append('../trimming_and_af_correction/')

# Import the trim_and_af_correct file
import trim_and_af_correct as trimaf

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


### ANALYSIS SPECIFICATION ###

### 2020_03_13 note: Now when you change the below variables 'channel_w' and 'channel' to FL3, the analysis should properly process the FL3 data! -EO

# White cell specification
fc_data_path_w = '../test data/experiment_output.xlsx'
sample_ids_w = ['S0002'] # White cell sample in experiment output file
channel_w = 'FL1' #Set to 'FL1' or 'FL3'

# Sample specification
fc_data_path = 'experiment_output.xlsx'
channel = 'FL1' #Set to 'FL1' or 'FL3'

## Sample IDs can be directly specified:
# sample_ids = ['S0001', 'S0002', 'S0011'] # uncomment for specified file names only

## Or sample IDs can be determined by the fc output file:
df = pd.read_excel(fc_data_path, sheet_name='Samples') # comment for specified file names only
sample_ids = df['ID'].values # comment for specified file names only

# Output specification
trim_plots_folder = 'plt/'

# Trimming algorithm parameters
## White cells
thresh_w = 0.005
bw_w = 0.05
## Samples
thresh = 0.005
bw = 0.05

### ANALYSIS SCRIPT ###
# Import white cell data
data_w = trimaf.import_hists(fc_data_path_w, sample_ids_w, channel_w)

# Import sample data
data = trimaf.import_hists(fc_data_path, sample_ids, channel)

# Trim white cell data
## Default values for thresh and bandwidth are 0.005 and 0.05, but can be overwritten
trimmed_w = trimaf.trim(data_w, thresh=thresh_w, bw=bw_w)

# Save white-cell trimming plots
## Prefix is used for the white cells to distinguish them from the other samples
trimaf.plot(trimmed_w, trim_plots_folder, prefix='white')

# Trim sample data
trimmed = trimaf.trim(data, thresh=thresh, bw=bw)

# Save sample trimming plots
trimaf.plot(trimmed, trim_plots_folder)

# Calculate af-corrected sample summary statistics
af_corrected = trimaf.af_correct(trimmed, trimmed_w)

# Save results
trimaf.save(af_corrected, fc_data_path)
