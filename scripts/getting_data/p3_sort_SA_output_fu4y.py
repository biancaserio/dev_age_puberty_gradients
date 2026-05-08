#!/usr/bin/env python3


# Computes total surface area by extracting data output of script p3_area2mgh.sh (surface area- available data in /pt_02667/ directory at 4Y FOLLOW-UP) and summing across hemispheres


# don't forget to activate the conda environment in order to load packages when running via terminal . /data/u_serio_software/miniforge3/etc/profile.d/conda.sh;conda activate base


########################################
### Load packages
########################################

print("----- Loading packages -----")

# General
import numpy as np
import pandas as pd
import os

# Computing
import scipy.io  # loadmat
import fnmatch  # for comparing patterns of syntax

# Neuroimaging
import nibabel as nib


########################################
### Define directories
########################################

print("----- Defining directories -----")

codedir = os.path.abspath('')  # obtain current direction from which script is runnning

datadir = '/data/pt_02667/data/ABCD/'
datadir_local = '/data/p_02667/development/data/'

resdir = '/data/p_02667/development/results/'
resdir_fig = '/data/p_02667/development/results/figures/'





########################################
### Available data in directory
########################################

print("----- Sorting through available data in fu4y directory -----")


# fetch data available in directory post datalad download
path_list_fu4y_SA = os.listdir(datadir+'ABCD_fu4y_SA')
path_list_fu4y_SA.sort()  # this contains html symlinks as well as other files that are not subject files on top

    
sub_list_fu4y_SA = []  # contains all subjects with SA output


for sub in path_list_fu4y_SA:
    
    # only considering subject files and filtering past the html broken link duplicates of each subject folder
    if sub.startswith("sub-") and "html" not in sub:
        sub_list_fu4y_SA.append(sub)




# define list that will contain the total surface area data
total_SA_fu4y = []

counter = 0

for sub in sub_list_fu4y_SA:

    counter += 1

    print(f"Subject {counter} out of {len(sub_list_fu4y_SA)}")

    sub_path_to_SA_data = datadir+'ABCD_fu4y_SA/'+sub

    # Load the .mgh file
    mgh_file_lh = nib.load(sub_path_to_SA_data+"/lh.native.area.mgh")
    mgh_file_rh = nib.load(sub_path_to_SA_data+"/rh.native.area.mgh")
    
    # Extract the image data as a NumPy array
    sa_lh = mgh_file_lh.get_fdata()
    sa_rh = mgh_file_rh.get_fdata()

    # Sum surface within and between hemispheres
    tot_sa = float(sum(sa_lh) + sum(sa_rh))

    # Append to list
    total_SA_fu4y.append(tot_sa)






########################################
### Export as csv
########################################

print(f"----- Exporting as csv at {datadir}ABCD_total_SA_fu4y.csv -----")

# Create a DataFrame
df_tot_SA_fu4y = pd.DataFrame({"src_subject_id_fmt": sub_list_fu4y_SA, "tot_SA_fu4y": total_SA_fu4y})

# Export
df_tot_SA_fu4y.to_csv(datadir+'ABCD_total_SA_fu4y.csv', header = True, index = False)


print(f"Number of subjects with Total SA output: {len(sub_list_fu4y_SA)}")
print(f"Are there any values of total SA that are NaN? -> {pd.isna(total_SA_fu4y).any()}")



