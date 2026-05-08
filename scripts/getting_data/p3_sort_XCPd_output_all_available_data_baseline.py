#!/usr/bin/env python3


# Sorts XCP-D output for all available data in /pt_02667/ directory
# code for baseline


# don't forget to activate the conda environment in order to load packages when running via terminal . /data/u_serio_software/miniforge3/etc/profile.d/conda.sh;conda activate base


########################################
### Load packages
########################################

print("----- Loading packages -----")

# General
import numpy as np
#import pandas as pd
import os

# Computing
import scipy.io  # loadmat
import fnmatch  # for comparing patterns of syntax
import hdf5storage  # hdf5storage to read and write in HDF5 format (instead of .mat file, when matrix is too large)

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

print("----- Sorting through available data in baseline directory -----")


# fetch data available in directory post datalad download
path_list_XCPd = os.listdir(datadir+'XCP-D_output')
path_list_XCPd.sort()  # this contains html symlinks as well as other files that are not subject files on top

    
sub_list_XCPd = []  # contains all subjects with XCP-d output


for sub in path_list_XCPd:
    
    # only considering subject files and filtering past the html broken link duplicates of each subject folder
    if sub.startswith("sub-") and "html" not in sub:
        sub_list_XCPd.append(sub)



## Subjects that passed QC (minimum 2 viable runs) vs those who did not

# lists and dict to troubleshoot
sub_list_XCPd_missing_data = []  # should be empty
qc_not_passed = {'sub': [], 'number_runs': []}  # zero or 1 usable runs (I think only 1 possible but check)

# final lists with usable data
XCPd_sub_list_with_mean_corr_matrix = []
XCPd_mean_corr_matrices = []  # mean correlation matrices by subject (if minimum 2 runs) -> check that it's all the subjects, or else problem with who it is ->
XCPd_number_runs_for_mean_matrix_computation = []


# defining the pattern of file name that corresponds to XCP-d output correlation matrix
fmri_data_pattern = "*_task-rest_run-*_space-fsLR_seg-4S456Parcels_den-91k_stat-pearsoncorrelation_boldmap.pconn.nii"


counter = 0

for sub in sub_list_XCPd:

    counter += 1

    print(f"Subject {counter} out of {len(sub_list_XCPd)}")

    sub_path_to_func_data = datadir+'XCP-D_output/'+sub+'/ses-baselineYear1Arm1/func/'

    sub_corr_matrices = []
    sub_runs = []
    
    # if directory of subject contains baseline data & functional data (-> if the nested directory ".../baseline/func/" exists)
    if os.path.isdir(sub_path_to_func_data):

        for filename in os.listdir(sub_path_to_func_data):

            # if the file is a correlation matrix file
            if fnmatch.fnmatch(filename, fmri_data_pattern):

                # get the data and save it to that subject's list of matrices
                corr_matrix = np.asarray(nib.load(sub_path_to_func_data+filename).get_fdata())
                sub_corr_matrices.append(corr_matrix)

        
        # count how many runs with correlation matrix
        sub_runs = len(sub_corr_matrices)
        
        if sub_runs >= 2:
    
            # make array of correlation matrices
            sub_corr_matrices = np.array(sub_corr_matrices)

            # Compute mean correlation matrix for subject (mean matrix - i.e., averaged across runs)
            sub_mean_corr_matrix = np.mean(sub_corr_matrices, axis=0)
        
            # save subject's mean matrix to the overall list of matrices
            XCPd_mean_corr_matrices.append(np.array(sub_mean_corr_matrix))

            # save subject's name to list (corresponding to mean matrices
            XCPd_sub_list_with_mean_corr_matrix.append(sub)

            XCPd_number_runs_for_mean_matrix_computation.append(sub_runs)

        else:
            qc_not_passed['sub'].append(sub)
            qc_not_passed['number_runs'].append(sub_runs)
                
    
    else:
        sub_list_XCPd_missing_data.append(sub)  # I don't think any XCP directories would have missing baseline/func directory but just in case



# Make final arrays containing mean correlation matrices
XCPd_mean_corr_matrices_with_subcortex = np.array(XCPd_mean_corr_matrices)  # shape 456x456 (with subcortex)
XCPd_mean_corr_matrices = np.array(XCPd_mean_corr_matrices)[:, :400, :400]  # shape 400x400 (only cortical)


########################################
### Export as matfile 
########################################

#print("----- Exporting as matfile -----")

#mdict = {'fc_matrices': XCPd_mean_corr_matrices, 'sub_ID': XCPd_sub_list_with_mean_corr_matrix, 'runs_count': #XCPd_number_runs_for_mean_matrix_computation}

#scipy.io.savemat(datadir+'ABCD_baseline_fc_matrices_XCP-d_output.mat', mdict)


print("----- Exporting as HDF5 file (matrix too large for Matlab 5 format (.mat)  -----")

mdict = {
    'fc_matrices': XCPd_mean_corr_matrices,
    'sub_ID': XCPd_sub_list_with_mean_corr_matrix,
    'runs_count': XCPd_number_runs_for_mean_matrix_computation
}

hdf5storage.write(mdict, filename=datadir+'ABCD_baseline_fc_matrices_XCP-d_output.h5', format='7.3')


print(f"Number of subjects with XCP-d output: {len(sub_list_XCPd)}")
print(f"Number of subjects with minimum 2 runs of QC'd data: {len(XCPd_mean_corr_matrices)}")
print(f"Number of subjects who did NOT pass QC (only 1 viable run) -> discard subject (no mean matrix): {len(qc_not_passed)}")

