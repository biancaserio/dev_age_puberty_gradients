#!/usr/bin/env python3


# Script that gets final sample (as reported in data trees)
# - Gets data from XCPd sorted output
# - Finds and removes NaNs in FC matrices
# - Finds and removes subjects with missing parental PDS score from clean final sample
# - Adds CBCL scores (total, internalizing, externalizing) to clean final sample (even if they may have NaNs - they will be excluded for the CBCL analyses later)
# - Adds FES scores (fam_conflict_y) to clean final sample (even if they may have NaNs - they will be excluded for the FES analyses later)
# - Computes Fisher r-to-z transfomration so that there is also that option in FC matrices files
# Does this for all 3 timepoints

# Also created "full" fc matrices array and demographics dataframe (concatenated across timepoints)

# don't forget to activate the conda environment in order to load packages when running via terminal . /data/u_serio_software/miniforge3/etc/profile.d/conda.sh;conda activate base


########################################
### Load packages
########################################

print("\n\n\n\n")

print("----- Loading packages -----")

# General
import numpy as np
import pandas as pd
import os

# Computing
import scipy.io  # loadmat
import hdf5storage  # hdf5storage to write in HDF5 format (instead of .mat file, when matrix is too large) -> hdf5storage.write() hdf5storage.read()

print("\n\n\n\n")



########################################
### Define directories
########################################

print("----- Defining directories -----")

codedir = os.path.abspath('')  # obtain current direction from which script is runnning

datadir = '/data/pt_02667/data/ABCD/'
datadir_local = '/data/p_02667/development/data/'

resdir = '/data/p_02667/development/results/'
resdir_fig = '/data/p_02667/development/results/figures/'

print("\n\n\n\n")




########################################
### FC Matrices
########################################

print("----- FC Matrices  -----")

print("\n\n### Baseline ###\n\n")

# Loading HDF file (file was too large for matfile), with HDF reader e.g., hdf5storage
hdf_file_fc_matrices_baseline = hdf5storage.read(filename=datadir+'ABCD_baseline_fc_matrices_XCP-d_output.h5')


fc_matrices_baseline = hdf_file_fc_matrices_baseline['fc_matrices']
sub_ID_baseline = hdf_file_fc_matrices_baseline['sub_ID']
#runs_count_baseline = hdf_file_fc_matrices_baseline['runs_count']


print(f"Shape of FC matrices from XCP-d output at baseline: {fc_matrices_baseline.shape}")


# Identify subjects with any NaNs in their matrices
nan_subjects = [i for i in range(fc_matrices_baseline.shape[0]) if np.isnan(fc_matrices_baseline[i]).any()]

print(f"Number of subjects with NaNs in their FC matrices: {len(nan_subjects)}")
print(f"Subjects with NaNs in their FC matrices: {nan_subjects}")


# Deleting from the array the subjects with NaNs
fc_matrices_baseline = np.delete(fc_matrices_baseline, nan_subjects, axis=0)
sub_ID_baseline = np.delete(sub_ID_baseline, nan_subjects, axis=0)


print(f"Shape FC matrices at baseline: {fc_matrices_baseline.shape}")



print("\n\n### 2y follow-up ###\n\n")

# loading matfile exports
fc_matrices_fu2y = scipy.io.loadmat(datadir+'ABCD_fu2y_fc_matrices_XCP-d_output.mat')['fc_matrices'] 
sub_ID_fu2y = scipy.io.loadmat(datadir+'ABCD_fu2y_fc_matrices_XCP-d_output.mat')['sub_ID']
#runs_count_fu2y = scipy.io.loadmat(datadir+'ABCD_fu2y_fc_matrices_XCP-d_output.mat')['runs_count']


print(f"Shape of FC matrices from XCP-d output at fu2y: {fc_matrices_fu2y.shape}")


# Identify subjects with any NaNs in their matrices
nan_subjects = [i for i in range(fc_matrices_fu2y.shape[0]) if np.isnan(fc_matrices_fu2y[i]).any()]

print(f"Number of subjects with NaNs in their FC matrices: {len(nan_subjects)}")
print(f"Subjects with NaNs in their FC matrices: {nan_subjects}")


# Deleting from the array the subjects with NaNs
fc_matrices_fu2y = np.delete(fc_matrices_fu2y, nan_subjects, axis=0)
sub_ID_fu2y = np.delete(sub_ID_fu2y, nan_subjects, axis=0)


print(f"Shape FC matrices at fu2y: {fc_matrices_fu2y.shape}")


# intersecting subjects with FC data for baseline and 2y follow-up
common = set(sub_ID_baseline) & set(sub_ID_fu2y)
print(f"Number of subjects with both baseline and 2y follow-up FC data: {len(common)}")



print("\n\n### 4y follow-up ###\n\n")

# loading matfile exports
fc_matrices_fu4y = scipy.io.loadmat(datadir+'ABCD_fu4y_fc_matrices_XCP-d_output.mat')['fc_matrices'] 
sub_ID_fu4y = scipy.io.loadmat(datadir+'ABCD_fu4y_fc_matrices_XCP-d_output.mat')['sub_ID']
#runs_count_fu4y = scipy.io.loadmat(datadir+'ABCD_fu4y_fc_matrices_XCP-d_output.mat')['runs_count']


print(f"Shape of FC matrices from XCP-d output at fu4y: {fc_matrices_fu4y.shape}")


# Identify subjects with any NaNs in their matrices
nan_subjects = [i for i in range(fc_matrices_fu4y.shape[0]) if np.isnan(fc_matrices_fu4y[i]).any()]

print(f"Number of subjects with NaNs in their FC matrices: {len(nan_subjects)}")
print(f"Subjects with NaNs in their FC matrices: {nan_subjects}")


# Deleting from the array the subjects with NaNs
fc_matrices_fu4y = np.delete(fc_matrices_fu4y, nan_subjects, axis=0)
sub_ID_fu4y = np.delete(sub_ID_fu4y, nan_subjects, axis=0)


print(f"Shape FC matrices at fu4y: {fc_matrices_fu4y.shape}")


# intersecting subjects with FC data for baseline, 2y follow-up and 4y follow-up
common = set(sub_ID_baseline) & set(sub_ID_fu2y) & set(sub_ID_fu4y)
print(f"Number of subjects with baseline, 2y and 4y follow-up FC data: {len(common)}")


print("\n\n\n\n")




########################################
### Demographics and Puberty
########################################

print("----- Demographics and Puberty  -----")

print("\n\n### Baseline ###\n\n")


# Load puberty and demographics data
abcd_puberty_baseline = pd.read_csv(datadir_local+'abcd_puberty_baseline.csv')


## take a temporary subset of the dataframe including subjects with available FC data & in order correspnding to FC data

# subset the DataFrame to include only rows where subject id from dataframe is in subject id of FC matrices
abcd_puberty_baseline_temp = abcd_puberty_baseline[abcd_puberty_baseline['src_subject_id_fmt'].isin(sub_ID_baseline)].copy()

# Reorder the rows of `df_subset` to match the order in sub_ID_baseline
abcd_puberty_baseline_temp = abcd_puberty_baseline_temp.set_index('src_subject_id_fmt')  # Set 'src_subject_id_fmt' as the index
abcd_puberty_baseline_temp = abcd_puberty_baseline_temp.reindex(sub_ID_baseline)  # Reorder the rows based on sub_ID_baseline

# Reset the index if you want 'subject_id' back as a column
abcd_puberty_baseline_temp = abcd_puberty_baseline_temp.reset_index()

print(f"{np.sum(abcd_puberty_baseline_temp.sex == 'F')} females, {np.sum(abcd_puberty_baseline_temp.sex == 'M')} males")


## PDS strictly parental report -> find NaN values

nan_indices_p_pds = np.where(np.isnan(abcd_puberty_baseline_temp.pds_p_score))[0]
nan_indices_p_pds_M = np.where(np.isnan(abcd_puberty_baseline_temp[abcd_puberty_baseline_temp['sex'] == 'M'].pds_p_score))[0]
nan_indices_p_pds_F = np.where(np.isnan(abcd_puberty_baseline_temp[abcd_puberty_baseline_temp['sex'] == 'F'].pds_p_score))[0]

print(f"There are {len(nan_indices_p_pds)} NaN values for the parental PDS report in total - indices: {nan_indices_p_pds}, of which:")
print(f"    - {len(nan_indices_p_pds_M)} are male")
print(f"    - {len(nan_indices_p_pds_F)} are female")


print("\n...Excluding these subjects with no parental PDS score...\n")


### Excluding these subjects with no parental PDS score

# Dataframe
abcd_puberty_baseline_temp = abcd_puberty_baseline_temp.drop(nan_indices_p_pds, axis=0).reset_index(drop=True)  # really important to reset index everytime you drop, otherwise dropping via index != dropping via rows for later steps

# FC matrices
fc_matrices_baseline = np.delete(fc_matrices_baseline, nan_indices_p_pds, axis=0)

# Subject IDs
sub_ID_baseline = np.delete(sub_ID_baseline, nan_indices_p_pds, axis=0)


print(f"Number of subjects at baseline: {len(fc_matrices_baseline)}")



print("\n\n### 2y follow-up ###\n\n")


# Load puberty and demographics data
abcd_puberty_fu2y = pd.read_csv(datadir_local+'abcd_puberty_fu2y.csv')


## take a temporary subset of the dataframe including subjects with available FC data & in order correspnding to FC data

# subset the DataFrame to include only rows where subject id from dataframe is in subject id of FC matrices
abcd_puberty_fu2y_temp = abcd_puberty_fu2y[abcd_puberty_fu2y['src_subject_id_fmt'].isin(sub_ID_fu2y)].copy()

# Reorder the rows of `df_subset` to match the order in sub_ID_fu2y
abcd_puberty_fu2y_temp = abcd_puberty_fu2y_temp.set_index('src_subject_id_fmt')  # Set 'src_subject_id_fmt' as the index
abcd_puberty_fu2y_temp = abcd_puberty_fu2y_temp.reindex(sub_ID_fu2y)  # Reorder the rows based on sub_ID_fu2y

# Reset the index if you want 'subject_id' back as a column
abcd_puberty_fu2y_temp = abcd_puberty_fu2y_temp.reset_index()

print(f"{np.sum(abcd_puberty_fu2y_temp.sex == 'F')} females, {np.sum(abcd_puberty_fu2y_temp.sex == 'M')} males")


## PDS strictly parental report -> find NaN values

nan_indices_p_pds = np.where(np.isnan(abcd_puberty_fu2y_temp.pds_p_score))[0]
nan_indices_p_pds_M = np.where(np.isnan(abcd_puberty_fu2y_temp[abcd_puberty_fu2y_temp['sex'] == 'M'].pds_p_score))[0]
nan_indices_p_pds_F = np.where(np.isnan(abcd_puberty_fu2y_temp[abcd_puberty_fu2y_temp['sex'] == 'F'].pds_p_score))[0]

print(f"There are {len(nan_indices_p_pds)} NaN values for the parental PDS report in total - indices: {nan_indices_p_pds}, of which:")
print(f"    - {len(nan_indices_p_pds_M)} are male")
print(f"    - {len(nan_indices_p_pds_F)} are female")


print("\n...Excluding these subjects with no parental PDS score...\n")


### Excluding these subjects with no parental PDS score

# Dataframe
abcd_puberty_fu2y_temp = abcd_puberty_fu2y_temp.drop(nan_indices_p_pds, axis=0).reset_index(drop=True)  # really important to reset index everytime you drop, otherwise dropping via index != dropping via rows for later steps

# FC matrices
fc_matrices_fu2y = np.delete(fc_matrices_fu2y, nan_indices_p_pds, axis=0)

# Subject IDs
sub_ID_fu2y = np.delete(sub_ID_fu2y, nan_indices_p_pds, axis=0)


print(f"Number of subjects at fu2y: {len(fc_matrices_fu2y)}")


print("\n\n### 4y follow-up ###\n\n")


# Load puberty and demographics data
abcd_puberty_fu4y = pd.read_csv(datadir_local+'abcd_puberty_fu4y.csv')


## take a temporary subset of the dataframe including subjects with available FC data & in order correspnding to FC data

# subset the DataFrame to include only rows where subject id from dataframe is in subject id of FC matrices
abcd_puberty_fu4y_temp = abcd_puberty_fu4y[abcd_puberty_fu4y['src_subject_id_fmt'].isin(sub_ID_fu4y)].copy()

# Reorder the rows of `df_subset` to match the order in sub_ID_fu4y
abcd_puberty_fu4y_temp = abcd_puberty_fu4y_temp.set_index('src_subject_id_fmt')  # Set 'src_subject_id_fmt' as the index
abcd_puberty_fu4y_temp = abcd_puberty_fu4y_temp.reindex(sub_ID_fu4y)  # Reorder the rows based on sub_ID_fu4y

# Reset the index if you want 'subject_id' back as a column
abcd_puberty_fu4y_temp = abcd_puberty_fu4y_temp.reset_index()

print(f"{np.sum(abcd_puberty_fu4y_temp.sex == 'F')} females, {np.sum(abcd_puberty_fu4y_temp.sex == 'M')} males")



## PDS strictly parental report -> find NaN values

nan_indices_p_pds = np.where(np.isnan(abcd_puberty_fu4y_temp.pds_p_score))[0]
nan_indices_p_pds_M = np.where(np.isnan(abcd_puberty_fu4y_temp[abcd_puberty_fu4y_temp['sex'] == 'M'].pds_p_score))[0]
nan_indices_p_pds_F = np.where(np.isnan(abcd_puberty_fu4y_temp[abcd_puberty_fu4y_temp['sex'] == 'F'].pds_p_score))[0]

print(f"There are {len(nan_indices_p_pds)} NaN values for the parental PDS report in total - indices: {nan_indices_p_pds}, of which:")
print(f"    - {len(nan_indices_p_pds_M)} are male")
print(f"    - {len(nan_indices_p_pds_F)} are female")


print("\n...Excluding these subjects with no parental PDS score...\n")


### Excluding these subjects with no parental PDS score

# Dataframe
abcd_puberty_fu4y_temp = abcd_puberty_fu4y_temp.drop(nan_indices_p_pds, axis=0).reset_index(drop=True)  # really important to reset index everytime you drop, otherwise dropping via index != dropping via rows for later steps

# FC matrices
fc_matrices_fu4y = np.delete(fc_matrices_fu4y, nan_indices_p_pds, axis=0)

# Subject IDs
sub_ID_fu4y = np.delete(sub_ID_fu4y, nan_indices_p_pds, axis=0)


print(f"Number of subjects at fu4y: {len(fc_matrices_fu4y)}")


print("\n\n\n\n")


########################################
### Total Surface Area
########################################

print("----- Total Surface Area  -----")
# Checking for missing values (there should be none) and adding tot SA to demographics dataframe...


print("\n\n### Baseline ###\n\n")

# Load total SA data
total_SA_baseline = pd.read_csv(datadir+'ABCD_total_SA_baseline.csv')


# checking missing values 
missing_from_FS = set(sub_ID_baseline) - set(total_SA_baseline.src_subject_id_fmt)

print(f"Number of missing values (there should be none): {len(missing_from_FS)}")


if len(missing_from_FS)>0:

    ### Excluding these subjects with no total surface area
    print("Excluding subjects with no total surface area")

    # Find out form the df the indices of the IDs to drop
    indices_to_drop = abcd_puberty_baseline_temp[abcd_puberty_baseline_temp['src_subject_id_fmt'].isin(missing_from_FS)].index

    # Dataframe
    # Drop rows where src_subject_id_fmt is in the list
    abcd_puberty_baseline_temp = abcd_puberty_baseline_temp.drop(indices_to_drop.to_numpy(), axis=0).reset_index(drop=True)  # really important to reset index everytime you drop, otherwise dropping via index != dropping via rows for later steps
    
    # FC matrices
    fc_matrices_baseline = np.delete(fc_matrices_baseline, indices_to_drop.to_numpy(), axis=0)
    
    # Subject IDs
    sub_ID_baseline = np.delete(sub_ID_baseline, indices_to_drop.to_numpy(), axis=0)


# adding the total SA values to the puberty dataframe containing all data
abcd_demo_baseline_temp = abcd_puberty_baseline_temp.merge(total_SA_baseline, on="src_subject_id_fmt", how="left")
abcd_demo_baseline_temp = abcd_demo_baseline_temp.rename(columns={"tot_SA_baseline": "tot_SA"})


print(f"\n\nCheck that subject IDs in FC data and demo tables match (baseline): {np.array_equal(np.array(sub_ID_baseline), np.array(abcd_demo_baseline_temp.src_subject_id_fmt))}")



print("\n\n### Fu2y ###\n\n")

# Load total SA data
total_SA_fu2y = pd.read_csv(datadir+'ABCD_total_SA_fu2y.csv')


# checking missing values 
missing_from_FS = set(sub_ID_fu2y) - set(total_SA_fu2y.src_subject_id_fmt)
print(f"Number of missing values (there should be none): {len(missing_from_FS)}")


if len(missing_from_FS)>0:

    ### Excluding these subjects with no total surface area
    print("Excluding subjects with no total surface area")

    # Find out form the df the indices of the IDs to drop
    indices_to_drop = abcd_puberty_fu2y_temp[abcd_puberty_fu2y_temp['src_subject_id_fmt'].isin(missing_from_FS)].index
    
    # Dataframe
    # Drop rows where src_subject_id_fmt is in the list
    abcd_puberty_fu2y_temp = abcd_puberty_fu2y_temp.drop(indices_to_drop, axis=0).reset_index(drop=True)  # really important to reset index everytime you drop, otherwise dropping via index != dropping via rows for later steps
    
    # FC matrices
    fc_matrices_fu2y = np.delete(fc_matrices_fu2y, indices_to_drop.to_numpy(), axis=0)
    
    # Subject IDs
    sub_ID_fu2y = np.delete(sub_ID_fu2y, indices_to_drop.to_numpy(), axis=0)


# adding the total SA values to the puberty dataframe containing all data
abcd_demo_fu2y_temp = abcd_puberty_fu2y_temp.merge(total_SA_fu2y, on="src_subject_id_fmt", how="left")
abcd_demo_fu2y_temp = abcd_demo_fu2y_temp.rename(columns={"tot_SA_fu2y": "tot_SA"})


print(f"\n\nCheck that subject IDs in FC data and demo tables match (fu2y): {np.array_equal(np.array(sub_ID_fu2y), np.array(abcd_demo_fu2y_temp.src_subject_id_fmt))}")


# intersecting subjects with full data available for baseline and 2y follow-up 
common = set(sub_ID_baseline) & set(sub_ID_fu2y)
print(f"\nNumber of subjects with both baseline and 2y follow-up FC and PDS data: {len(common)}")


print("\n\n### Fu4y ###\n\n")

# Load total SA data
total_SA_fu4y = pd.read_csv(datadir+'ABCD_total_SA_fu4y.csv')


# checking missing values 
missing_from_FS = set(sub_ID_fu4y) - set(total_SA_fu4y.src_subject_id_fmt)
print(f"Number of missing values (there should be none): {len(missing_from_FS)}")


if len(missing_from_FS)>0:

    ### Excluding these subjects with no total surface area
    print("Excluding subjects with no total surface area")

    # Find out form the df the indices of the IDs to drop
    indices_to_drop = abcd_puberty_fu4y_temp[abcd_puberty_fu4y_temp['src_subject_id_fmt'].isin(missing_from_FS)].index
    
    # Dataframe
    # Drop rows where src_subject_id_fmt is in the list
    abcd_puberty_fu4y_temp = abcd_puberty_fu4y_temp.drop(indices_to_drop, axis=0).reset_index(drop=True)  # really important to reset index everytime you drop, otherwise dropping via index != dropping via rows for later steps
    
    # FC matrices
    fc_matrices_fu4y = np.delete(fc_matrices_fu4y, indices_to_drop.to_numpy(), axis=0)
    
    # Subject IDs
    sub_ID_fu4y = np.delete(sub_ID_fu4y, indices_to_drop.to_numpy(), axis=0)


# adding the total SA values to the puberty dataframe containing all data
abcd_demo_fu4y_temp = abcd_puberty_fu4y_temp.merge(total_SA_fu4y, on="src_subject_id_fmt", how="left")
abcd_demo_fu4y_temp = abcd_demo_fu4y_temp.rename(columns={"tot_SA_fu4y": "tot_SA"})


print(f"\n\nCheck that subject IDs in FC data and demo tables match (fu4y): {np.array_equal(np.array(sub_ID_fu4y), np.array(abcd_demo_fu4y_temp.src_subject_id_fmt))}")



# intersecting subjects with FC & PDS data for baseline, 2y follow-up and 4y follow-up
common = set(sub_ID_baseline) & set(sub_ID_fu2y) & set(sub_ID_fu4y)
print(f"\nNumber of subjects with both baseline, 2y and 4y follow-up FC and PDS data: {len(common)}")


print("\n\n\n\n")



########################################
### CBCL
########################################


print("----- CBCL  -----")
# The only thing I am doing is merging the available data for the above final sample (ie not checking for NaNs, because analyses involving CBCL are secondary and will be in a separate subsample
# So in the current merge (left join), I am preserving the subjects of abcd_demo_XXX_temp dataframes (not the subjects of CBCL), and for subjects for whom there is no rows in CBCL, this gets set to NaN


print("\n\n### Baseline ###\n\n")

# Load CBCL data
abcd_cbcl_baseline = pd.read_csv(datadir_local+'abcd_cbcl_baseline.csv')

# Merge
abcd_demo_baseline_temp2 = abcd_demo_baseline_temp.merge(abcd_cbcl_baseline, on=["src_subject_id_fmt", "eventname"], how="left")



print("\n\n### Fu2y ###\n\n")

# Load CBCL data
abcd_cbcl_fu2y = pd.read_csv(datadir_local+'abcd_cbcl_fu2y.csv')

# Merge
abcd_demo_fu2y_temp2 = abcd_demo_fu2y_temp.merge(abcd_cbcl_fu2y, on=["src_subject_id_fmt", "eventname"], how="left")



print("\n\n### Fu4y ###\n\n")

# Load CBCL data
abcd_cbcl_fu4y = pd.read_csv(datadir_local+'abcd_cbcl_fu4y.csv')

# Merge
abcd_demo_fu4y_temp2 = abcd_demo_fu4y_temp.merge(abcd_cbcl_fu4y, on=["src_subject_id_fmt", "eventname"], how="left")



print("\n\n\n\n")




########################################
### FES
########################################


print("----- FES  -----")
# The only thing I am doing is merging the available data for the above final sample (ie not checking for NaNs, because analyses involving FES are secondary and will be in a separate subsample
# So in the current merge (left join), I am preserving the subjects of abcd_demo_XXX_temp dataframes (not the subjects of FES), and for subjects for whom there is no rows in FES, this gets set to NaN


print("\n\n### Baseline ###\n\n")

# Load FES data
abcd_fes_baseline = pd.read_csv(datadir_local+'abcd_fes_baseline.csv')

# Merge
abcd_demo_baseline = abcd_demo_baseline_temp2.merge(abcd_fes_baseline, on=["src_subject_id_fmt", "eventname"], how="left")



print("\n\n### Fu2y ###\n\n")

# Load FES data
abcd_fes_fu2y = pd.read_csv(datadir_local+'abcd_fes_fu2y.csv')

# Merge
abcd_demo_fu2y = abcd_demo_fu2y_temp2.merge(abcd_fes_fu2y, on=["src_subject_id_fmt", "eventname"], how="left")



print("\n\n### Fu4y ###\n\n")

# Load FES data
abcd_fes_fu4y = pd.read_csv(datadir_local+'abcd_fes_fu4y.csv')

# Merge
abcd_demo_fu4y = abcd_demo_fu4y_temp2.merge(abcd_fes_fu4y, on=["src_subject_id_fmt", "eventname"], how="left")



print("\n\n\n\n")


########################################
### Creating r-to-z transformed fc matrices in order to also have that option in arrays to export
########################################

print("----- Fisher r-to-z transformation of FC matrices  -----")

print("\n\n### Baseline ###\n\n")


fc_matrices_baseline_z = np.arctanh(fc_matrices_baseline)
fc_matrices_baseline_z[np.isinf(fc_matrices_baseline_z)] = 0 


print("\n\n### fu2y ###\n\n")

fc_matrices_fu2y_z = np.arctanh(fc_matrices_fu2y)
fc_matrices_fu2y_z[np.isinf(fc_matrices_fu2y_z)] = 0 


print("\n\n### fu4y ###\n\n")

fc_matrices_fu4y_z = np.arctanh(fc_matrices_fu4y)
fc_matrices_fu4y_z[np.isinf(fc_matrices_fu4y_z)] = 0 





########################################
### Merging fc matrices arrays and demographics dataframes to create full sample ones
########################################



print("----- Merging fc matrices arrays and demographics dataframes to create full sample ones  -----")

# merge demographics dataframes across timepoints
abcd_demo_full = pd.concat(
    [abcd_demo_baseline, abcd_demo_fu2y, abcd_demo_fu4y],
    axis=0,   # stack rows
    ignore_index=True  # reset the index
)



# merge the fc matrices / corresponding subject IDs across timepoints

fc_matrices_full = np.concatenate(
    [fc_matrices_baseline, fc_matrices_fu2y, fc_matrices_fu4y],
    axis=0
)


sub_ID_full = np.concatenate(
    [sub_ID_baseline, sub_ID_fu2y, sub_ID_fu4y],
    axis=0
)





########################################
### Final Samples Cleaned
########################################

print("\n")

print(f"Demographics of sample at baseline: N = {len(sub_ID_baseline)} ({np.sum(abcd_demo_baseline.sex == 'F')} females, {np.sum(abcd_demo_baseline.sex == 'M')} males)")

print(f"Demographics of sample at fu2y: N = {len(sub_ID_fu2y)} ({np.sum(abcd_demo_fu2y.sex == 'F')} females, {np.sum(abcd_demo_fu2y.sex == 'M')} males)")

print(f"Demographics of sample at fu4y: N = {len(sub_ID_fu4y)} ({np.sum(abcd_demo_fu4y.sex == 'F')} females, {np.sum(abcd_demo_fu4y.sex == 'M')} males)")

print(f"Demographics of full sample: N = {len(sub_ID_full)} ({np.sum(abcd_demo_full.sex == 'F')} females, {np.sum(abcd_demo_full.sex == 'M')} males)")


print("\n\n\n\n")








########################################
### Exporting data
########################################

print("----- Exporting Cleaned Data  -----")

print("Demographics as csv\n")

print(f"-> Baseline demographics at: {datadir_local}abcd_demo_baseline_clean.csv")
abcd_demo_baseline.to_csv(datadir_local+'abcd_demo_baseline_clean.csv', header = True, index = False)


print(f"-> Fu2y demographics at: {datadir_local}abcd_demo_fu2y_clean.csv")
abcd_demo_fu2y.to_csv(datadir_local+'abcd_demo_fu2y_clean.csv', header = True, index = False)


print(f"-> Fu4y demographics at: {datadir_local}abcd_demo_fu4y_clean.csv")
abcd_demo_fu4y.to_csv(datadir_local+'abcd_demo_fu4y_clean.csv', header = True, index = False)


print(f"-> Full demographics at: {datadir_local}abcd_demo_full_clean.csv")
abcd_demo_full.to_csv(datadir_local+'abcd_demo_full_clean.csv', header = True, index = False)







print("\n\n")



print("FC as HDF5 file\n")

print(f"-> Baseline FC at: {datadir_local}abcd_fc_matrices_baseline_clean.h5")


mdict = {
    'fc_matrices_baseline': fc_matrices_baseline,
    'fc_matrices_baseline_z': fc_matrices_baseline_z,
    'sub_ID_baseline': sub_ID_baseline
}

# !!! will probably throw error to overwrite - should be able to bypass this by adding the parameter: truncate_existing=True   # <-- ensures existing file is overwritten
hdf5storage.write(mdict, filename=datadir_local+'abcd_fc_matrices_baseline_clean.h5', format='7.3')



print(f"-> Fu2y FC at: {datadir_local}abcd_fc_matrices_fu2y_clean.h5")

mdict = {
    'fc_matrices_fu2y': fc_matrices_fu2y,
    'fc_matrices_fu2y_z': fc_matrices_fu2y_z,
    'sub_ID_fu2y': sub_ID_fu2y
}


hdf5storage.write(mdict, filename=datadir_local+'abcd_fc_matrices_fu2y_clean.h5', format='7.3')




print(f"-> Fu4y FC at: {datadir_local}abcd_fc_matrices_fu4y_clean.h5")


mdict = {
    'fc_matrices_fu4y': fc_matrices_fu4y,
    'fc_matrices_fu4y_z': fc_matrices_fu4y_z,
    'sub_ID_fu4y': sub_ID_fu4y
}


hdf5storage.write(mdict, filename=datadir_local+'abcd_fc_matrices_fu4y_clean.h5', format='7.3')




print(f"-> Full FC at: {datadir_local}abcd_fc_matrices_full_clean.h5")

mdict = {
    'fc_matrices_full': np.array(fc_matrices_full),
    'sub_ID_full': np.array(sub_ID_full.tolist(), dtype=object)
}


hdf5storage.write(mdict, filename=datadir_local+'abcd_fc_matrices_full_clean.h5', format='7.3')



