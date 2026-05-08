# Age and pubertal stage effects on functional brain organization

### This is the repository for the preprint:
Bianca Serio, Lars Dinkelbach, Mylla Marsiglia, Laura Waite, Felix Hoffstaedter, Daniel S. Margulies, Simon B. Eickhoff, & Sofie L. Valk (2026) **Individual differences reveal distinct age and pubertal contributions to the refinement of the functional cortical hierarchy during adolescence**. https://doi.org/10.64898/2026.05.07.723547. 

## Scripts

**1. MRI processing commands**
- `mri_processing/fmriprep.sh`
- `mri_processing/freesurfer_recon-all.sh`
- `mri_processing/xcp-d_postprocessing.sh` 
- `mri_processing/custom_confounds_24P_csf_wm.yaml` 

**2. Getting & preparing the data**
- `getting_data/get_ABCD_data_fmriprep_xcp-d_fc_specifiedsubs_*.sh` gets and postprocesses ABCD fMRIprep files from DataLad using XCP-D - outputs FC matrices  (* scripts per ABCD timepoint)
- `getting_data/p3_sort_XCPd_output_all_available_data_*.py` sorts XCP-D output for all available data in /pt_02667/ directory  (* scripts per ABCD timepoint)
- `getting_data/get_ABCD_data_freesurfer_specified_files_all_available_subs_*.sh` downloades ABCD data freesurfer files from DataLad (* scripts per ABCD timepoint)
- `getting_data/p3_area2mgh.sh` converts native FreeSurfer surfaces .area to .mgh format (for total surface area computation)
- `getting_data/p3_sort_SA_output_*.py` computes total surface area by extracting data output of script p3_area2mgh.sh and summing across hemispheres (* scripts per ABCD timepoint)

- `p3_for_gd_create_midthickness_fsaverage5.sh` create a template midthickness surface in fsaverage space (in order to later compute a template Schaefer-400 GD matrix)
- `p3_computing_template_gd.sh` create a template GD matrix in Schaefer-400 space using midthickness surface in fsaverage space
- `p3_compute_mean_gd.ipynb` computes mean geodesic distances at the subject level based on fixed Schaefer-400 template distances
- `p3_cleaning_covariates.ipynb` cleans covariates to feed to p3_cleaning_final_sample.py
- `p3_cleaning_final_sample.py` gets final sample based on shared data availability (as reported in data trees)

**2. Main analyses**
- `p3_main.ipynb` computes S-A axis and features of S-A axis development, plotting Figures 1 & 2, FC strength computation
- `p3_SAexpansion_outlier_removal.ipynb` computes S-A axis expansion computation, outlier removal, plotting Figures 3 & 4
  
- `p3_bam_rcoef.R` BAM analysis - similarity to adult S-A axis
- `p3_bam_gradflip.R` BAM analysis - gradient flip
- `p3_bam_SA_expansion_net_gd_rm_outliers.R` BAM analysis - S-A axis expansion (network-level, removed outliers)

- `p3_supplementary_bam_rcoef_gradflip_relationship.R` BAM analysis - supplementary - relationship between similarity to adult S-A axis and gradient flip
- `p3_supplementary_bam_FC.R` BAM analysis - supplementary - FC strength by network
- `p3_supplementary_bam_FC_SA_bin.R` BAM analysis - supplementary - FC strength by S-A axis bin
  
- `p3_sensitivity_bam_SA_expansion_net_full_sample.R` BAM analysis - sensitivity - S-A axis expansion (network-level, full sample)
- `p3_sensitivity_bam_SA_expansion_node_rm_outliers.R` BAM analysis - sensitivity - S-A axis expansion (region-level, removed outliers)
- `p3_sensitivity_computing_BMIz.R` Computing BMI (standardized by age and sex) for sensitivity analysis
- `p3_sensitivity_bam_rcoef_BMI.R` BAM analysis - sensitivity - BMI as covariate on similarity to adult S-A axis
- `p3_sensitivity_bam_gradflip_BMI.R` BAM analysis - sensitivity - BMI as covariate on gradient flip
- `p3_sensitivity_bam_SA_expansion_net_gd_rm_outliers_BMI.R` BAM analysis - sensitivity - BMI as covariate on S-A axis expansion (network-level, removed outliers)
- `p3_sensitivity_bam_rcoef_SES.R` BAM analysis - sensitivity - SES as covariate on similarity to adult S-A axis
- `p3_sensitivity_bam_gradflip_SES.R` BAM analysis - sensitivity - SES as covariate on gradient flip
- `p3_sensitivity_bam_SA_expansion_net_gd_rm_outliers_SES.R` BAM analysis - sensitivity - SES as covariate on S-A axis expansion (network-level, removed outliers)


**3. Functions**
- `p3_myfunctions.ipynb` contains functions used for main analyses


## Data
- Sample includes 9-16 year-olds from [(ABCD study)](https://abcdstudy.org/)
- `sublists` folder includes subject lists for baseline, 2-year follow-up, and 4-year follow-up


## Support
Please address any questions about the analyses or code to [Bianca Serio](mailto:serio@cbs.mpg.de)

---

### Research poster to be presented at:
- Annual Meeting of the Organization for Human Brain Mapping (OHBM), Bordeaux 2026
- Annual meeting of the Organization for the Study of Sex Differences (OSSD), Hawaii 2026


