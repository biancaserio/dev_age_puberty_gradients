# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Development

# Age and PDS on FC strength by S-A axis bin (based on adult S-A axis)
# MAIN ANALYSIS - with OUTLIERS REMOVED 

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SET UP 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Load packages

# General
library(rstudioapi) # gets path of script  
library(R.matlab)
require(tidyverse) # Note: longCombat currently will not run if tidyverse suite is loaded. This may be fixed in the future. For now, please run longCombat before loading tidyverse.
library(plyr)  # ddply
library(ggplot2)  # plotting

library(mgcv)  # Mixed GAM Computation Vehicle with Automatic Smoothness Estimation (Generalized additive (mixed) models)
#library(nlme)  # apparently required by mgcv
library(splines)  # For testing the presence of nonlinear effects

library(gratia)  # draw function (for gamm_fit$gam results -> smooth terms)

require(lme4)  # for lmer (linear mixed effects model)
require(lmerTest)  # to obtain p-values for lmer - this actually overrides lme4's lmer() and prints the p-values for the fixed effects (which aren't present in lme4's lmer())

library(svglite)  # To save figures as vector files
library(viridis)  # to define colors for bins

### Clear environment
rm(list = ls())


#### set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/pt_02667/data/ABCD/'
datadir_local = '/data/p_02667/development/data/'
resdir = '/data/p_02667/development/results/supplementary/fc/SA_bin/'
resdir_fig = '/data/p_02667/development/results/supplementary/fc/SA_bin/figures/'



#### set directory to path of current script
setwd(codedir) 


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### mean FC strength (of top 10% connections) averaged across 10 bins (based on S-A axis loading of seed region)


## Baseline 
mat_SA_expansion_metrics_baseline <- readMat(paste(datadir_local, 'SA_expansion_metrics_baseline_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_baseline)  # print contents of matfile

fc_binned_baseline <- mat_SA_expansion_metrics_baseline$fc.binned.baseline  # shape: 3950, 1

# sub ID 
sub_ID_baseline <- mat_SA_expansion_metrics_baseline$sub.ID.baseline   # shape: 3950, 1

dim(fc_binned_baseline)
dim(sub_ID_baseline)


## 2y follow-up 
mat_SA_expansion_metrics_fu2y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu2y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu2y)  # print contents of matfile

fc_binned_fu2y <- mat_SA_expansion_metrics_fu2y$fc.binned.fu2y  # shape: 1252, 1

# sub ID 
sub_ID_fu2y <- mat_SA_expansion_metrics_fu2y$sub.ID.fu2y   # shape: 1252, 1

dim(fc_binned_fu2y)
dim(sub_ID_fu2y)


## 4y follow-up 
mat_SA_expansion_metrics_fu4y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu4y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu4y)  # print contents of matfile

fc_binned_fu4y <- mat_SA_expansion_metrics_fu4y$fc.binned.fu4y  # shape: 906, 1

# sub ID 
sub_ID_fu4y <- mat_SA_expansion_metrics_fu4y$sub.ID.fu4y   # shape: 906, 1

dim(fc_binned_fu4y)
dim(sub_ID_fu4y)


### S-A axis expansion metrics (computed at the network-level)

## Baseline 
mat_SA_expansion_metrics_baseline <- readMat(paste(datadir_local, 'SA_expansion_metrics_baseline_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_baseline)  # print contents of matfile

SA_expansion_net_baseline <- mat_SA_expansion_metrics_baseline$SA.expansion.net.baseline  # shape: 1, 3950
SA_expansion_net_baseline <- t(SA_expansion_net_baseline)  #transposing to get the shape: 3950, 1


## 2y follow-up
mat_SA_expansion_metrics_fu2y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu2y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu2y)  # print contents of matfile

SA_expansion_net_fu2y <- mat_SA_expansion_metrics_fu2y$SA.expansion.net.fu2y  # shape: 1, 1253
SA_expansion_net_fu2y <- t(SA_expansion_net_fu2y)  #transposing to get the shape: 1253, 1


## 4y follow-up
mat_SA_expansion_metrics_fu4y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu4y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu4y)  # print contents of matfile

SA_expansion_net_fu4y <- mat_SA_expansion_metrics_fu4y$SA.expansion.net.fu4y  # shape: 1, 906
SA_expansion_net_fu4y <- t(SA_expansion_net_fu4y)  #transposing to get the shape: 906, 1



### Demographics with covariates

# baseline
abcd_demo_baseline = read.csv(paste(datadir_local, 'abcd_demo_baseline_clean_rm_outliers.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

dim(abcd_demo_baseline)


# 2y
abcd_demo_fu2y = read.csv(paste(datadir_local, 'abcd_demo_fu2y_clean_rm_outliers.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

dim(abcd_demo_fu2y)


# 4y
abcd_demo_fu4y = read.csv(paste(datadir_local, 'abcd_demo_fu4y_clean_rm_outliers.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

dim(abcd_demo_fu4y)





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Data checks and computations

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


##### Making _long dataframes containing all timepoints

### Before concatenating variables and demographics across timepoints, need to check that subjects correspond across brain and demographic data

identical(as.vector(sub_ID_baseline), abcd_demo_baseline$src_subject_id_fmt)  # MUST BE TRUE to proceed
identical(as.vector(sub_ID_fu2y), abcd_demo_fu2y$src_subject_id_fmt)  # MUST BE TRUE to proceed
identical(as.vector(sub_ID_fu4y), abcd_demo_fu4y$src_subject_id_fmt)  # MUST BE TRUE to proceed


### Create lists of FC by bins and concatenate lists and demo dataframes for longitudinal analyses

#  Defining bin lists
fc_1_baseline   <- fc_binned_baseline[, 1]
fc_2_baseline  <- fc_binned_baseline[, 2]
fc_3_baseline  <- fc_binned_baseline[, 3]
fc_4_baseline  <- fc_binned_baseline[, 4]
fc_5_baseline   <- fc_binned_baseline[, 5]
fc_6_baseline  <- fc_binned_baseline[, 6]
fc_7_baseline <- fc_binned_baseline[, 7]
fc_8_baseline <- fc_binned_baseline[, 8]
fc_9_baseline <- fc_binned_baseline[, 9]
fc_10_baseline <- fc_binned_baseline[, 10]

fc_1_fu2y   <- fc_binned_fu2y[, 1]
fc_2_fu2y  <- fc_binned_fu2y[, 2]
fc_3_fu2y  <- fc_binned_fu2y[, 3]
fc_4_fu2y  <- fc_binned_fu2y[, 4]
fc_5_fu2y   <- fc_binned_fu2y[, 5]
fc_6_fu2y  <- fc_binned_fu2y[, 6]
fc_7_fu2y <- fc_binned_fu2y[, 7]
fc_8_fu2y <- fc_binned_fu2y[, 8]
fc_9_fu2y <- fc_binned_fu2y[, 9]
fc_10_fu2y <- fc_binned_fu2y[, 10]

fc_1_fu4y   <- fc_binned_fu4y[, 1]
fc_2_fu4y  <- fc_binned_fu4y[, 2]
fc_3_fu4y  <- fc_binned_fu4y[, 3]
fc_4_fu4y  <- fc_binned_fu4y[, 4]
fc_5_fu4y   <- fc_binned_fu4y[, 5]
fc_6_fu4y  <- fc_binned_fu4y[, 6]
fc_7_fu4y <- fc_binned_fu4y[, 7]
fc_8_fu4y <- fc_binned_fu4y[, 8]
fc_9_fu4y <- fc_binned_fu4y[, 9]
fc_10_fu4y <- fc_binned_fu4y[, 10]


bin1_mean_FC_long = c(fc_1_baseline, fc_1_fu2y, fc_1_fu4y)
bin2_mean_FC_long = c(fc_2_baseline, fc_2_fu2y, fc_2_fu4y)
bin3_mean_FC_long = c(fc_3_baseline, fc_3_fu2y, fc_3_fu4y)
bin4_mean_FC_long = c(fc_4_baseline, fc_4_fu2y, fc_4_fu4y)
bin5_mean_FC_long = c(fc_5_baseline, fc_5_fu2y, fc_5_fu4y)
bin6_mean_FC_long = c(fc_6_baseline, fc_6_fu2y, fc_6_fu4y)
bin7_mean_FC_long = c(fc_7_baseline, fc_7_fu2y, fc_7_fu4y)
bin8_mean_FC_long = c(fc_8_baseline, fc_8_fu2y, fc_8_fu4y)
bin9_mean_FC_long = c(fc_9_baseline, fc_9_fu2y, fc_9_fu4y)
bin10_mean_FC_long = c(fc_10_baseline, fc_10_fu2y, fc_10_fu4y)

str(bin10_mean_FC_long)



## S-A axis expansion
SA_expansion_net_long = c(SA_expansion_net_baseline, SA_expansion_net_fu2y, SA_expansion_net_fu4y)

str(SA_expansion_net_long)



# demo 
df_list <- list(abcd_demo_baseline, abcd_demo_fu2y, abcd_demo_fu4y)
abcd_demo_long <- bind_rows(df_list)

str(abcd_demo_long)






#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS -> age and pds effects on mean FC strength per network

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare the data
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Prepare the covariates

# select a subset of variables from the demo dataframe
covar_df <- abcd_demo_long %>%
  
  dplyr::select(src_subject_id_fmt, age_months, pds_p_score, sex, tot_SA, family_id, site_id, eventname) %>% 
  dplyr::rename("subject_id" = src_subject_id_fmt, "age" = age_months, "pds" = pds_p_score)


# set data type 
covar_df$age <- as.numeric(covar_df$age)
covar_df$pds <- as.numeric(covar_df$pds)
covar_df$sex <- as.factor(covar_df$sex)
covar_df$tot_SA <- as.numeric(covar_df$tot_SA)
covar_df$family_id <- as.factor(covar_df$family_id)
covar_df$site_id <- as.factor(covar_df$site_id)
covar_df$subject_id <- as.factor(covar_df$subject_id)
covar_df$eventname <- as.factor(covar_df$eventname)


### For BAM: Determining the random effects variables (interactions) needed to code sub nested in fam nested in site

any(duplicated(paste(covar_df$site_id, covar_df$family_id)))  # TRUE
# TRUE  => family labels repeat across sites (not globally unique)
# FALSE => family labels are unique per site combination

any(duplicated(paste(covar_df$site_id, covar_df$family_id, covar_df$subject_id)))  # TRUE
# TRUE  => family and subject labels repeat across sites (not globally unique)
# FALSE => family and subject labels are unique per site combination

# fyi
length(unique(covar_df$site_id))  # 22
length(unique(interaction(covar_df$site_id, covar_df$family_id, drop=TRUE)))  # 4424
length(unique(interaction(covar_df$site_id, covar_df$family_id, covar_df$subject_id, drop=TRUE)))  # 4933


## -> Given that we get TRUE: not globally unique, we need to specify 3 random effects
# s(site_id, bs = "re") +
# s(site_id.family_id, bs = "re") +
# s(site_id.family_id.subject_id, bs = "re")

## Therefore creating the relevant interaction terms:
# Nested effect of family within site
covar_df$site.family  <- interaction(covar_df$site_id, covar_df$family_id, drop = TRUE)

# Nested effect of subject within family within site
covar_df$site.family.subject <- interaction(covar_df$site_id, covar_df$family_id, covar_df$subject_id, drop = TRUE)

str(covar_df)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM analyses: Mean FC (net) ~ age + pds + sex + totSA + (1 | sub+fam+site)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### bin1

bam_FC_bin1 <- bam(
  bin1_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin1 = summary(bam_FC_bin1)

# Save summary
capture.output(summary_bam_FC_bin1, file = file.path(resdir, "summary_bam_FC_bin1_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin1)

# Extract the data for plotting
df_bam_FC_bin1_age <- smooth_estimates(bam_FC_bin1, select = "s(age)")
df_bam_FC_bin1_pds <- smooth_estimates(bam_FC_bin1, select = "s(pds)")



### bin2

bam_FC_bin2 <- bam(
  bin2_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin2 = summary(bam_FC_bin2)

# Save summary
capture.output(summary_bam_FC_bin2, file = file.path(resdir, "summary_bam_FC_bin2_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin2)

# Extract the data for plotting
df_bam_FC_bin2_age <- smooth_estimates(bam_FC_bin2, select = "s(age)")
df_bam_FC_bin2_pds <- smooth_estimates(bam_FC_bin2, select = "s(pds)")


### bin3

bam_FC_bin3 <- bam(
  bin3_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin3 = summary(bam_FC_bin3)

# Save summary
capture.output(summary_bam_FC_bin3, file = file.path(resdir, "summary_bam_FC_bin3_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin3)

# Extract the data for plotting
df_bam_FC_bin3_age <- smooth_estimates(bam_FC_bin3, select = "s(age)")
df_bam_FC_bin3_pds <- smooth_estimates(bam_FC_bin3, select = "s(pds)")


### bin4

bam_FC_bin4 <- bam(
  bin4_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin4 = summary(bam_FC_bin4)

# Save summary
capture.output(summary_bam_FC_bin4, file = file.path(resdir, "summary_bam_FC_bin4_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin4)

# Extract the data for plotting
df_bam_FC_bin4_age <- smooth_estimates(bam_FC_bin4, select = "s(age)")
df_bam_FC_bin4_pds <- smooth_estimates(bam_FC_bin4, select = "s(pds)")


### bin5

bam_FC_bin5 <- bam(
  bin5_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin5 = summary(bam_FC_bin5)

# Save summary
capture.output(summary_bam_FC_bin5, file = file.path(resdir, "summary_bam_FC_bin5_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin5)


# Extract the data for plotting
df_bam_FC_bin5_age <- smooth_estimates(bam_FC_bin5, select = "s(age)")
df_bam_FC_bin5_pds <- smooth_estimates(bam_FC_bin5, select = "s(pds)")


### bin6

bam_FC_bin6 <- bam(
  bin6_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin6 = summary(bam_FC_bin6)

# Save summary
capture.output(summary_bam_FC_bin6, file = file.path(resdir, "summary_bam_FC_bin6_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin6)

# Extract the data for plotting
df_bam_FC_bin6_age <- smooth_estimates(bam_FC_bin6, select = "s(age)")
df_bam_FC_bin6_pds <- smooth_estimates(bam_FC_bin6, select = "s(pds)")


### bin7

bam_FC_bin7 <- bam(
  bin7_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin7 = summary(bam_FC_bin7)

# Save summary
capture.output(summary_bam_FC_bin7, file = file.path(resdir, "summary_bam_FC_bin7_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin7)

# Extract the data for plotting
df_bam_FC_bin7_age <- smooth_estimates(bam_FC_bin7, select = "s(age)")
df_bam_FC_bin7_pds <- smooth_estimates(bam_FC_bin7, select = "s(pds)")


### bin8

bam_FC_bin8 <- bam(
  bin8_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin8 = summary(bam_FC_bin8)

# Save summary
capture.output(summary_bam_FC_bin8, file = file.path(resdir, "summary_bam_FC_bin8_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin8)

# Extract the data for plotting
df_bam_FC_bin8_age <- smooth_estimates(bam_FC_bin8, select = "s(age)")
df_bam_FC_bin8_pds <- smooth_estimates(bam_FC_bin8, select = "s(pds)")


### bin9

bam_FC_bin9 <- bam(
  bin9_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin9 = summary(bam_FC_bin9)

# Save summary
capture.output(summary_bam_FC_bin9, file = file.path(resdir, "summary_bam_FC_bin9_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin9)

# Extract the data for plotting
df_bam_FC_bin9_age <- smooth_estimates(bam_FC_bin9, select = "s(age)")
df_bam_FC_bin9_pds <- smooth_estimates(bam_FC_bin9, select = "s(pds)")


### bin10

bam_FC_bin10 <- bam(
  bin10_mean_FC_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_FC_bin10 = summary(bam_FC_bin10)

# Save summary
capture.output(summary_bam_FC_bin10, file = file.path(resdir, "summary_bam_FC_bin10_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_bin10)

# Extract the data for plotting
df_bam_FC_bin10_age <- smooth_estimates(bam_FC_bin10, select = "s(age)")
df_bam_FC_bin10_pds <- smooth_estimates(bam_FC_bin10, select = "s(pds)")


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of age and pds effects on Mean FC per network -> all in one
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Create a vector of 10 hex codes spanning the palette
bin_colors <- viridis(10)
bin_colors

# To see what they look like immediately:
scales::show_col(bin_colors)


# Define your custom colors 
bin_colors <- c(
  "bin1" = "#FDE725FF",
  "bin2" = "#B4DE2CFF",
  "bin3" = "#6DCD59FF",
  "bin4" = "#35B779FF",
  "bin5" = "#1F9E89FF",
  "bin6" = "#26828EFF",
  "bin7" = "#31688EFF",
  "bin8" = "#3E4A89FF",
  "bin9" = "#482878FF",
  "bin10" = "#440154FF"
)

# Define your desired order
bin_order <- c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9","bin10")

# To harmonize y-axes in panel figures (both age and pds, all networks): 
# Find the global range for smooth estimates (including confidence intervals)
smes_ylim <- range(
  c(df_bam_FC_bin1_age$.estimate - 2 * df_bam_FC_bin1_age$.se, df_bam_FC_bin1_age$.estimate + 2 * df_bam_FC_bin1_age$.se,
    df_bam_FC_bin1_pds$.estimate - 2 * df_bam_FC_bin1_pds$.se, df_bam_FC_bin1_pds$.estimate + 2 * df_bam_FC_bin1_pds$.se,
    df_bam_FC_bin2_age$.estimate - 2 * df_bam_FC_bin2_age$.se, df_bam_FC_bin2_age$.estimate + 2 * df_bam_FC_bin2_age$.se,
    df_bam_FC_bin2_pds$.estimate - 2 * df_bam_FC_bin2_pds$.se, df_bam_FC_bin2_pds$.estimate + 2 * df_bam_FC_bin2_pds$.se,
    df_bam_FC_bin3_age$.estimate - 2 * df_bam_FC_bin3_age$.se, df_bam_FC_bin3_age$.estimate + 2 * df_bam_FC_bin3_age$.se,
    df_bam_FC_bin3_pds$.estimate - 2 * df_bam_FC_bin3_pds$.se, df_bam_FC_bin3_pds$.estimate + 2 * df_bam_FC_bin3_pds$.se,
    df_bam_FC_bin4_age$.estimate - 2 * df_bam_FC_bin4_age$.se, df_bam_FC_bin4_age$.estimate + 2 * df_bam_FC_bin4_age$.se,
    df_bam_FC_bin4_pds$.estimate - 2 * df_bam_FC_bin4_pds$.se, df_bam_FC_bin4_pds$.estimate + 2 * df_bam_FC_bin4_pds$.se,
    df_bam_FC_bin5_age$.estimate - 2 * df_bam_FC_bin5_age$.se, df_bam_FC_bin5_age$.estimate + 2 * df_bam_FC_bin5_age$.se,
    df_bam_FC_bin5_pds$.estimate - 2 * df_bam_FC_bin5_pds$.se, df_bam_FC_bin5_pds$.estimate + 2 * df_bam_FC_bin5_pds$.se,
    df_bam_FC_bin6_age$.estimate - 2 * df_bam_FC_bin6_age$.se, df_bam_FC_bin6_age$.estimate + 2 * df_bam_FC_bin6_age$.se,
    df_bam_FC_bin6_pds$.estimate - 2 * df_bam_FC_bin6_pds$.se, df_bam_FC_bin6_pds$.estimate + 2 * df_bam_FC_bin6_pds$.se,
    df_bam_FC_bin7_age$.estimate - 2 * df_bam_FC_bin7_age$.se, df_bam_FC_bin7_age$.estimate + 2 * df_bam_FC_bin7_age$.se,
    df_bam_FC_bin7_pds$.estimate - 2 * df_bam_FC_bin7_pds$.se, df_bam_FC_bin7_pds$.estimate + 2 * df_bam_FC_bin7_pds$.se,
    df_bam_FC_bin8_age$.estimate - 2 * df_bam_FC_bin8_age$.se, df_bam_FC_bin8_age$.estimate + 2 * df_bam_FC_bin8_age$.se,
    df_bam_FC_bin8_pds$.estimate - 2 * df_bam_FC_bin8_pds$.se, df_bam_FC_bin8_pds$.estimate + 2 * df_bam_FC_bin8_pds$.se,
    df_bam_FC_bin9_age$.estimate - 2 * df_bam_FC_bin9_age$.se, df_bam_FC_bin9_age$.estimate + 2 * df_bam_FC_bin9_age$.se,
    df_bam_FC_bin9_pds$.estimate - 2 * df_bam_FC_bin9_pds$.se, df_bam_FC_bin9_pds$.estimate + 2 * df_bam_FC_bin9_pds$.se,
    df_bam_FC_bin10_age$.estimate - 2 * df_bam_FC_bin10_age$.se, df_bam_FC_bin10_age$.estimate + 2 * df_bam_FC_bin10_age$.se,
    df_bam_FC_bin10_pds$.estimate - 2 * df_bam_FC_bin10_pds$.se, df_bam_FC_bin10_pds$.estimate + 2 * df_bam_FC_bin10_pds$.se)
  ) 


### Age

# Combine and set factor levels
df_age_all <- bind_rows(
  df_bam_FC_bin1_age  %>% mutate(Bin = "bin1"),
  df_bam_FC_bin2_age  %>% mutate(Bin = "bin2"),
  df_bam_FC_bin3_age  %>% mutate(Bin = "bin3"),
  df_bam_FC_bin4_age  %>% mutate(Bin = "bin4"),
  df_bam_FC_bin5_age  %>% mutate(Bin = "bin5"),
  df_bam_FC_bin6_age  %>% mutate(Bin = "bin6"),
  df_bam_FC_bin7_age  %>% mutate(Bin = "bin7"),
  df_bam_FC_bin8_age  %>% mutate(Bin = "bin8"),
  df_bam_FC_bin9_age  %>% mutate(Bin = "bin9"),
  df_bam_FC_bin10_age  %>% mutate(Bin = "bin10")
) %>%
  mutate(Bin = factor(Bin, levels = bin_order))  # <-- set legend order here


# Plot

plot_smes_age_all <- ggplot(df_age_all, aes(x = age, y = .estimate, color = Bin, fill = Bin)) +
  #geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
  #alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = bin_colors) +
  scale_fill_manual(values = bin_colors) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on mean strength of\nfunctional connections",  # functional connectivity strength
    color = "Bin",
    fill = "Bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_age_all



### PDS

# Combine all smooth estimate data frames into one
df_pds_all <- bind_rows(
  df_bam_FC_bin1_pds  %>% mutate(Bin = "bin1"),
  df_bam_FC_bin2_pds  %>% mutate(Bin = "bin2"),
  df_bam_FC_bin3_pds  %>% mutate(Bin = "bin3"),
  df_bam_FC_bin4_pds  %>% mutate(Bin = "bin4"),
  df_bam_FC_bin5_pds  %>% mutate(Bin = "bin5"),
  df_bam_FC_bin6_pds  %>% mutate(Bin = "bin6"),
  df_bam_FC_bin7_pds  %>% mutate(Bin = "bin7"),
  df_bam_FC_bin8_pds  %>% mutate(Bin = "bin8"),
  df_bam_FC_bin9_pds  %>% mutate(Bin = "bin9"),
  df_bam_FC_bin10_pds  %>% mutate(Bin = "bin10")
) %>%
  mutate(Bin = factor(Bin, levels = bin_order))  # <-- set legend order here


# Plot
plot_smes_pds_all <- ggplot(df_pds_all, aes(x = pds, y = .estimate, color = Bin, fill = Bin)) +
  #geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
  #alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = bin_colors) +
  scale_fill_manual(values = bin_colors) +
  labs(
    x = "Pubertal stage (PDS score)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Bin",
    fill = "Bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_pds_all

# before saving decide if you want with or without CIs (and remove or add "_CIs" accordingly)
#ggsave(paste(resdir_fig, 'bam_smes_age_on_FC_per_network_all_rm_outliers.svg', sep = ''), plot = plot_smes_age_all, width = 6, height = 5, dpi = 300)
#ggsave(paste(resdir_fig, 'bam_smes_pds_on_FC_per_network_all_rm_outliers.svg', sep = ''), plot = plot_smes_pds_all, width = 6, height = 5, dpi = 300)


# To harmonize y-axes: Apply ylim to plots
plot_smes_age_all <- plot_smes_age_all + coord_cartesian(ylim = smes_ylim)
plot_smes_pds_all <- plot_smes_pds_all + coord_cartesian(ylim = smes_ylim)

# Combine plots in one figure
plot_smes_panel <- (plot_smes_age_all + patchwork::plot_spacer() + plot_smes_pds_all) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1)) & # Adjust 0.1 to increase/decrease the gap
  theme(
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# be careful at current CIs when plotting - only need to plot this without CIs
ggsave(paste(resdir_fig, 'bam_smes_panel_ageEpds_on_FC_per_network_all_rm_outliers.svg', sep = ''), plot = plot_smes_panel, width = 11, height = 5, dpi = 300)



### only displaying significant results (manually selecting them)

### Age


# Define your desired order (to comment out)
bin_order_sig_age <- c(
  "bin1",
  "bin2", 
  #"bin3", 
  #"bin4", 
  #"bin5", 
  #"bin6", 
  "bin7"
  #"bin8", 
  #"bin9",
  #"bin10"
  )


# Combine and set factor levels (to comment out)
df_age_all_sig <- bind_rows(
  df_bam_FC_bin1_age  %>% mutate(Bin = "bin1"),
  df_bam_FC_bin2_age  %>% mutate(Bin = "bin2"),
  #df_bam_FC_bin3_age  %>% mutate(Bin = "bin3"),
  #df_bam_FC_bin4_age  %>% mutate(Bin = "bin4"),
  #df_bam_FC_bin5_age  %>% mutate(Bin = "bin5"),
  #df_bam_FC_bin6_age  %>% mutate(Bin = "bin6"),
  df_bam_FC_bin7_age  %>% mutate(Bin = "bin7")
  #df_bam_FC_bin8_age  %>% mutate(Bin = "bin8"),
  #df_bam_FC_bin9_age  %>% mutate(Bin = "bin9"),
  #df_bam_FC_bin10_age  %>% mutate(Bin = "bin10")
) %>%
  mutate(Bin = factor(Bin, levels = bin_order))  # <-- set legend order here



# Plot
plot_smes_age_all_sig <- ggplot(df_age_all_sig, aes(x = age, y = .estimate, color = Bin, fill = Bin)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = bin_colors) +
  scale_fill_manual(values = bin_colors) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Bin",
    fill = "Bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_age_all_sig



### PDS

# Define your desired order (to comment out)
bin_order_sig_pds <- c(
  "bin1",
  #"bin2", 
  #"bin3", 
  #"bin4", 
  "bin5"
  #"bin6", 
  #"bin7", 
  #"bin8", 
  #"bin9",
  #"bin10"
  )


# Combine all smooth estimate data frames into one (to comment out)
df_pds_all_sig <- bind_rows(
  df_bam_FC_bin1_pds  %>% mutate(Bin = "bin1"),
  #df_bam_FC_bin2_pds  %>% mutate(Bin = "bin2"),
  #df_bam_FC_bin3_pds  %>% mutate(Bin = "bin3"),
  #df_bam_FC_bin4_pds  %>% mutate(Bin = "bin4"),
  df_bam_FC_bin5_pds  %>% mutate(Bin = "bin5")
  #df_bam_FC_bin6_pds  %>% mutate(Bin = "bin6"),
  #df_bam_FC_bin7_pds  %>% mutate(Bin = "bin7"),
  #df_bam_FC_bin8_pds  %>% mutate(Bin = "bin8"),
  #df_bam_FC_bin9_pds  %>% mutate(Bin = "bin9"),
  #df_bam_FC_bin10_pds  %>% mutate(Bin = "bin10")
) %>%
  mutate(Bin = factor(Bin, levels = bin_order))  # <-- set legend order here



# Plot
plot_smes_pds_all_sig <- ggplot(df_pds_all_sig, aes(x = pds, y = .estimate, color = Bin, fill = Bin)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = bin_colors) +
  scale_fill_manual(values = bin_colors) +
  labs(
    x = "Pubertal stage (PDS score)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Bin",
    fill = "Bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_pds_all_sig

# before saving decide if you want with or without CIs (and remove or add "_CIs" accordingly)
#ggsave(paste(resdir_fig, 'bam_smes_age_on_FC_per_network_all_sig_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_age_all_sig, width = 6, height = 5, dpi = 300)
#ggsave(paste(resdir_fig, 'bam_smes_pds_on_FC_per_network_all_sig_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_pds_all_sig, width = 6, height = 5, dpi = 300)


# To harmonize y-axes: Apply ylim to plots
plot_smes_age_all_sig <- plot_smes_age_all_sig + coord_cartesian(ylim = smes_ylim)
plot_smes_pds_all_sig <- plot_smes_pds_all_sig + coord_cartesian(ylim = smes_ylim)

# Combine plots in one figure
plot_smes_panel_all_sig <- (plot_smes_age_all_sig + patchwork::plot_spacer() + plot_smes_pds_all_sig) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1)) & # Adjust 0.1 to increase/decrease the gap
  theme(
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# be careful at current CIs when plotting - only need to plot this with CIs
ggsave(paste(resdir_fig, 'bam_smes_panel_ageEpds_on_FC_per_network_all_sig_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_panel_all_sig, width = 11, height = 5, dpi = 300)





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS -> mean FCs effects on S-A axis expansion

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Bin-level (10 Bins)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM analysis: S-A axis expansion ~ age + pds + FCs(of each bin) + sex + totSA + (1 | sub+fam+site)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bam_SA_exp_FCs_all_bin <- bam(
  SA_expansion_net_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    s(bin1_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin2_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin3_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin4_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin5_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin6_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin7_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin8_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin9_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(bin10_mean_FC_long, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE   # allows automatic penalization of unneeded terms
)

summary_bam_SA_exp_FCs_all_bin_rm_outliers = summary(bam_SA_exp_FCs_all_bin)

# Save summary
capture.output(summary_bam_SA_exp_FCs_all_bin_rm_outliers, file = file.path(resdir, "summary_bam_SA_exp_FCs_all_bin_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_FCs_all_bin)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of effects (on SA exp) of mean FCs per network -> all in one
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define your custom colors 
bin_colors <- c(
  "bin1" = "#FDE725FF",
  "bin2" = "#B4DE2CFF",
  "bin3" = "#6DCD59FF",
  "bin4" = "#35B779FF",
  "bin5" = "#1F9E89FF",
  "bin6" = "#26828EFF",
  "bin7" = "#31688EFF",
  "bin8" = "#3E4A89FF",
  "bin9" = "#482878FF",
  "bin10" = "#440154FF"
)

# Define your desired order
bin_order <- c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9","bin10")

# Extract the data for more control
df_age <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(age)")
df_pds <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(pds)")

df_bin1 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin1_mean_FC_long)")
df_bin2 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin2_mean_FC_long)")
df_bin3 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin3_mean_FC_long)")
df_bin4 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin4_mean_FC_long)")
df_bin5 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin5_mean_FC_long)")
df_bin6 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin6_mean_FC_long)")
df_bin7 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin7_mean_FC_long)")
df_bin8 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin8_mean_FC_long)")
df_bin9 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin9_mean_FC_long)")
df_bin10 <- smooth_estimates(bam_SA_exp_FCs_all_bin, select = "s(bin10_mean_FC_long)")



# Combine and set factor levels (by also renaming the predictor column to "mean_FC" for all networks before binding)

df_all_bin <- bind_rows(
  df_bin1  %>% dplyr::rename(mean_FC = bin1_mean_FC_long)   %>% mutate(Bin = "bin1"),
  df_bin2  %>% dplyr::rename(mean_FC = bin2_mean_FC_long)   %>% mutate(Bin = "bin2"),
  df_bin3  %>% dplyr::rename(mean_FC = bin3_mean_FC_long)   %>% mutate(Bin = "bin3"),
  df_bin4  %>% dplyr::rename(mean_FC = bin4_mean_FC_long)   %>% mutate(Bin = "bin4"),
  df_bin5  %>% dplyr::rename(mean_FC = bin5_mean_FC_long)   %>% mutate(Bin = "bin5"),
  df_bin6  %>% dplyr::rename(mean_FC = bin6_mean_FC_long)   %>% mutate(Bin = "bin6"),
  df_bin7  %>% dplyr::rename(mean_FC = bin7_mean_FC_long)   %>% mutate(Bin = "bin7"),
  df_bin8  %>% dplyr::rename(mean_FC = bin8_mean_FC_long)   %>% mutate(Bin = "bin8"),
  df_bin9  %>% dplyr::rename(mean_FC = bin9_mean_FC_long)   %>% mutate(Bin = "bin9"),
  df_bin10  %>% dplyr::rename(mean_FC = bin10_mean_FC_long)   %>% mutate(Bin = "bin10")
) %>%
  mutate(Bin = factor(Bin, levels = bin_order))  # set legend order here


# Plot (effects of FCs of all networks on S-A axis expansion)
plot_smes_all_bin <- ggplot(df_all_bin, aes(x = mean_FC, y = .estimate, color = Bin, fill = Bin)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = bin_colors) +
  scale_fill_manual(values = bin_colors) +
  labs(
    x = "Mean strength of functional connections",
    y = "Partial effect on S-A axis expansion",
    color = "Bin",
    fill = "Bin"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_all_bin


# before saving decide if you want with or without CIs (and remove or add "_CIs" accordingly)
ggsave(paste(resdir_fig, 'bam_smes_all_bins_FCs_on_SAexp_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_all_bin, width = 6, height = 5, dpi = 300)




plot_smes_age = ggplot(df_age, aes(x = age, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "honeydew3", alpha = 0.25) +
  geom_line(color = "honeydew4", linewidth = 1.2) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on S-A axis expansion",
    #title = unique(df_age$.smooth)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    #axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

plot_smes_pds = ggplot(df_pds, aes(x = pds, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "thistle", alpha = 0.25) +
  geom_line(color = "thistle4", linewidth = 1.2) +
  labs(
    x = "Pubertal stage (PDS score)",
    y = "Partial effect on S-A axis expansion",
    #title = unique(df_pds$.smooth)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    #axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )


ggsave(paste(resdir_fig, 'bam_smes_all_bins_age_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_all_bins_pds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)






