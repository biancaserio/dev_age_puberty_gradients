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


### Clear environment
rm(list = ls())


#### set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/pt_02667/data/ABCD/'
datadir_local = '/data/p_02667/development/data/'
resdir = '/data/p_02667/development/results/supplementary/fc/'
resdir_fig = '/data/p_02667/development/results/supplementary/fc/figures/'



#### set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### mean FC strength (of top 10% connections) averaged across network of seed region


## Baseline 
mat_SA_expansion_metrics_baseline <- readMat(paste(datadir_local, 'SA_expansion_metrics_baseline_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_baseline)  # print contents of matfile

fc_network_baseline <- mat_SA_expansion_metrics_baseline$fc.network.baseline  # shape: 3950, 1

# sub ID 
sub_ID_baseline <- mat_SA_expansion_metrics_baseline$sub.ID.baseline   # shape: 3950, 1

dim(fc_network_baseline)
dim(sub_ID_baseline)


## 2y follow-up 
mat_SA_expansion_metrics_fu2y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu2y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu2y)  # print contents of matfile

fc_network_fu2y <- mat_SA_expansion_metrics_fu2y$fc.network.fu2y  # shape: 1252, 1

# sub ID 
sub_ID_fu2y <- mat_SA_expansion_metrics_fu2y$sub.ID.fu2y   # shape: 1252, 1

dim(fc_network_fu2y)
dim(sub_ID_fu2y)


## 4y follow-up 
mat_SA_expansion_metrics_fu4y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu4y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu4y)  # print contents of matfile

fc_network_fu4y <- mat_SA_expansion_metrics_fu4y$fc.network.fu4y  # shape: 906, 1

# sub ID 
sub_ID_fu4y <- mat_SA_expansion_metrics_fu4y$sub.ID.fu4y   # shape: 906, 1

dim(fc_network_fu4y)
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


### Create lists of FC by network and concatenate lists and demo dataframes for longitudinal analyses

#  Defining network lists
fc_V_baseline   <- fc_network_baseline[, 1]
fc_SM_baseline  <- fc_network_baseline[, 2]
fc_DA_baseline  <- fc_network_baseline[, 3]
fc_VA_baseline  <- fc_network_baseline[, 4]
fc_L_baseline   <- fc_network_baseline[, 5]
fc_FP_baseline  <- fc_network_baseline[, 6]
fc_DMN_baseline <- fc_network_baseline[, 7]

fc_V_fu2y   <- fc_network_fu2y[, 1]
fc_SM_fu2y  <- fc_network_fu2y[, 2]
fc_DA_fu2y  <- fc_network_fu2y[, 3]
fc_VA_fu2y  <- fc_network_fu2y[, 4]
fc_L_fu2y   <- fc_network_fu2y[, 5]
fc_FP_fu2y  <- fc_network_fu2y[, 6]
fc_DMN_fu2y <- fc_network_fu2y[, 7]

fc_V_fu4y   <- fc_network_fu4y[, 1]
fc_SM_fu4y  <- fc_network_fu4y[, 2]
fc_DA_fu4y  <- fc_network_fu4y[, 3]
fc_VA_fu4y  <- fc_network_fu4y[, 4]
fc_L_fu4y   <- fc_network_fu4y[, 5]
fc_FP_fu4y  <- fc_network_fu4y[, 6]
fc_DMN_fu4y <- fc_network_fu4y[, 7]


V_mean_FC_long = c(fc_V_baseline, fc_V_fu2y, fc_V_fu4y)
SM_mean_FC_long = c(fc_SM_baseline, fc_SM_fu2y, fc_SM_fu4y)
DA_mean_FC_long = c(fc_DA_baseline, fc_DA_fu2y, fc_DA_fu4y)
VA_mean_FC_long = c(fc_VA_baseline, fc_VA_fu2y, fc_VA_fu4y)
L_mean_FC_long = c(fc_L_baseline, fc_L_fu2y, fc_L_fu4y)
FP_mean_FC_long = c(fc_FP_baseline, fc_FP_fu2y, fc_FP_fu4y)
DMN_mean_FC_long = c(fc_DMN_baseline, fc_DMN_fu2y, fc_DMN_fu4y)

str(V_mean_FC_long)


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

### V

bam_FC_V <- bam(
  V_mean_FC_long ~
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

summary_bam_FC_V = summary(bam_FC_V)

# Save summary
capture.output(summary_bam_FC_V, file = file.path(resdir, "summary_bam_FC_V_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_V)

# Extract the data for plotting
df_bam_FC_V_age <- smooth_estimates(bam_FC_V, select = "s(age)")
df_bam_FC_V_pds <- smooth_estimates(bam_FC_V, select = "s(pds)")



### SM

bam_FC_SM <- bam(
  SM_mean_FC_long ~
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

summary_bam_FC_SM = summary(bam_FC_SM)

# Save summary
capture.output(summary_bam_FC_SM, file = file.path(resdir, "summary_bam_FC_SM_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_SM)

# Extract the data for plotting
df_bam_FC_SM_age <- smooth_estimates(bam_FC_SM, select = "s(age)")
df_bam_FC_SM_pds <- smooth_estimates(bam_FC_SM, select = "s(pds)")


### DA

bam_FC_DA <- bam(
  DA_mean_FC_long ~
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

summary_bam_FC_DA = summary(bam_FC_DA)

# Save summary
capture.output(summary_bam_FC_DA, file = file.path(resdir, "summary_bam_FC_DA_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_DA)

# Extract the data for plotting
df_bam_FC_DA_age <- smooth_estimates(bam_FC_DA, select = "s(age)")
df_bam_FC_DA_pds <- smooth_estimates(bam_FC_DA, select = "s(pds)")


### VA

bam_FC_VA <- bam(
  VA_mean_FC_long ~
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

summary_bam_FC_VA = summary(bam_FC_VA)

# Save summary
capture.output(summary_bam_FC_VA, file = file.path(resdir, "summary_bam_FC_VA_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_VA)

# Extract the data for plotting
df_bam_FC_VA_age <- smooth_estimates(bam_FC_VA, select = "s(age)")
df_bam_FC_VA_pds <- smooth_estimates(bam_FC_VA, select = "s(pds)")


### L

bam_FC_L <- bam(
  L_mean_FC_long ~
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

summary_bam_FC_L = summary(bam_FC_L)

# Save summary
capture.output(summary_bam_FC_L, file = file.path(resdir, "summary_bam_FC_L_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_L)


# Extract the data for plotting
df_bam_FC_L_age <- smooth_estimates(bam_FC_L, select = "s(age)")
df_bam_FC_L_pds <- smooth_estimates(bam_FC_L, select = "s(pds)")


### FP

bam_FC_FP <- bam(
  FP_mean_FC_long ~
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

summary_bam_FC_FP = summary(bam_FC_FP)

# Save summary
capture.output(summary_bam_FC_FP, file = file.path(resdir, "summary_bam_FC_FP_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_FP)

# Extract the data for plotting
df_bam_FC_FP_age <- smooth_estimates(bam_FC_FP, select = "s(age)")
df_bam_FC_FP_pds <- smooth_estimates(bam_FC_FP, select = "s(pds)")


### DMN

bam_FC_DMN <- bam(
  DMN_mean_FC_long ~
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

summary_bam_FC_DMN = summary(bam_FC_DMN)

# Save summary
capture.output(summary_bam_FC_DMN, file = file.path(resdir, "summary_bam_FC_DMN_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_FC_DMN)

# Extract the data for plotting
df_bam_FC_DMN_age <- smooth_estimates(bam_FC_DMN, select = "s(age)")
df_bam_FC_DMN_pds <- smooth_estimates(bam_FC_DMN, select = "s(pds)")



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of age and pds effects on Mean FC per network -> all in one
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Define your custom colors
network_colors <- c(
  "V" = "darkorchid",
  "SM" = "steelblue",
  "DA" = "forestgreen",
  "VA" = "orchid",
  "L" = "lemonchiffon",
  "FP" = "orange",
  "DMN" = "indianred"
)

# Define your desired order
network_order <- c("V", "SM", "DA", "VA", "L", "FP", "DMN")

# To harmonize y-axes in panel figures (both age and pds, all networks): 
# Find the global range for smooth estimates (including confidence intervals)
smes_ylim <- range(
  c(df_bam_FC_V_age$.estimate - 2 * df_bam_FC_V_age$.se, df_bam_FC_V_age$.estimate + 2 * df_bam_FC_V_age$.se,
    df_bam_FC_V_pds$.estimate - 2 * df_bam_FC_V_pds$.se, df_bam_FC_V_pds$.estimate + 2 * df_bam_FC_V_pds$.se,
    df_bam_FC_SM_age$.estimate - 2 * df_bam_FC_SM_age$.se, df_bam_FC_SM_age$.estimate + 2 * df_bam_FC_SM_age$.se,
    df_bam_FC_SM_pds$.estimate - 2 * df_bam_FC_SM_pds$.se, df_bam_FC_SM_pds$.estimate + 2 * df_bam_FC_SM_pds$.se,
    df_bam_FC_DA_age$.estimate - 2 * df_bam_FC_DA_age$.se, df_bam_FC_DA_age$.estimate + 2 * df_bam_FC_DA_age$.se,
    df_bam_FC_DA_pds$.estimate - 2 * df_bam_FC_DA_pds$.se, df_bam_FC_DA_pds$.estimate + 2 * df_bam_FC_DA_pds$.se,
    df_bam_FC_VA_age$.estimate - 2 * df_bam_FC_VA_age$.se, df_bam_FC_VA_age$.estimate + 2 * df_bam_FC_VA_age$.se,
    df_bam_FC_VA_pds$.estimate - 2 * df_bam_FC_VA_pds$.se, df_bam_FC_VA_pds$.estimate + 2 * df_bam_FC_VA_pds$.se,
    df_bam_FC_L_age$.estimate - 2 * df_bam_FC_L_age$.se, df_bam_FC_L_age$.estimate + 2 * df_bam_FC_L_age$.se,
    df_bam_FC_L_pds$.estimate - 2 * df_bam_FC_L_pds$.se, df_bam_FC_L_pds$.estimate + 2 * df_bam_FC_L_pds$.se,
    df_bam_FC_FP_age$.estimate - 2 * df_bam_FC_FP_age$.se, df_bam_FC_FP_age$.estimate + 2 * df_bam_FC_FP_age$.se,
    df_bam_FC_FP_pds$.estimate - 2 * df_bam_FC_FP_pds$.se, df_bam_FC_FP_pds$.estimate + 2 * df_bam_FC_FP_pds$.se,
    df_bam_FC_DMN_age$.estimate - 2 * df_bam_FC_DMN_age$.se, df_bam_FC_DMN_age$.estimate + 2 * df_bam_FC_DMN_age$.se,
    df_bam_FC_DMN_pds$.estimate - 2 * df_bam_FC_DMN_pds$.se, df_bam_FC_DMN_pds$.estimate + 2 * df_bam_FC_DMN_pds$.se)
)


### Age

# Combine and set factor levels
df_age_all <- bind_rows(
  df_bam_FC_V_age  %>% mutate(Network = "V"),
  df_bam_FC_SM_age %>% mutate(Network = "SM"),
  df_bam_FC_DA_age %>% mutate(Network = "DA"),
  df_bam_FC_VA_age %>% mutate(Network = "VA"),
  df_bam_FC_L_age  %>% mutate(Network = "L"),
  df_bam_FC_FP_age %>% mutate(Network = "FP"),
  df_bam_FC_DMN_age %>% mutate(Network = "DMN")
) %>%
  mutate(Network = factor(Network, levels = network_order))  # <-- set legend order here


# Plot

plot_smes_age_all <- ggplot(df_age_all, aes(x = age, y = .estimate, color = Network, fill = Network)) +
  #geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
  #alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = network_colors) +
  scale_fill_manual(values = network_colors) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Network",
    fill = "Network"
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
  df_bam_FC_V_pds  %>% mutate(Network = "V"),
  df_bam_FC_SM_pds %>% mutate(Network = "SM"),
  df_bam_FC_DA_pds %>% mutate(Network = "DA"),
  df_bam_FC_VA_pds %>% mutate(Network = "VA"),
  df_bam_FC_L_pds  %>% mutate(Network = "L"),
  df_bam_FC_FP_pds %>% mutate(Network = "FP"),
  df_bam_FC_DMN_pds %>% mutate(Network = "DMN")
) %>%
  mutate(Network = factor(Network, levels = network_order))  # <-- set legend order here


# Plot
plot_smes_pds_all <- ggplot(df_pds_all, aes(x = pds, y = .estimate, color = Network, fill = Network)) +
  #geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
  #alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = network_colors) +
  scale_fill_manual(values = network_colors) +
  labs(
    x = "Pubertal stage (PDS score)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Network",
    fill = "Network"
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

# Define your desired order
network_order_sig_age <- c(
  #"V", 
  #"SM", 
  #"DA", 
  #"VA", 
  "L" 
  #"FP", 
  #"DMN"
)


# Combine and set factor levels
df_age_all_sig <- bind_rows(
  #df_bam_FC_V_age  %>% mutate(Network = "V"),
  #df_bam_FC_SM_age %>% mutate(Network = "SM"),
  #df_bam_FC_DA_age %>% mutate(Network = "DA"),
  #df_bam_FC_VA_age %>% mutate(Network = "VA"),
  df_bam_FC_L_age  %>% mutate(Network = "L"),
  #df_bam_FC_FP_age %>% mutate(Network = "FP"),
  #df_bam_FC_DMN_age %>% mutate(Network = "DMN")
) %>%
  mutate(Network = factor(Network, levels = network_order_sig_age))  # <-- set legend order here


# Plot
plot_smes_age_all_sig <- ggplot(df_age_all_sig, aes(x = age, y = .estimate, color = Network, fill = Network)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = network_colors) +
  scale_fill_manual(values = network_colors) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Network",
    fill = "Network"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_age_all_sig



### PDS

# Define your desired order
network_order_sig_pds <- c(
  #"V", 
  #"SM", 
  #"DA", 
  #"VA", 
  #"L", 
  "FP"
  #"DMN"
)

# Combine all smooth estimate data frames into one
df_pds_all_sig <- bind_rows(
  #df_bam_FC_V_pds  %>% mutate(Network = "V"),
  #df_bam_FC_SM_pds %>% mutate(Network = "SM"),
  #df_bam_FC_DA_pds %>% mutate(Network = "DA"),
  #df_bam_FC_VA_pds %>% mutate(Network = "VA"),
  #df_bam_FC_L_pds  %>% mutate(Network = "L")
  df_bam_FC_FP_pds %>% mutate(Network = "FP"),
  #df_bam_FC_DMN_pds %>% mutate(Network = "DMN")
) %>%
  mutate(Network = factor(Network, levels = network_order_sig_pds))  # <-- set legend order here


# Plot
plot_smes_pds_all_sig <- ggplot(df_pds_all_sig, aes(x = pds, y = .estimate, color = Network, fill = Network)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = network_colors) +
  scale_fill_manual(values = network_colors) +
  labs(
    x = "Pubertal stage (PDS score)",
    y = "Partial effect on mean strength of\nfunctional connections",
    color = "Network",
    fill = "Network"
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





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of age and pds effects (predicted marginal trajectories) on Mean FC -> one plot per network (and effect)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



### 1. Build a prediction grid (it will be used for all networks)

# for Age

age_grid<- seq(
  min(covar_df$age, na.rm=TRUE),
  max(covar_df$age, na.rm=TRUE),
  length.out=200
)

newdata_age <- expand.grid(
  age = age_grid,
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# for PDS
pds_grid <- seq(
  min(covar_df$pds, na.rm = TRUE),
  max(covar_df$pds, na.rm = TRUE),
  length.out = 200
)

newdata_pds <- expand.grid(
  pds = pds_grid,
  sex = levels(covar_df$sex),
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold totSA at mean
  site_id = covar_df$site_id[1],
  site.family = covar_df$site.family[1],
  site.family.subject = covar_df$site.family.subject[1]
)




### V

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_V, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_V, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the visual network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_V_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the visual network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_V_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)



### SM

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_SM, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_SM, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the somatomotor network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_SM_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the somatomotor network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_SM_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)


### DA

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_DA, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_DA, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the dorsal attention network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_DA_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the dorsal attention network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_DA_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)


### VA

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_VA, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_VA, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the ventral attention network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_VA_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the ventral attention network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_VA_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)


### L

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_L, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_L, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the limbic network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_L_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the limbic network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_L_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)


### FP

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_FP, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_FP, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the frontoparietal network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_FP_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the frontoparietal network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_FP_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)


### DMN

## 2. Generate the predictions

# Age

smooth_preds_age <- fitted_values(
  bam_FC_DMN, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# PDS

smooth_preds_pds <- fitted_values(
  bam_FC_DMN, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


## 3. Plot marginal trajectories

# Age

plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Age (in months)", y = "Mean strength of functional connections\nin the default-mode network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_age_on_FC_DMN_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


# PDS

plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Mean strength of functional connections\nin the default-mode network", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_pds_on_FC_DMN_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)







#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS -> mean FCs effects on S-A axis expansion

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Network-level (Yeo 7)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM analysis: S-A axis expansion ~ age + pds + FCs(of each network) + sex + totSA + (1 | sub+fam+site)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bam_SA_exp_FCs_all_net <- bam(
  SA_expansion_net_long ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    s(V_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(SM_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(DA_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(VA_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(L_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(FP_mean_FC_long, k = 10, fx = FALSE) +  # default
    s(DMN_mean_FC_long, k = 10, fx = FALSE) +  # default
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

summary_bam_SA_exp_FCs_all_net_rm_outliers = summary(bam_SA_exp_FCs_all_net)

# Save summary
capture.output(summary_bam_SA_exp_FCs_all_net_rm_outliers, file = file.path(resdir, "summary_bam_SA_exp_FCs_all_net_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_FCs_all_net)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of effects (on SA exp) of mean FCs per network -> all in one
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define your custom colors
network_colors <- c(
  "V" = "darkorchid",
  "SM" = "steelblue",
  "DA" = "forestgreen",
  "VA" = "orchid",
  "L" = "lemonchiffon",
  "FP" = "orange",
  "DMN" = "indianred"
)

# Define your desired order
network_order <- c("V", "SM", "DA", "VA", "L", "FP", "DMN")


# Extract the data for more control
df_age <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(age)")
df_pds <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(pds)")

df_V <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(V_mean_FC_long)")
df_SM <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(SM_mean_FC_long)")
df_DA <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(DA_mean_FC_long)")
df_VA <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(VA_mean_FC_long)")
df_L <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(L_mean_FC_long)")
df_FP <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(FP_mean_FC_long)")
df_DMN <- smooth_estimates(bam_SA_exp_FCs_all_net, select = "s(DMN_mean_FC_long)")


# Combine and set factor levels (by also renaming the predictor column to "mean_FC" for all networks before binding)

df_all_net <- bind_rows(
  #df_V  %>% dplyr::rename(mean_FC = V_mean_FC_long)   %>% mutate(Network = "V"),
  df_SM %>% dplyr::rename(mean_FC = SM_mean_FC_long)  %>% mutate(Network = "SM"),
  df_DA %>% dplyr::rename(mean_FC = DA_mean_FC_long)  %>% mutate(Network = "DA"),
  df_VA %>% dplyr::rename(mean_FC = VA_mean_FC_long)  %>% mutate(Network = "VA"),
  #df_L  %>% dplyr::rename(mean_FC = L_mean_FC_long)   %>% mutate(Network = "L"),
  df_FP %>% dplyr::rename(mean_FC = FP_mean_FC_long)  %>% mutate(Network = "FP"),
  df_DMN %>% dplyr::rename(mean_FC = DMN_mean_FC_long) %>% mutate(Network = "DMN")
) %>%
  mutate(Network = factor(Network, levels = network_order))  # set legend order here


# Plot (effects of FCs of all networks on S-A axis expansion)
plot_smes_all_net <- ggplot(df_all_net, aes(x = mean_FC, y = .estimate, color = Network, fill = Network)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = network_colors) +
  scale_fill_manual(values = network_colors) +
  labs(
    x = "Mean strength of functional connections",
    y = "Partial effect on S-A axis expansion",
    color = "Network",
    fill = "Network"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_all_net


# before saving decide if you want with or without CIs (and remove or add "_CIs" accordingly) + if I only want to plot significance
#ggsave(paste(resdir_fig, 'bam_smes_all_networks_FCs_on_SAexp_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_all_net, width = 6, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_all_networks_FCs_on_SAexp_rm_outliers_CIs_onlysig.svg', sep = ''), plot = plot_smes_all_net, width = 6, height = 5, dpi = 300)




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


ggsave(paste(resdir_fig, 'bam_smes_all_networks_age_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_all_networks_pds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)






# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of effects of Mean FCs by network (predicted marginal trajectories) on S-A axis expansion -> one plot per network
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## V

# 1. Build prediction grid

V_grid<- seq(
  min(V_mean_FC_long, na.rm=TRUE),
  max(V_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_V <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = V_grid,  # these are the predictions that we're plotting
  SM_mean_FC_long = mean(SM_mean_FC_long, na.rm = TRUE),  # hold at mean
  DA_mean_FC_long = mean(DA_mean_FC_long, na.rm = TRUE),  # hold at mean
  VA_mean_FC_long = mean(VA_mean_FC_long, na.rm = TRUE),  # hold at mean
  L_mean_FC_long = mean(L_mean_FC_long, na.rm = TRUE),  # hold at mean
  FP_mean_FC_long = mean(FP_mean_FC_long, na.rm = TRUE),  # hold at mean
  DMN_mean_FC_long = mean(DMN_mean_FC_long, na.rm = TRUE),  # hold at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_V <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_V, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_V <- ggplot(
  smooth_preds_V, 
  aes(x = V_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the visual network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_V_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_V, width = 6, height = 5, dpi = 300)



## SM

# 1. Build prediction grid

SM_grid<- seq(
  min(SM_mean_FC_long, na.rm=TRUE),
  max(SM_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_SM <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = mean(V_mean_FC_long, na.rm = TRUE),  # hold at mean
  SM_mean_FC_long = SM_grid,  # these are the predictions that we're plotting
  DA_mean_FC_long = mean(DA_mean_FC_long, na.rm = TRUE),  # hold at mean
  VA_mean_FC_long = mean(VA_mean_FC_long, na.rm = TRUE),  # hold at mean
  L_mean_FC_long = mean(L_mean_FC_long, na.rm = TRUE),  # hold at mean
  FP_mean_FC_long = mean(FP_mean_FC_long, na.rm = TRUE),  # hold at mean
  DMN_mean_FC_long = mean(DMN_mean_FC_long, na.rm = TRUE),  # hold at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_SM <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_SM, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_SM <- ggplot(
  smooth_preds_SM, 
  aes(x = SM_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the somatomotor network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_SM_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_SM, width = 6, height = 5, dpi = 300)



## DA

# 1. Build prediction grid

DA_grid<- seq(
  min(DA_mean_FC_long, na.rm=TRUE),
  max(DA_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_DA <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = mean(V_mean_FC_long, na.rm = TRUE),  # hold at mean
  SM_mean_FC_long = mean(SM_mean_FC_long, na.rm = TRUE),  # hold at mean
  DA_mean_FC_long = DA_grid,  # these are the predictions that we're plotting
  VA_mean_FC_long = mean(VA_mean_FC_long, na.rm = TRUE),  # hold at mean
  L_mean_FC_long = mean(L_mean_FC_long, na.rm = TRUE),  # hold at mean
  FP_mean_FC_long = mean(FP_mean_FC_long, na.rm = TRUE),  # hold at mean
  DMN_mean_FC_long = mean(DMN_mean_FC_long, na.rm = TRUE),  # hold at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_DA <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_DA, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_DA <- ggplot(
  smooth_preds_DA, 
  aes(x = DA_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the dorsal attention network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_DA_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_DA, width = 6, height = 5, dpi = 300)



## VA

# 1. Build prediction grid

VA_grid<- seq(
  min(VA_mean_FC_long, na.rm=TRUE),
  max(VA_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_VA <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = mean(V_mean_FC_long, na.rm = TRUE),  # hold at mean
  SM_mean_FC_long = mean(SM_mean_FC_long, na.rm = TRUE),  # hold at mean
  DA_mean_FC_long = mean(DA_mean_FC_long, na.rm = TRUE),  # hold at mean
  VA_mean_FC_long = VA_grid,  # these are the predictions that we're plotting
  L_mean_FC_long = mean(L_mean_FC_long, na.rm = TRUE),  # hold at mean
  FP_mean_FC_long = mean(FP_mean_FC_long, na.rm = TRUE),  # hold at mean
  DMN_mean_FC_long = mean(DMN_mean_FC_long, na.rm = TRUE),  # hold at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_VA <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_VA, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_VA <- ggplot(
  smooth_preds_VA, 
  aes(x = VA_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the ventral attention network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_VA_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_VA, width = 6, height = 5, dpi = 300)



## L

# 1. Build prediction grid

L_grid<- seq(
  min(L_mean_FC_long, na.rm=TRUE),
  max(L_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_L <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = mean(V_mean_FC_long, na.rm = TRUE),  # hold at mean
  SM_mean_FC_long = mean(SM_mean_FC_long, na.rm = TRUE),  # hold at mean
  DA_mean_FC_long = mean(DA_mean_FC_long, na.rm = TRUE),  # hold at mean
  VA_mean_FC_long = mean(VA_mean_FC_long, na.rm = TRUE),  # hold at mean
  L_mean_FC_long = L_grid,  # these are the predictions that we're plotting
  FP_mean_FC_long = mean(FP_mean_FC_long, na.rm = TRUE),  # hold at mean
  DMN_mean_FC_long = mean(DMN_mean_FC_long, na.rm = TRUE),  # hold at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_L <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_L, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_L <- ggplot(
  smooth_preds_L, 
  aes(x = L_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the limbic network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_L_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_L, width = 6, height = 5, dpi = 300)



## FP

# 1. Build prediction grid

FP_grid<- seq(
  min(FP_mean_FC_long, na.rm=TRUE),
  max(FP_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_FP <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = mean(V_mean_FC_long, na.rm = TRUE),  # hold at mean
  SM_mean_FC_long = mean(SM_mean_FC_long, na.rm = TRUE),  # hold at mean
  DA_mean_FC_long = mean(DA_mean_FC_long, na.rm = TRUE),  # hold at mean
  VA_mean_FC_long = mean(VA_mean_FC_long, na.rm = TRUE),  # hold at mean
  L_mean_FC_long = mean(L_mean_FC_long, na.rm = TRUE),  # hold at mean
  FP_mean_FC_long = FP_grid,  # these are the predictions that we're plotting
  DMN_mean_FC_long = mean(DMN_mean_FC_long, na.rm = TRUE),  # hold at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_FP <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_FP, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_FP <- ggplot(
  smooth_preds_FP, 
  aes(x = FP_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the frontoparietal network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_FP_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_FP, width = 6, height = 5, dpi = 300)



## DMN

# 1. Build prediction grid

DMN_grid<- seq(
  min(DMN_mean_FC_long, na.rm=TRUE),
  max(DMN_mean_FC_long, na.rm=TRUE),
  length.out=200
)

newdata_DMN <- expand.grid(
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  V_mean_FC_long = mean(V_mean_FC_long, na.rm = TRUE),  # hold at mean
  SM_mean_FC_long = mean(SM_mean_FC_long, na.rm = TRUE),  # hold at mean
  DA_mean_FC_long = mean(DA_mean_FC_long, na.rm = TRUE),  # hold at mean
  VA_mean_FC_long = mean(VA_mean_FC_long, na.rm = TRUE),  # hold at mean
  L_mean_FC_long = mean(L_mean_FC_long, na.rm = TRUE),  # hold at mean
  FP_mean_FC_long = mean(FP_mean_FC_long, na.rm = TRUE),  # hold at mean
  DMN_mean_FC_long = DMN_grid,   # these are the predictios that we're plotting
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_DMN <- fitted_values(
  bam_SA_exp_FCs_all_net, 
  data = newdata_DMN, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_DMN <- ggplot(
  smooth_preds_DMN, 
  aes(x = DMN_mean_FC_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "Mean strength of functional connections\nin the default-mode network", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_FC_DMN_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_DMN, width = 6, height = 5, dpi = 300)


