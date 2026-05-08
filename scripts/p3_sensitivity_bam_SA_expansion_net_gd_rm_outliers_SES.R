
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Development

# Supplementary analyses: adding SES (parental education) as covariate in analyses:
# S-A axis expansion analyses 
# (S-A expansion metric computed at the NETWORK-level) 
# (cf MAIN ANALYSIS) - with OUTLIERS REMOVED 

# Geodesic Distances analyses -> "mean length of functional connections"?
# (GDs = mean GD of template Schaefer 400 GD matrix, but individual top 10% functional connections)


# Testing for the effects of age/pds on S-A axis expansion (leg1) and of S-A axis expansion on rcoef/grad_flip (leg2) 
# Testing for the effects of age/pds on mean GD per network (leg1)
# Testing for the effects of mean GD per network on S-A axis expansion (network-level) <- trying to show S<-0->A

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
resdir = '/data/p_02667/development/results/supplementary/SES/'
resdir_fig = '/data/p_02667/development/results/supplementary/SES/figures/'


#### set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### S-A axis expansion metrics (computed at the network-level)

## Baseline 
mat_SA_expansion_metrics_baseline <- readMat(paste(datadir_local, 'SA_expansion_metrics_baseline_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_baseline)  # print contents of matfile

SA_expansion_net_baseline <- mat_SA_expansion_metrics_baseline$SA.expansion.net.baseline  # shape: 1, 3950
SA_expansion_net_baseline <- t(SA_expansion_net_baseline)  #transposing to get the shape: 3950, 1

# sub ID 
sub_ID_baseline <- mat_SA_expansion_metrics_baseline$sub.ID.baseline   # shape: 3950, 1

dim(SA_expansion_net_baseline)
dim(sub_ID_baseline)


## 2y follow-up
mat_SA_expansion_metrics_fu2y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu2y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu2y)  # print contents of matfile

SA_expansion_net_fu2y <- mat_SA_expansion_metrics_fu2y$SA.expansion.net.fu2y  # shape: 1, 1253
SA_expansion_net_fu2y <- t(SA_expansion_net_fu2y)  #transposing to get the shape: 1253, 1

# sub ID 
sub_ID_fu2y <- mat_SA_expansion_metrics_fu2y$sub.ID.fu2y   # shape: 1253, 1

dim(SA_expansion_net_fu2y)
dim(sub_ID_fu2y)


## 4y follow-up
mat_SA_expansion_metrics_fu4y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu4y_rm_outliers.mat', sep = ''))
names(mat_SA_expansion_metrics_fu4y)  # print contents of matfile

SA_expansion_net_fu4y <- mat_SA_expansion_metrics_fu4y$SA.expansion.net.fu4y  # shape: 1, 906
SA_expansion_net_fu4y <- t(SA_expansion_net_fu4y)  #transposing to get the shape: 906, 1

# sub ID 
sub_ID_fu4y <- mat_SA_expansion_metrics_fu4y$sub.ID.fu4y   # shape: 906, 1

dim(SA_expansion_net_fu4y)
dim(sub_ID_fu4y)



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




### Brain organization metrics

## baseline
# r coefs
list_SA_axis_corr_to_ref_baseline <- mat_SA_expansion_metrics_baseline$list.SA.axis.corr.to.ref.baseline  # shape: 1, 3950
list_SA_axis_corr_to_ref_baseline <- t(list_SA_axis_corr_to_ref_baseline)  #transposing to get the shape: 3950, 1

# gradient number corresponding to S-A axis
list_SA_axis_grad_num_baseline <- mat_SA_expansion_metrics_baseline$list.SA.axis.grad.num.baseline  # shape: 1, 3950
list_SA_axis_grad_num_baseline <- t(list_SA_axis_grad_num_baseline)  #transposing to get the shape: 3950, 1


## 2y follow-up
# r coefs
list_SA_axis_corr_to_ref_fu2y <- mat_SA_expansion_metrics_fu2y$list.SA.axis.corr.to.ref.fu2y  # shape: 1, 1253
list_SA_axis_corr_to_ref_fu2y <- t(list_SA_axis_corr_to_ref_fu2y)  #transposing to get the shape: 1253, 1

# gradient number corresponding to S-A axis
list_SA_axis_grad_num_fu2y <- mat_SA_expansion_metrics_fu2y$list.SA.axis.grad.num.fu2y  # shape: 1, 1253
list_SA_axis_grad_num_fu2y <- t(list_SA_axis_grad_num_fu2y)  #transposing to get the shape: 1253, 1


## 4y follow-up
# r coefs
list_SA_axis_corr_to_ref_fu4y <- mat_SA_expansion_metrics_fu4y$list.SA.axis.corr.to.ref.fu4y  # shape: 1, 906
list_SA_axis_corr_to_ref_fu4y <- t(list_SA_axis_corr_to_ref_fu4y)  #transposing to get the shape: 906, 1

# gradient number corresponding to S-A axis
list_SA_axis_grad_num_fu4y <- mat_SA_expansion_metrics_fu4y$list.SA.axis.grad.num.fu4y  # shape: 1, 906
list_SA_axis_grad_num_fu4y <- t(list_SA_axis_grad_num_fu4y)  #transposing to get the shape: 906, 1




### Mean Geodesic Distances 

## Baseline 

# V
V_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$V.mean.gd.baseline  # shape: 1, 3950
V_mean_gd_baseline <- t(V_mean_gd_baseline)  #transposing to get the shape: 3950, 1

# SM
SM_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$SM.mean.gd.baseline  # shape: 1, 3950
SM_mean_gd_baseline <- t(SM_mean_gd_baseline)  #transposing to get the shape: 3950, 1

# DA
DA_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$DA.mean.gd.baseline  # shape: 1, 3950
DA_mean_gd_baseline <- t(DA_mean_gd_baseline)  #transposing to get the shape: 3950, 1

# VA
VA_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$VA.mean.gd.baseline  # shape: 1, 3950
VA_mean_gd_baseline <- t(VA_mean_gd_baseline)  #transposing to get the shape: 3950, 1

# L
L_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$L.mean.gd.baseline  # shape: 1, 3950
L_mean_gd_baseline <- t(L_mean_gd_baseline)  #transposing to get the shape: 3950, 1

# FP
FP_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$FP.mean.gd.baseline  # shape: 1, 3950
FP_mean_gd_baseline <- t(FP_mean_gd_baseline)  #transposing to get the shape: 3950, 1

# DMN
DMN_mean_gd_baseline <- mat_SA_expansion_metrics_baseline$DMN.mean.gd.baseline  # shape: 1, 3950
DMN_mean_gd_baseline <- t(DMN_mean_gd_baseline)  #transposing to get the shape: 3950, 1



## Fu2y 

# V
V_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$V.mean.gd.fu2y  # shape: 1, 1253
V_mean_gd_fu2y <- t(V_mean_gd_fu2y)  #transposing to get the shape: 1253, 1

# SM
SM_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$SM.mean.gd.fu2y  # shape: 1, 1253
SM_mean_gd_fu2y <- t(SM_mean_gd_fu2y)  #transposing to get the shape: 1253, 1

# DA
DA_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$DA.mean.gd.fu2y  # shape: 1, 1253
DA_mean_gd_fu2y <- t(DA_mean_gd_fu2y)  #transposing to get the shape: 1253, 1

# VA
VA_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$VA.mean.gd.fu2y  # shape: 1, 1253
VA_mean_gd_fu2y <- t(VA_mean_gd_fu2y)  #transposing to get the shape: 1253, 1

# L
L_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$L.mean.gd.fu2y  # shape: 1, 1253
L_mean_gd_fu2y <- t(L_mean_gd_fu2y)  #transposing to get the shape: 1253, 1

# FP
FP_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$FP.mean.gd.fu2y  # shape: 1, 1253
FP_mean_gd_fu2y <- t(FP_mean_gd_fu2y)  #transposing to get the shape: 1253, 1

# DMN
DMN_mean_gd_fu2y <- mat_SA_expansion_metrics_fu2y$DMN.mean.gd.fu2y  # shape: 1, 1253
DMN_mean_gd_fu2y <- t(DMN_mean_gd_fu2y)  #transposing to get the shape: 1253, 1



## Fu4y 

# V
V_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$V.mean.gd.fu4y  # shape: 1, 906
V_mean_gd_fu4y <- t(V_mean_gd_fu4y)  #transposing to get the shape: 906, 1

# SM
SM_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$SM.mean.gd.fu4y  # shape: 1, 906
SM_mean_gd_fu4y <- t(SM_mean_gd_fu4y)  #transposing to get the shape: 906, 1

# DA
DA_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$DA.mean.gd.fu4y  # shape: 1, 906
DA_mean_gd_fu4y <- t(DA_mean_gd_fu4y)  #transposing to get the shape: 906, 1

# VA
VA_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$VA.mean.gd.fu4y  # shape: 1, 906
VA_mean_gd_fu4y <- t(VA_mean_gd_fu4y)  #transposing to get the shape: 906, 1

# L
L_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$L.mean.gd.fu4y  # shape: 1, 906
L_mean_gd_fu4y <- t(L_mean_gd_fu4y)  #transposing to get the shape: 906, 1

# FP
FP_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$FP.mean.gd.fu4y  # shape: 1, 906
FP_mean_gd_fu4y <- t(FP_mean_gd_fu4y)  #transposing to get the shape: 906, 1

# DMN
DMN_mean_gd_fu4y <- mat_SA_expansion_metrics_fu4y$DMN.mean.gd.fu4y  # shape: 1, 906
DMN_mean_gd_fu4y <- t(DMN_mean_gd_fu4y)  #transposing to get the shape: 906, 1





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


### Concatenate lists for longitudinal analyses

## S-A axis expansion
SA_expansion_net_long = c(SA_expansion_net_baseline, SA_expansion_net_fu2y, SA_expansion_net_fu4y)

str(SA_expansion_net_long)


## demo 
df_list <- list(abcd_demo_baseline, abcd_demo_fu2y, abcd_demo_fu4y)
abcd_demo_long <- bind_rows(df_list)

str(abcd_demo_long)


## r coefs
list_SA_axis_corr_to_ref_long = c(list_SA_axis_corr_to_ref_baseline, list_SA_axis_corr_to_ref_fu2y, list_SA_axis_corr_to_ref_fu4y)

str(list_SA_axis_corr_to_ref_long)


## S-A axis gradient numbers (for grad_flip)
list_SA_axis_grad_num_long = c(list_SA_axis_grad_num_baseline, list_SA_axis_grad_num_fu2y, list_SA_axis_grad_num_fu4y)

str(list_SA_axis_grad_num_long)


## Geodesic Distances
# by yeo network
V_mean_gd_long = c(V_mean_gd_baseline, V_mean_gd_fu2y, V_mean_gd_fu4y)
SM_mean_gd_long = c(SM_mean_gd_baseline, SM_mean_gd_fu2y, SM_mean_gd_fu4y)
DA_mean_gd_long = c(DA_mean_gd_baseline, DA_mean_gd_fu2y, DA_mean_gd_fu4y)
VA_mean_gd_long = c(VA_mean_gd_baseline, VA_mean_gd_fu2y, VA_mean_gd_fu4y)
L_mean_gd_long = c(L_mean_gd_baseline, L_mean_gd_fu2y, L_mean_gd_fu4y)
FP_mean_gd_long = c(FP_mean_gd_baseline, FP_mean_gd_fu2y, FP_mean_gd_fu4y)
DMN_mean_gd_long = c(DMN_mean_gd_baseline, DMN_mean_gd_fu2y, DMN_mean_gd_fu4y)



##### Transforming S-A axis expansion metric due to skew? -> currently haven't done it

# Plot histogram
hist(SA_expansion_net_long,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Distribution of S-A axis expansion values",
     xlab = "DMN_disp",
     ylab = "Frequency")




##### Fisher r-to-z transforming skewed data rcoef

list_SA_axis_corr_to_ref_long_z <- 0.5 * log((1 + list_SA_axis_corr_to_ref_long) / (1 - list_SA_axis_corr_to_ref_long))


##### Making the gradient flip variable binary (has the flip occurred or not) 
# 1: S-A axis is gradient 1; flip has occurred
# 0: S-A axis not gradient 1 (mostly gradient 2 but goes up to 6); flip has not yet occurred

grad_flip <- ifelse(list_SA_axis_grad_num_long == 1, 1, 0)



##### Mean GDs as short medium long

## by network

# confirmed mean GDs by network
mean(V_mean_gd_long)  # short
mean(SM_mean_gd_long)  # short
mean(DA_mean_gd_long)  # medium
mean(VA_mean_gd_long)  # medium
mean(L_mean_gd_long)  # long
mean(FP_mean_gd_long)  # long
mean(DMN_mean_gd_long)  # long

mean_gd_short_net <- (V_mean_gd_long + SM_mean_gd_long) / 2
mean_gd_medium_net <- (DA_mean_gd_long + VA_mean_gd_long) / 2
mean_gd_long_net <- (L_mean_gd_long + FP_mean_gd_long + DMN_mean_gd_long) / 3

mean(mean_gd_short_net)
mean(mean_gd_medium_net)
mean(mean_gd_long_net)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS -> S-A axis expansion analyses (legs 1 and 2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare the data
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Prepare the covariates

# select a subset of variables from the demo dataframe
covar_df <- abcd_demo_long %>%
  
  dplyr::select(src_subject_id_fmt, age_months, pds_p_score, sex, tot_SA, edu_parent, family_id, site_id, eventname) %>% 
  dplyr::rename("subject_id" = src_subject_id_fmt, "age" = age_months, "pds" = pds_p_score)


# set data type 
covar_df$age <- as.numeric(covar_df$age)
covar_df$pds <- as.numeric(covar_df$pds)
covar_df$sex <- as.factor(covar_df$sex)
covar_df$tot_SA <- as.numeric(covar_df$tot_SA)
covar_df$edu_parent <- as.numeric(covar_df$edu_parent)
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
# RUN BAM analyses for SA_expansion_net_long
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Leg 1

# S-A axis expansion ~ age + pds + sex + totSA + (1 | sub+fam+site)
# + edu parent

bam_SA_exp_leg1 <- bam(
  SA_expansion_net_long ~
    sex +
    edu_parent +
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

summary_bam_SA_exp_leg1_rm_outliers = summary(bam_SA_exp_leg1)

# Save summary
capture.output(summary_bam_SA_exp_leg1_rm_outliers, file = file.path(resdir, "summary_bam_SA_exp_leg1_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_leg1)



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### Plot

### Smooth Estimates

df_age <- smooth_estimates(bam_SA_exp_leg1, select = "s(age)")
df_pds <- smooth_estimates(bam_SA_exp_leg1, select = "s(pds)")


# To harmonize y-axes in panel figure: 
# Find the global range for smooth estimates (including confidence intervals)
smes_ylim <- range(
  c(df_age$.estimate - 2 * df_age$.se, df_age$.estimate + 2 * df_age$.se,
    df_pds$.estimate - 2 * df_pds$.se, df_pds$.estimate + 2 * df_pds$.se)
)

plot_smes_age = ggplot(df_age, aes(x = age, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "honeydew3", alpha = 0.25) +
  geom_line(color = "honeydew4", linewidth = 1.2) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on S-A axis expansion",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )

plot_smes_pds = ggplot(df_pds, aes(x = pds, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "thistle", alpha = 0.25) +
  geom_line(color = "thistle4", linewidth = 1.2) +
  labs(
    x = "Pubertal stage (PDS score)",
    y = "Partial effect on S-A axis expansion",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )

ggsave(paste(resdir_fig, 'bam_smes_leg1_age_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_leg1_pds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)


# To harmonize y-axes: Apply ylim to plots
plot_smes_age <- plot_smes_age + coord_cartesian(ylim = smes_ylim)
plot_smes_pds <- plot_smes_pds + coord_cartesian(ylim = smes_ylim)

# Combine plots in one figure
plot_smes_panel <- (plot_smes_age + patchwork::plot_spacer() + plot_smes_pds) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1)) & # Adjust 0.1 to increase/decrease the gap
  theme(
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(paste(resdir_fig, 'bam_smes_panel_leg1_ageEpds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_smes_panel, width = 11, height = 5, dpi = 300)



### Predicted marginal trajectories (using fitted_values() from gratia package)

## Age

# 1. Build prediction grid

age_grid<- seq(
  min(covar_df$age, na.rm=TRUE),
  max(covar_df$age, na.rm=TRUE),
  length.out=200
)

newdata_age <- expand.grid(
  age = age_grid,
  sex = levels(covar_df$sex),
  edu_parent = mean(covar_df$edu_parent, na.rm = TRUE),  # hold edu_parent at mean
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_age <- fitted_values(
  bam_SA_exp_leg1, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
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
  
  labs(x = "Age (in months)", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_leg1_age_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


## PDS

# 1. Build prediction grid for PDS
pds_grid <- seq(
  min(covar_df$pds, na.rm = TRUE),
  max(covar_df$pds, na.rm = TRUE),
  length.out = 200
)

newdata_pds <- expand.grid(
  pds = pds_grid,
  sex = levels(covar_df$sex),
  edu_parent = mean(covar_df$edu_parent, na.rm = TRUE),  # hold edu_parent at mean
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold totSA at mean
  site_id = covar_df$site_id[1],
  site.family = covar_df$site.family[1],
  site.family.subject = covar_df$site.family.subject[1]
)

# 2. Generate predictions using gratia
smooth_preds_pds <- fitted_values(
  bam_SA_exp_leg1, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)

# 3. Plot marginal trajectories for PDS
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
  
  labs(x = "Pubertal stage (PDS score)", y = "S-A axis expansion", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_leg1_pds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)



# To harmonize y-axes in panel figure: 
# Find the global range for marginal trajectories
marg_ylim <- range(
  c(smooth_preds_age$.lower_ci, smooth_preds_age$.upper_ci,
    smooth_preds_pds$.lower_ci, smooth_preds_pds$.upper_ci)
)

# Update Age Plot
plot_margtraj_age <- plot_margtraj_age + 
  scale_y_continuous(
    limits = marg_ylim, # Set the shared limits here
  )

# Update PDS Plot
plot_margtraj_pds <- plot_margtraj_pds + 
  scale_y_continuous(
    limits = marg_ylim, # Use the exact same limits
  )

# Combine plots
plot_margtraj_panel <- (plot_margtraj_age + patchwork::plot_spacer() + plot_margtraj_pds) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1), guides = "collect") &   # Adjust 0.1 to increase/decrease the gap
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_panel_leg1_ageEpds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_panel, width = 11, height = 5, dpi = 300)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



### Leg 2

# rcoef ~ S-A axis expansion + age + pds + sex + totSA + (1 | sub+fam+site)
# + edu parent

bam_SA_exp_leg2_rcoef <- bam(
  list_SA_axis_corr_to_ref_long_z ~
    sex +
    edu_parent +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    s(SA_expansion_net_long, k = 10, fx = FALSE) +  # default
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

summary_bam_SA_exp_leg2_rcoef_rm_outliers = summary(bam_SA_exp_leg2_rcoef)

# Save summary
capture.output(summary_bam_SA_exp_leg2_rcoef_rm_outliers, file = file.path(resdir, "summary_bam_SA_exp_leg2_rcoef_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_leg2_rcoef)


##### Plot

### Smooth Estimates

df_SAexp <- smooth_estimates(bam_SA_exp_leg2_rcoef, select = "s(SA_expansion_net_long)")

plot_smes_SAexp = ggplot(df_SAexp, aes(x = SA_expansion_net_long, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "gold", alpha = 0.25) +
  geom_line(color = "gold2", linewidth = 1.2) +
  labs(
    x = "S-A axis expansion",
    y = "Partial effect on similarity to adult S-A axis\n(Δ Fisher z)",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )

ggsave(paste(resdir_fig, 'bam_smes_leg2_SAexp_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_smes_SAexp, width = 5, height = 5, dpi = 300)



### Predicted marginal trajectories (using fitted_values() from gratia package)

# 1. Build prediction grid
SAexp_grid<- seq(
  min(SA_expansion_net_long, na.rm=TRUE),
  max(SA_expansion_net_long, na.rm=TRUE),
  length.out=200
)

newdata_SAexp <- expand.grid(
  SA_expansion_net_long = SAexp_grid,
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  edu_parent = mean(covar_df$edu_parent, na.rm = TRUE),  # hold edu_parent at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)

# 2. Generate the predictions
smooth_preds_SAexp <- fitted_values(
  bam_SA_exp_leg2_rcoef, 
  data = newdata_SAexp, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)

# 3. Plot marginal trajectories
plot_margtraj_SAexp <- ggplot(
  smooth_preds_SAexp, 
  aes(x = SA_expansion_net_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  # Dual Axis Configuration
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., # No transformation of the scale itself (stay in alignment)
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), # Transform the labels only (and format as string with fitted paddings)
      breaks = derive() # This ensures it uses the exact same tick positions as the Z axis
    )) +
  
  labs(x = "S-A axis expansion", y = "Similarity to adult S-A axis (Fisher z)", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_leg2_SAexp_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_margtraj_SAexp, width = 6, height = 5, dpi = 300)




##### grad_flip ~ S-A axis expansion + pds + sex + totSA + (1 | sub+fam+site)
# + edu parent

bam_SA_exp_leg2_gradflip <- bam(
  grad_flip ~ 
    sex +
    edu_parent +
    s(age, k = 10, fx = FALSE) +  # default (including age because -although it wasnt a significant predictor in gradflip main analysis- it was for M)
    s(pds, k = 10, fx = FALSE) +  # default
    s(SA_expansion_net_long, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df,
  family = binomial(link = "logit"),
  discrete = TRUE,        # big speed-up for large datasets
  method = "fREML",        # recommended for models with random effects
  nthreads = 36,            # number of cores for parallel computation -> can also try: parallel::detectCores()
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_SA_exp_leg2_gradflip_rm_outliers = summary(bam_SA_exp_leg2_gradflip)

# Save summary
capture.output(summary_bam_SA_exp_leg2_gradflip_rm_outliers, file = file.path(resdir, "summary_bam_SA_exp_leg2_gradflip_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_leg2_gradflip)



##### Plot

### Smooth Estimates

df_SAexp <- smooth_estimates(bam_SA_exp_leg2_gradflip, select = "s(SA_expansion_net_long)")

plot_smes_SAexp = ggplot(df_SAexp, aes(x = SA_expansion_net_long, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "gold", alpha = 0.25) +
  geom_line(color = "gold2", linewidth = 1.2) +
  labs(
    x = "S-A axis expansion",
    y = "Partial effect on gradient flip (log-odds)",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )

ggsave(paste(resdir_fig, 'bam_smes_leg2_SAexp_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_smes_SAexp, width = 5, height = 5, dpi = 300)



### Predicted marginal trajectories (using fitted_values() from gratia package)

## S-A axis expansion

# 1. Build prediction grid
SAexp_grid<- seq(
  min(SA_expansion_net_long, na.rm=TRUE),
  max(SA_expansion_net_long, na.rm=TRUE),
  length.out=200
)

newdata_SAexp <- expand.grid(
  SA_expansion_net_long = SAexp_grid,
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  edu_parent = mean(covar_df$edu_parent, na.rm = TRUE),  # hold edu_parent at mean
  sex = levels(covar_df$sex),
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)

# 2. Generate the predictions
smooth_preds_SAexp <- fitted_values(
  bam_SA_exp_leg2_gradflip, 
  data = newdata_SAexp, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)

# 3. Plot marginal trajectories
plot_margtraj_SAexp <- ggplot(
  smooth_preds_SAexp, 
  aes(x = SA_expansion_net_long, y = .fitted, color = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci),  # Use .lower and .upper (instead of manual +/- 2*se) in order to get the limits of the credible interval on the fitted values, on the specified scale
    alpha = 0.25, 
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  scale_fill_manual(values = c("F" = "indianred3", "M" = "lightblue3")) + #, labels = c("F" = "Female", "M" = "Male")) +
  
  labs(x = "S-A axis expansion", y = "Probability of gradient flip", color = "Sex", fill = "Sex") +
  #scale_y_continuous(labels = scales::label_percent()) + # Optional: converts 0.5 to 50%
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_leg2_SAexp_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_margtraj_SAexp, width = 6, height = 5, dpi = 300)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS -> age and pds effects on mean GD per network

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM analyses: Mean GD (net) ~ age + pds + sex + totSA + (1 | sub+fam+site)
# + edu parent
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### V

bam_GD_V <- bam(
  V_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_V = summary(bam_GD_V)

# Save summary
capture.output(summary_bam_GD_V, file = file.path(resdir, "summary_bam_GD_V_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_V)

# Extract the data for plotting
df_bam_GD_V_age <- smooth_estimates(bam_GD_V, select = "s(age)")
df_bam_GD_V_pds <- smooth_estimates(bam_GD_V, select = "s(pds)")



### SM

bam_GD_SM <- bam(
  SM_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_SM = summary(bam_GD_SM)

# Save summary
capture.output(summary_bam_GD_SM, file = file.path(resdir, "summary_bam_GD_SM_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_SM)

# Extract the data for plotting
df_bam_GD_SM_age <- smooth_estimates(bam_GD_SM, select = "s(age)")
df_bam_GD_SM_pds <- smooth_estimates(bam_GD_SM, select = "s(pds)")


### DA

bam_GD_DA <- bam(
  DA_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_DA = summary(bam_GD_DA)

# Save summary
capture.output(summary_bam_GD_DA, file = file.path(resdir, "summary_bam_GD_DA_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_DA)

# Extract the data for plotting
df_bam_GD_DA_age <- smooth_estimates(bam_GD_DA, select = "s(age)")
df_bam_GD_DA_pds <- smooth_estimates(bam_GD_DA, select = "s(pds)")


### VA

bam_GD_VA <- bam(
  VA_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_VA = summary(bam_GD_VA)

# Save summary
capture.output(summary_bam_GD_VA, file = file.path(resdir, "summary_bam_GD_VA_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_VA)

# Extract the data for plotting
df_bam_GD_VA_age <- smooth_estimates(bam_GD_VA, select = "s(age)")
df_bam_GD_VA_pds <- smooth_estimates(bam_GD_VA, select = "s(pds)")


### L

bam_GD_L <- bam(
  L_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_L = summary(bam_GD_L)

# Save summary
capture.output(summary_bam_GD_L, file = file.path(resdir, "summary_bam_GD_L_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_L)


# Extract the data for plotting
df_bam_GD_L_age <- smooth_estimates(bam_GD_L, select = "s(age)")
df_bam_GD_L_pds <- smooth_estimates(bam_GD_L, select = "s(pds)")


### FP

bam_GD_FP <- bam(
  FP_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_FP = summary(bam_GD_FP)

# Save summary
capture.output(summary_bam_GD_FP, file = file.path(resdir, "summary_bam_GD_FP_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_FP)

# Extract the data for plotting
df_bam_GD_FP_age <- smooth_estimates(bam_GD_FP, select = "s(age)")
df_bam_GD_FP_pds <- smooth_estimates(bam_GD_FP, select = "s(pds)")


### DMN

bam_GD_DMN <- bam(
  DMN_mean_gd_long ~
    sex +
    edu_parent +
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

summary_bam_GD_DMN = summary(bam_GD_DMN)

# Save summary
capture.output(summary_bam_GD_DMN, file = file.path(resdir, "summary_bam_GD_DMN_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_GD_DMN)

# Extract the data for plotting
df_bam_GD_DMN_age <- smooth_estimates(bam_GD_DMN, select = "s(age)")
df_bam_GD_DMN_pds <- smooth_estimates(bam_GD_DMN, select = "s(pds)")



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of age and pds effects on Mean GD per network -> all in one
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
  c(df_bam_GD_V_age$.estimate - 2 * df_bam_GD_V_age$.se, df_bam_GD_V_age$.estimate + 2 * df_bam_GD_V_age$.se,
    df_bam_GD_V_pds$.estimate - 2 * df_bam_GD_V_pds$.se, df_bam_GD_V_pds$.estimate + 2 * df_bam_GD_V_pds$.se,
    df_bam_GD_SM_age$.estimate - 2 * df_bam_GD_SM_age$.se, df_bam_GD_SM_age$.estimate + 2 * df_bam_GD_SM_age$.se,
    df_bam_GD_SM_pds$.estimate - 2 * df_bam_GD_SM_pds$.se, df_bam_GD_SM_pds$.estimate + 2 * df_bam_GD_SM_pds$.se,
    df_bam_GD_DA_age$.estimate - 2 * df_bam_GD_DA_age$.se, df_bam_GD_DA_age$.estimate + 2 * df_bam_GD_DA_age$.se,
    df_bam_GD_DA_pds$.estimate - 2 * df_bam_GD_DA_pds$.se, df_bam_GD_DA_pds$.estimate + 2 * df_bam_GD_DA_pds$.se,
    df_bam_GD_VA_age$.estimate - 2 * df_bam_GD_VA_age$.se, df_bam_GD_VA_age$.estimate + 2 * df_bam_GD_VA_age$.se,
    df_bam_GD_VA_pds$.estimate - 2 * df_bam_GD_VA_pds$.se, df_bam_GD_VA_pds$.estimate + 2 * df_bam_GD_VA_pds$.se,
    df_bam_GD_L_age$.estimate - 2 * df_bam_GD_L_age$.se, df_bam_GD_L_age$.estimate + 2 * df_bam_GD_L_age$.se,
    df_bam_GD_L_pds$.estimate - 2 * df_bam_GD_L_pds$.se, df_bam_GD_L_pds$.estimate + 2 * df_bam_GD_L_pds$.se,
    df_bam_GD_FP_age$.estimate - 2 * df_bam_GD_FP_age$.se, df_bam_GD_FP_age$.estimate + 2 * df_bam_GD_FP_age$.se,
    df_bam_GD_FP_pds$.estimate - 2 * df_bam_GD_FP_pds$.se, df_bam_GD_FP_pds$.estimate + 2 * df_bam_GD_FP_pds$.se,
    df_bam_GD_DMN_age$.estimate - 2 * df_bam_GD_DMN_age$.se, df_bam_GD_DMN_age$.estimate + 2 * df_bam_GD_DMN_age$.se,
    df_bam_GD_DMN_pds$.estimate - 2 * df_bam_GD_DMN_pds$.se, df_bam_GD_DMN_pds$.estimate + 2 * df_bam_GD_DMN_pds$.se)
)


### Age

# Combine and set factor levels
df_age_all <- bind_rows(
  df_bam_GD_V_age  %>% mutate(Network = "V"),
  df_bam_GD_SM_age %>% mutate(Network = "SM"),
  df_bam_GD_DA_age %>% mutate(Network = "DA"),
  df_bam_GD_VA_age %>% mutate(Network = "VA"),
  df_bam_GD_L_age  %>% mutate(Network = "L"),
  df_bam_GD_FP_age %>% mutate(Network = "FP"),
  df_bam_GD_DMN_age %>% mutate(Network = "DMN")
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
    y = "Partial effect on mean length of\nfunctional connections",
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
  df_bam_GD_V_pds  %>% mutate(Network = "V"),
  df_bam_GD_SM_pds %>% mutate(Network = "SM"),
  df_bam_GD_DA_pds %>% mutate(Network = "DA"),
  df_bam_GD_VA_pds %>% mutate(Network = "VA"),
  df_bam_GD_L_pds  %>% mutate(Network = "L"),
  df_bam_GD_FP_pds %>% mutate(Network = "FP"),
  df_bam_GD_DMN_pds %>% mutate(Network = "DMN")
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
    y = "Partial effect on mean length of\nfunctional connections",
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
#ggsave(paste(resdir_fig, 'bam_smes_age_on_GD_per_network_all_rm_outliers.svg', sep = ''), plot = plot_smes_age_all, width = 6, height = 5, dpi = 300)
#ggsave(paste(resdir_fig, 'bam_smes_pds_on_GD_per_network_all_rm_outliers.svg', sep = ''), plot = plot_smes_pds_all, width = 6, height = 5, dpi = 300)


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
ggsave(paste(resdir_fig, 'bam_smes_panel_ageEpds_on_GD_per_network_all_rm_outliers.svg', sep = ''), plot = plot_smes_panel, width = 11, height = 5, dpi = 300)




### only displaying significant results (manually selecting them)

### Age

# Define your desired order
network_order_sig_age <- c(
  "V", 
  "SM", 
  #"DA", 
  "VA" 
  #"L", 
  #"FP", 
  #"DMN"
)


# Combine and set factor levels
df_age_all_sig <- bind_rows(
  df_bam_GD_V_age  %>% mutate(Network = "V"),
  df_bam_GD_SM_age %>% mutate(Network = "SM"),
  #df_bam_GD_DA_age %>% mutate(Network = "DA"),
  df_bam_GD_VA_age %>% mutate(Network = "VA"),
  #df_bam_GD_L_age  %>% mutate(Network = "L"),
  #df_bam_GD_FP_age %>% mutate(Network = "FP"),
  #df_bam_GD_DMN_age %>% mutate(Network = "DMN")
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
    y = "Partial effect on mean length of\nfunctional connections",
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
  "SM", 
  #"DA", 
  "VA", 
  "L", 
  "FP"
  #,
  #"DMN"
)

# Combine all smooth estimate data frames into one
df_pds_all_sig <- bind_rows(
  #df_bam_GD_V_pds  %>% mutate(Network = "V"),
  df_bam_GD_SM_pds %>% mutate(Network = "SM"),
  #df_bam_GD_DA_pds %>% mutate(Network = "DA"),
  df_bam_GD_VA_pds %>% mutate(Network = "VA"),
  df_bam_GD_L_pds  %>% mutate(Network = "L"),
  df_bam_GD_FP_pds %>% mutate(Network = "FP")
  #,
  #df_bam_GD_DMN_pds %>% mutate(Network = "DMN")
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
    y = "Partial effect on mean length of\nfunctional connections",
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
#ggsave(paste(resdir_fig, 'bam_smes_age_on_GD_per_network_all_sig_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_age_all_sig, width = 6, height = 5, dpi = 300)
#ggsave(paste(resdir_fig, 'bam_smes_pds_on_GD_per_network_all_sig_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_pds_all_sig, width = 6, height = 5, dpi = 300)


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
ggsave(paste(resdir_fig, 'bam_smes_panel_ageEpds_on_GD_per_network_all_sig_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_panel_all_sig, width = 11, height = 5, dpi = 300)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS -> mean GDs effects on S-A axis expansion

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Network-level (Yeo 7)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM analysis: S-A axis expansion ~ age + pds + GDs(of each network) + sex + totSA + (1 | sub+fam+site)
# + edu parent
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bam_SA_exp_GDs_all_net <- bam(
  SA_expansion_net_long ~
    sex +
    edu_parent +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    s(V_mean_gd_long, k = 10, fx = FALSE) +  # default
    s(SM_mean_gd_long, k = 10, fx = FALSE) +  # default
    s(DA_mean_gd_long, k = 10, fx = FALSE) +  # default
    s(VA_mean_gd_long, k = 10, fx = FALSE) +  # default
    s(L_mean_gd_long, k = 10, fx = FALSE) +  # default
    s(FP_mean_gd_long, k = 10, fx = FALSE) +  # default
    s(DMN_mean_gd_long, k = 10, fx = FALSE) +  # default
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

summary_bam_SA_exp_GDs_all_net_rm_outliers = summary(bam_SA_exp_GDs_all_net)

# Save summary
capture.output(summary_bam_SA_exp_GDs_all_net_rm_outliers, file = file.path(resdir, "summary_bam_SA_exp_GDs_all_net_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_GDs_all_net)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of effects (on SA exp) of mean GDs per network -> all in one
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
df_V <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(V_mean_gd_long)")
df_SM <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(SM_mean_gd_long)")
df_DA <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(DA_mean_gd_long)")
df_VA <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(VA_mean_gd_long)")
df_L <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(L_mean_gd_long)")
df_FP <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(FP_mean_gd_long)")
df_DMN <- smooth_estimates(bam_SA_exp_GDs_all_net, select = "s(DMN_mean_gd_long)")


# Combine and set factor levels (by also renaming the predictor column to "mean_gd" for all networks before binding)

df_all_net <- bind_rows(
  df_V  %>% dplyr::rename(mean_gd = V_mean_gd_long)   %>% mutate(Network = "V"),
  df_SM %>% dplyr::rename(mean_gd = SM_mean_gd_long)  %>% mutate(Network = "SM"),
  df_DA %>% dplyr::rename(mean_gd = DA_mean_gd_long)  %>% mutate(Network = "DA"),
  df_VA %>% dplyr::rename(mean_gd = VA_mean_gd_long)  %>% mutate(Network = "VA"),
  df_L  %>% dplyr::rename(mean_gd = L_mean_gd_long)   %>% mutate(Network = "L"),
  df_FP %>% dplyr::rename(mean_gd = FP_mean_gd_long)  %>% mutate(Network = "FP"),
  df_DMN %>% dplyr::rename(mean_gd = DMN_mean_gd_long) %>% mutate(Network = "DMN")
) %>%
  mutate(Network = factor(Network, levels = network_order))  # set legend order here


# Plot (effects of GDs of all networks on S-A axis expansion)
plot_smes_all_net <- ggplot(df_all_net, aes(x = mean_gd, y = .estimate, color = Network, fill = Network)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = network_colors) +
  scale_fill_manual(values = network_colors) +
  labs(
    x = "Mean length of functional connections",
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


# before saving decide if you want with or without CIs (and remove or add "_CIs" accordingly)
ggsave(paste(resdir_fig, 'bam_smes_all_networks_GDs_on_SAexp_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_all_net, width = 6, height = 5, dpi = 300)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Summarizing mean short/medium/long range connects (based on Yeo network: short (V,SM), medium (DA,VA), long (L,FP, DMN))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


str(mean_gd_short_net)
str(mean_gd_medium_net)
str(mean_gd_long_net)



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM analysis: S-A axis expansion ~ age + pds + GDs(sml) + sex + totSA + (1 | sub+fam+site)
# + edu parent
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


bam_SA_exp_GDs_sml_net <- bam(
  SA_expansion_net_long ~
    sex +
    edu_parent +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    s(mean_gd_short_net, k = 10, fx = FALSE) +  # default
    s(mean_gd_medium_net, k = 10, fx = FALSE) +  # default
    s(mean_gd_long_net, k = 10, fx = FALSE) +  # default
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

summary_bam_SA_exp_GDs_sml_net = summary(bam_SA_exp_GDs_sml_net)

# Save summary
capture.output(summary_bam_SA_exp_GDs_sml_net, file = file.path(resdir, "summary_bam_SA_exp_GDs_sml_net_rm_outliers.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_SA_exp_GDs_sml_net)



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting of effects (on SA exp) of mean GDs per network range length -> all in one
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define your custom colors
range_length_colors <- c(
  "Short" = "lightgoldenrod1",
  "Medium" = "olivedrab3",
  "Long" = "darkolivegreen"
)

# Define your desired order
connections_range_length <- c("Short", "Medium", "Long")


# Extract the data for more control
df_short <- smooth_estimates(bam_SA_exp_GDs_sml_net, select = "s(mean_gd_short_net)")
df_medium <- smooth_estimates(bam_SA_exp_GDs_sml_net, select = "s(mean_gd_medium_net)")
df_long <- smooth_estimates(bam_SA_exp_GDs_sml_net, select = "s(mean_gd_long_net)")


# Combine and set factor levels (by also renaming the predictor column to "mean_gd" for all networks before binding)

df_all_net <- bind_rows(
  df_short  %>% dplyr::rename(mean_gd = mean_gd_short_net)   %>% mutate(Length = "Short"),
  df_medium %>% dplyr::rename(mean_gd = mean_gd_medium_net)  %>% mutate(Length = "Medium"),
  df_long %>% dplyr::rename(mean_gd = mean_gd_long_net)  %>% mutate(Length = "Long")
) %>%
  mutate(Length = factor(Length, levels = connections_range_length))  # set legend order here


# Plot (effects of GDs of all networks on S-A axis expansion)
plot_smes_all_sml_net <- ggplot(df_all_net, aes(x = mean_gd, y = .estimate, color = Length, fill = Length)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = range_length_colors) +
  scale_fill_manual(values = range_length_colors) +
  labs(
    x = "Mean length of functional connections",
    y = "Partial effect on S-A axis expansion",
    color = "Network\nlength",
    fill = "Network\nlength"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
    #legend.position = "right"
  )

plot_smes_all_sml_net


# before saving decide if you want with or without CIs (and remove or add "_CIs" accordingly)
ggsave(paste(resdir_fig, 'bam_smes_sml_net_GDs_on_SAexp_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_all_sml_net, width = 6, height = 5, dpi = 300)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# To harmonize y-axes in panel figures (from 2 models: both from all net and sml): 
# Find the global range for smooth estimates (including confidence intervals)
smes_ylim <- range(
  c(df_V$.estimate - 2 * df_V$.se, df_V$.estimate + 2 * df_V$.se,
    df_SM$.estimate - 2 * df_SM$.se, df_SM$.estimate + 2 * df_SM$.se,
    df_DA$.estimate - 2 * df_DA$.se, df_DA$.estimate + 2 * df_DA$.se,
    df_VA$.estimate - 2 * df_VA$.se, df_VA$.estimate + 2 * df_VA$.se,
    df_L$.estimate - 2 * df_L$.se, df_L$.estimate + 2 * df_L$.se,
    df_FP$.estimate - 2 * df_FP$.se, df_FP$.estimate + 2 * df_FP$.se,
    df_DMN$.estimate - 2 * df_DMN$.se, df_DMN$.estimate + 2 * df_DMN$.se,
    df_short$.estimate - 2 * df_short$.se, df_short$.estimate + 2 * df_short$.se,
    df_medium$.estimate - 2 * df_medium$.se, df_medium$.estimate + 2 * df_medium$.se,
    df_long$.estimate - 2 * df_long$.se, df_long$.estimate + 2 * df_long$.se)
)


plot_smes_all_net <- plot_smes_all_net + coord_cartesian(ylim = smes_ylim)
plot_smes_all_sml_net <- plot_smes_all_sml_net + coord_cartesian(ylim = smes_ylim)

# Combine plots in one figure
plot_smes_panel_allnetGDs_sml <- (plot_smes_all_net + patchwork::plot_spacer() + plot_smes_all_sml_net) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1)) & # Adjust 0.1 to increase/decrease the gap
  theme(
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# be careful at current CIs when plotting - only need to plot this with CIs
ggsave(paste(resdir_fig, 'bam_smes_panel_allnetGDsEsml_net_GDs_on_SAexp_rm_outliers_CIs.svg', sep = ''), plot = plot_smes_panel_allnetGDs_sml, width = 11, height = 5, dpi = 300)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

