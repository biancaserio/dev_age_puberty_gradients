# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Development

# Supplementary analyses: LMER to check relationship between brain developmental regression (negative gradient flip and decreasing similarity to adult S-A axis) and CBCL total symptom scores

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

require(lme4)  # for lmer (linear mixed effects model)
require(lmerTest)  # to obtain p-values for lmer - this actually overrides lme4's lmer() and prints the p-values for the fixed effects (which aren't present in lme4's lmer())

library(svglite)  # To save figures as vector files


### Clear environment
rm(list = ls())


#### set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/p_02667/development/results/supplementary/grad_flip_understanding/'
datadir_local = '/data/p_02667/development/results/supplementary/grad_flip_understanding/'
resdir = '/data/p_02667/development/results/supplementary/grad_flip_understanding/'
resdir_fig = '/data/p_02667/development/results/supplementary/grad_flip_understanding/'



#### set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


abcd_long_gradflip_rcoef = read.csv(paste(datadir_local, 'abcd_long_gradflip_rcoef.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

str(abcd_long_gradflip_rcoef)




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# RUNNING MODELS

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare the data
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Prepare the covariates

# compute log of CBCL tot scores to dataframe -> log1p(x) computes log(1 + x) <- given that skew in CBCL symptoms in ABCD (healthy sample)
# The “+1” is needed because many children have a CBCL score of 0, and log(0) is undefined (−∞)
# Adding 1 first maps 0 to 0, and the transform remains monotonic, so higher scores still mean more symptoms
abcd_long_gradflip_rcoef$cbcl_tot_log <- log1p(abcd_long_gradflip_rcoef$cbcl_tot)

hist(abcd_long_gradflip_rcoef$cbcl_tot)
hist(abcd_long_gradflip_rcoef$cbcl_tot_log)

# select a subset of variables from the dataframe
covar_df <- abcd_long_gradflip_rcoef %>%

  dplyr::select(src_subject_id_fmt, age_months, pds_p_score, sex, cbcl_tot_log, family_id,
              site_id, eventname, neg_flip, rcoef_dec, n_visits)  %>% 
  dplyr::rename("subject_id" = src_subject_id_fmt, "age" = age_months, "pds" = pds_p_score)


# set data type 
covar_df$age <- as.numeric(covar_df$age)
covar_df$pds <- as.numeric(covar_df$pds)
covar_df$sex <- as.factor(covar_df$sex)
covar_df$cbcl_tot_log <- as.numeric(covar_df$cbcl_tot_log)
covar_df$family_id <- as.factor(covar_df$family_id)
covar_df$site_id <- as.factor(covar_df$site_id)
covar_df$subject_id <- as.factor(covar_df$subject_id)
covar_df$eventname <- as.factor(covar_df$eventname)

str(covar_df)


### Determining the random effects variables (interactions) needed to code sub nested in fam nested in site

any(duplicated(paste(covar_df$site_id, covar_df$family_id)))  # TRUE
# TRUE  => family labels repeat across sites (not globally unique)
# FALSE => family labels are unique per site combination

any(duplicated(paste(covar_df$site_id, covar_df$family_id, covar_df$subject_id)))  # TRUE
# TRUE  => family and subject labels repeat across sites (not globally unique)
# FALSE => family and subject labels are unique per site combination

# fyi
length(unique(covar_df$site_id)) 
length(unique(interaction(covar_df$site_id, covar_df$family_id, drop=TRUE)))  
length(unique(interaction(covar_df$site_id, covar_df$family_id, covar_df$subject_id, drop=TRUE)))  


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
# RUN Mixed Effects Regression on r coefficients
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Fit a linear mixed effects model: DV ~ Sex + PDS score + age + tot_SA + random nested effect(site/family relatedness/subject id)

# Standardize continuous predictors
covar_df$age_z     <- scale(covar_df$age)
covar_df$pds_z     <- scale(covar_df$pds)

covar_df$age_z <- as.numeric(scale(covar_df$age))
covar_df$pds_z <- as.numeric(scale(covar_df$pds))
covar_df$cbcl_tot_log_z <- as.numeric(scale(covar_df$cbcl_tot_log))


# also scaling r_coef, but doing this directly in the model because naming is confusing (would be _z_z)
# note: Fisher r-to-z and scaling steps are not redundant, they serve different purposes:
# Fisher r-to-z → statistical normality
# scale() → comparability and standardized coefficients (“per 1 SD” units)

site_id = covar_df$site_id
family_id = covar_df$family_id
subject_id = covar_df$subject_id

m <- lmer(cbcl_tot_log_z ~ neg_flip * rcoef_dec + sex + pds_z + age_z +
            (1 | site_id / family_id / subject_id),
          data = covar_df, REML = FALSE)

summary_m = summary(m)

m2 <- lmer(cbcl_tot_log_z ~ neg_flip + rcoef_dec + sex + pds_z + age_z +
            (1 | site_id / family_id / subject_id),
          data = covar_df, REML = FALSE)

summary_m2 = summary(m2)


# Save summary
capture.output(summary_m, file = file.path(resdir, "summary_lmer_cbcl_regression_gradflip_rcoef.txt"))






