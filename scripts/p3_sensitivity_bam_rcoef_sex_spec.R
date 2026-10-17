# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Development

# Supplementary analyses: Age and PDS on similarity to adult S-A axis (rcoef - continuous variabile - guassian) - Sex specific models 

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
datadir_local = '/data/p_02667/development/data/pds_cleaned/'
resdir = '/data/p_02667/development/results/supplementary/sex_spec/'
resdir_fig = '/data/p_02667/development/results/supplementary/sex_spec/figures/'



#### set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Similarity to adult

## Baseline (loading this because it has the subject ID list)
mat_SA_axis_NOT_aligned_baseline <- readMat(paste(datadir_local, 'array_SA_axis_NOT_aligned_baseline.mat', sep = ''))


# r coefs
list_SA_axis_corr_to_ref_baseline <- mat_SA_axis_NOT_aligned_baseline$list.SA.axis.corr.to.ref.baseline  # shape: 1, 4064
list_SA_axis_corr_to_ref_baseline <- t(list_SA_axis_corr_to_ref_baseline)  #transposing to get the shape: 4064, 1

# sub ID 
sub_ID_baseline <- mat_SA_axis_NOT_aligned_baseline$sub.ID.baseline   # shape: 4064, 1


dim(list_SA_axis_corr_to_ref_baseline)
dim(sub_ID_baseline)


## 2y follow-up
mat_SA_axis_NOT_aligned_fu2y <- readMat(paste(datadir_local, 'array_SA_axis_NOT_aligned_fu2y.mat', sep = ''))

# r coefs
list_SA_axis_corr_to_ref_fu2y <- mat_SA_axis_NOT_aligned_fu2y$list.SA.axis.corr.to.ref.fu2y  # shape: 1, 1296
list_SA_axis_corr_to_ref_fu2y <- t(list_SA_axis_corr_to_ref_fu2y)  #transposing to get the shape: 1296, 1

# sub ID 
sub_ID_fu2y <- mat_SA_axis_NOT_aligned_fu2y$sub.ID.fu2y   # shape: 1296, 1


dim(list_SA_axis_corr_to_ref_fu2y)
dim(sub_ID_fu2y)



## 4y follow-up
mat_SA_axis_NOT_aligned_fu4y <- readMat(paste(datadir_local, 'array_SA_axis_NOT_aligned_fu4y.mat', sep = ''))

# r coefs
list_SA_axis_corr_to_ref_fu4y <- mat_SA_axis_NOT_aligned_fu4y$list.SA.axis.corr.to.ref.fu4y  # shape: 1, 963
list_SA_axis_corr_to_ref_fu4y <- t(list_SA_axis_corr_to_ref_fu4y)  #transposing to get the shape: 963, 1

# sub ID 
sub_ID_fu4y <- mat_SA_axis_NOT_aligned_fu4y$sub.ID.fu4y   # shape: 963, 1


dim(list_SA_axis_corr_to_ref_fu4y)
dim(sub_ID_fu4y)



### Demographics with covariates

# baseline
abcd_demo_baseline = read.csv(paste(datadir_local, 'abcd_demo_baseline_clean.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

str(abcd_demo_baseline)


# 2y
abcd_demo_fu2y = read.csv(paste(datadir_local, 'abcd_demo_fu2y_clean.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')


str(abcd_demo_fu2y)


# 4y
abcd_demo_fu4y = read.csv(paste(datadir_local, 'abcd_demo_fu4y_clean.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

str(abcd_demo_fu4y)




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


### Concatenate r coefs lists and demo dataframes for longitudinal analyses

# r coefs
list_SA_axis_corr_to_ref_long = c(list_SA_axis_corr_to_ref_baseline, list_SA_axis_corr_to_ref_fu2y, list_SA_axis_corr_to_ref_fu4y)

str(list_SA_axis_corr_to_ref_long)


# demo 
df_list <- list(abcd_demo_baseline, abcd_demo_fu2y, abcd_demo_fu4y)
abcd_demo_long <- bind_rows(df_list)

str(abcd_demo_long)




##### Fisher r-to-z transforming skewed data 

### Distribution of r coefficients -> right skew

# Plot histogram
hist(list_SA_axis_corr_to_ref_long,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Distribution of Spearman's r Coefficients",
     xlab = "r",
     ylab = "Frequency")
abline(v = 0, col = "black", lty = 2)



### Apply Fisher r-to-z transformation

list_SA_axis_corr_to_ref_long_z <- 0.5 * log((1 + list_SA_axis_corr_to_ref_long) / (1 - list_SA_axis_corr_to_ref_long))

hist(list_SA_axis_corr_to_ref_long_z,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Distribution of Fisher z-Transformed r",
     xlab = "Fisher z",
     ylab = "Frequency")
abline(v = 0, col = "black", lty = 2)




##### Splitting the demographics dataframe and grad flip variables into male and female (for by sex analyses)

# Split dataframe
abcd_demo_long_M <- subset(abcd_demo_long, sex == "M")
abcd_demo_long_F <- subset(abcd_demo_long, sex == "F")

# Get logical indices for male and female
male_idx <- abcd_demo_long$sex == "M"
female_idx <- abcd_demo_long$sex == "F"

# Split list_SA_axis_corr_to_ref_long_z_M using the same indices
list_SA_axis_corr_to_ref_long_z_M <- list_SA_axis_corr_to_ref_long_z[male_idx]
list_SA_axis_corr_to_ref_long_z_F <- list_SA_axis_corr_to_ref_long_z[female_idx]

str(list_SA_axis_corr_to_ref_long_z_M)
str(list_SA_axis_corr_to_ref_long_z_F)


# for descriptive info info
mean(list_SA_axis_corr_to_ref_long[male_idx])
mean(list_SA_axis_corr_to_ref_long[female_idx])

median(list_SA_axis_corr_to_ref_long[male_idx])
median(list_SA_axis_corr_to_ref_long[female_idx])




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

## Males

# select a subset of variables from the demo dataframe
covar_df_M <- abcd_demo_long_M %>%
  
  dplyr::select(src_subject_id_fmt, age_months, pds_p_score, sex, tot_SA, family_id, site_id, eventname) %>% 
  dplyr::rename("subject_id" = src_subject_id_fmt, "age" = age_months, "pds" = pds_p_score)


# set data type 
covar_df_M$age <- as.numeric(covar_df_M$age)
covar_df_M$pds <- as.numeric(covar_df_M$pds)
covar_df_M$sex <- as.factor(covar_df_M$sex)
covar_df_M$tot_SA <- as.numeric(covar_df_M$tot_SA)
covar_df_M$family_id <- as.factor(covar_df_M$family_id)
covar_df_M$site_id <- as.factor(covar_df_M$site_id)
covar_df_M$subject_id <- as.factor(covar_df_M$subject_id)
covar_df_M$eventname <- as.factor(covar_df_M$eventname)


# Nested effect of family within site
covar_df_M$site.family  <- interaction(covar_df_M$site_id, covar_df_M$family_id, drop = TRUE)

# Nested effect of subject within family within site
covar_df_M$site.family.subject <- interaction(covar_df_M$site_id, covar_df_M$family_id, covar_df_M$subject_id, drop = TRUE)



### For BAM: Determining the random effects variables (interactions) needed to code sub nested in fam nested in site

any(duplicated(paste(covar_df_M$site_id, covar_df_M$family_id)))  # TRUE
# TRUE  => family labels repeat across sites (not globally unique)
# FALSE => family labels are unique per site combination

any(duplicated(paste(covar_df_M$site_id, covar_df_M$family_id, covar_df_M$subject_id)))  # TRUE
# TRUE  => family and subject labels repeat across sites (not globally unique)
# FALSE => family and subject labels are unique per site combination

# fyi
length(unique(covar_df_M$site_id)) 
length(unique(interaction(covar_df_M$site_id, covar_df_M$family_id, drop=TRUE))) 
length(unique(interaction(covar_df_M$site_id, covar_df_M$family_id, covar_df_M$subject_id, drop=TRUE))) 


## -> Given that we get TRUE: not globally unique, we need to specify 3 random effects
# s(site_id, bs = "re") +
# s(site_id.family_id, bs = "re") +
# s(site_id.family_id.subject_id, bs = "re")

## Therefore creating the relevant interaction terms:
# Nested effect of family within site
covar_df_M$site.family  <- interaction(covar_df_M$site_id, covar_df_M$family_id, drop = TRUE)

# Nested effect of subject within family within site
covar_df_M$site.family.subject <- interaction(covar_df_M$site_id, covar_df_M$family_id, covar_df_M$subject_id, drop = TRUE)

str(covar_df_M)





## Females 

# select a subset of variables from the demo dataframe
covar_df_F <- abcd_demo_long_F %>%
  
  dplyr::select(src_subject_id_fmt, age_months, pds_p_score, sex, tot_SA, family_id, site_id, eventname) %>% 
  dplyr::rename("subject_id" = src_subject_id_fmt, "age" = age_months, "pds" = pds_p_score)


# set data type 
covar_df_F$age <- as.numeric(covar_df_F$age)
covar_df_F$pds <- as.numeric(covar_df_F$pds)
covar_df_F$sex <- as.factor(covar_df_F$sex)
covar_df_F$tot_SA <- as.numeric(covar_df_F$tot_SA)
covar_df_F$family_id <- as.factor(covar_df_F$family_id)
covar_df_F$site_id <- as.factor(covar_df_F$site_id)
covar_df_F$subject_id <- as.factor(covar_df_F$subject_id)
covar_df_F$eventname <- as.factor(covar_df_F$eventname)

# Nested effect of family within site
covar_df_F$site.family  <- interaction(covar_df_F$site_id, covar_df_F$family_id, drop = TRUE)

# Nested effect of subject within family within site
covar_df_F$site.family.subject <- interaction(covar_df_F$site_id, covar_df_F$family_id, covar_df_F$subject_id, drop = TRUE)


### For BAM: Determining the random effects variables (interactions) needed to code sub nested in fam nested in site

any(duplicated(paste(covar_df_F$site_id, covar_df_F$family_id)))  # TRUE
# TRUE  => family labels repeat across sites (not globally unique)
# FALSE => family labels are unique per site combination

any(duplicated(paste(covar_df_F$site_id, covar_df_F$family_id, covar_df_F$subject_id)))  # TRUE
# TRUE  => family and subject labels repeat across sites (not globally unique)
# FALSE => family and subject labels are unique per site combination

# fyi
length(unique(covar_df_F$site_id)) 
length(unique(interaction(covar_df_F$site_id, covar_df_F$family_id, drop=TRUE))) 
length(unique(interaction(covar_df_F$site_id, covar_df_F$family_id, covar_df_F$subject_id, drop=TRUE))) 


## -> Given that we get TRUE: not globally unique, we need to specify 3 random effects
# s(site_id, bs = "re") +
# s(site_id.family_id, bs = "re") +
# s(site_id.family_id.subject_id, bs = "re")

## Therefore creating the relevant interaction terms:
# Nested effect of family within site
covar_df_F$site.family  <- interaction(covar_df_F$site_id, covar_df_F$family_id, drop = TRUE)

# Nested effect of subject within family within site
covar_df_F$site.family.subject <- interaction(covar_df_F$site_id, covar_df_F$family_id, covar_df_F$subject_id, drop = TRUE)

str(covar_df_F)




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM on r coefficients in Males
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#### BAM

bam_main_rcoef_M <- bam(
  list_SA_axis_corr_to_ref_long_z_M ~
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df_M,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_main_rcoef_M <- summary(bam_main_rcoef_M)

# Save summary
capture.output(summary_bam_main_rcoef_M, file = file.path(resdir, "summary_bam_main_rcoef_M.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_main_rcoef_M)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### Plot

### Smooth Estimates

df_age <- smooth_estimates(bam_main_rcoef_M, select = "s(age)")
df_pds <- smooth_estimates(bam_main_rcoef_M, select = "s(pds)")

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
  scale_x_continuous(
    breaks = seq(108, 192, by = 12),          # ticks show 9 to 15 years - in 12 months
    limits = c(106, 190),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on similarity to adult S-A axis\n(Δ Fisher z)",
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
    y = "Partial effect on similarity to adult S-A axis\n(Δ Fisher z)",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )


ggsave(paste(resdir_fig, 'bam_smes_main_age_on_rcoef_M.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_main_pds_on_rcoef_M.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)


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

ggsave(paste(resdir_fig, 'bam_smes_panel_main_ageEpds_on_rcoef_M.svg', sep = ''), plot = plot_smes_panel, width = 11, height = 5, dpi = 300)



### Predicted marginal trajectories (using fitted_values() from gratia package)

## Age

# 1. Build prediction grid

age_grid<- seq(
  min(covar_df_M$age, na.rm=TRUE),
  max(covar_df_M$age, na.rm=TRUE),
  length.out=200
)

newdata_age <- expand.grid(
  age = age_grid,
  pds = mean(covar_df_M$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df_M$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df_M$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df_M$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df_M$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_age <- fitted_values(
  bam_main_rcoef_M, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci), 
    alpha = 0.25, 
    fill = "lightblue3", 
    color = NA
  ) +
  geom_line(
    linewidth = 1.2, 
    color = "lightblue3"
  ) +
  scale_x_continuous(
    breaks = seq(108, 192, by = 12),          # ticks show 9 to 15 years - in 12 months
    limits = c(106, 190),
    expand = expansion(mult = 0)
  ) +
  # Dual Axis Configuration - Changed sec_axis to sec.axis below
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), 
      breaks = derive() 
    )) +
  
  labs(
    x = "Age (in months)", 
    y = "Similarity to adult S-A axis (Fisher z)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_main_age_on_rcoef_M.svg', sep = ''), 
       plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


ggsave(paste(resdir_fig, 'bam_margtraj_main_age_on_rcoef_M.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


## PDS

# 1. Build prediction grid for PDS
pds_grid <- seq(
  min(covar_df_M$pds, na.rm = TRUE),
  max(covar_df_M$pds, na.rm = TRUE),
  length.out = 200
)

newdata_pds <- expand.grid(
  pds = pds_grid,
  age = mean(covar_df_M$age, na.rm = TRUE),  # hold age at mean
  tot_SA = mean(covar_df_M$tot_SA, na.rm = TRUE), # FIXED: added $tot_SA
  site_id = covar_df_M$site_id[1],
  site.family = covar_df_M$site.family[1],
  site.family.subject = covar_df_M$site.family.subject[1]
)

# 2. Generate predictions using gratia
smooth_preds_pds <- fitted_values(
  bam_main_rcoef_M, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)

# 3. Plot marginal trajectories for PDS
plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci), 
    alpha = 0.25, 
    fill = "lightblue3", 
    color = NA
  ) +
  geom_line(linewidth = 1.2, color = "lightblue3") +
  
  # Dual Axis Configuration
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), 
      breaks = derive() 
    )) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Similarity to adult S-A axis (Fisher z)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

# --- Harmonizing y-axes ---

# Find the global range for marginal trajectories
marg_ylim <- range(
  c(smooth_preds_age$.lower_ci, smooth_preds_age$.upper_ci,
    smooth_preds_pds$.lower_ci, smooth_preds_pds$.upper_ci)
)

# Update Age Plot with shared limits
plot_margtraj_age <- plot_margtraj_age + 
  scale_y_continuous(
    limits = marg_ylim, 
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)),
      breaks = derive()
    )
  )

# Update PDS Plot with shared limits
plot_margtraj_pds <- plot_margtraj_pds + 
  scale_y_continuous(
    limits = marg_ylim, 
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)),
      breaks = derive()
    )
  )

# Save the updated PDS plot
ggsave(paste(resdir_fig, 'bam_margtraj_main_pds_on_rcoef_M.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)



# Combine plots
plot_margtraj_panel <- (plot_margtraj_age + patchwork::plot_spacer() + plot_margtraj_pds) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1), guides = "collect") &   # Adjust 0.1 to increase/decrease the gap
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_panel_main_ageEpds_on_rcoef_M.svg', sep = ''), plot = plot_margtraj_panel, width = 11, height = 5, dpi = 300)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUN BAM on r coefficients in Females
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#### BAM

bam_main_rcoef_F <- bam(
  list_SA_axis_corr_to_ref_long_z_F ~
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    tot_SA +
    s(site_id, bs = "re") +
    s(site.family, bs = "re") +
    s(site.family.subject, bs = "re"),
  data = covar_df_F,
  method = "fREML",
  nthreads = 36,
  discrete = TRUE,  # key for speed/memory with large n
  select = TRUE            # allows automatic penalization of unneeded terms
)

summary_bam_main_rcoef_F <- summary(bam_main_rcoef_F)

# Save summary
capture.output(summary_bam_main_rcoef_F, file = file.path(resdir, "summary_bam_main_rcoef_F.txt"))

# Check for adequate basis dimension and model fit
mgcv::k.check(bam_main_rcoef_F)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### Plot

### Smooth Estimates

df_age <- smooth_estimates(bam_main_rcoef_F, select = "s(age)")
df_pds <- smooth_estimates(bam_main_rcoef_F, select = "s(pds)")

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
  scale_x_continuous(
    breaks = seq(108, 192, by = 12),          # ticks show 9 to 15 years - in 12 months
    limits = c(106, 190),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on similarity to adult S-A axis\n(Δ Fisher z)",
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
    y = "Partial effect on similarity to adult S-A axis\n(Δ Fisher z)",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )


ggsave(paste(resdir_fig, 'bam_smes_main_age_on_rcoef_F.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_main_pds_on_rcoef_F.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)


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

ggsave(paste(resdir_fig, 'bam_smes_panel_main_ageEpds_on_rcoef_F.svg', sep = ''), plot = plot_smes_panel, width = 11, height = 5, dpi = 300)



### Predicted marginal trajectories (using fitted_values() from gratia package)

## Age

# 1. Build prediction grid

age_grid<- seq(
  min(covar_df_F$age, na.rm=TRUE),
  max(covar_df_F$age, na.rm=TRUE),
  length.out=200
)

newdata_age <- expand.grid(
  age = age_grid,
  pds = mean(covar_df_F$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df_F$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df_F$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df_F$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df_F$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_age <- fitted_values(
  bam_main_rcoef_F, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)


# 3. Plot marginal trajectories
plot_margtraj_age <- ggplot(
  smooth_preds_age, 
  aes(x = age, y = .fitted)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci), 
    alpha = 0.25, 
    fill = "indianred3",  # Fixed fill color
    color = NA
  ) +
  geom_line(
    linewidth = 1.2, 
    color = "indianred3"   # Fixed line color
  ) +
  scale_x_continuous(
    breaks = seq(108, 192, by = 12),          # ticks show 9 to 15 years - in 12 months
    limits = c(106, 190),
    expand = expansion(mult = 0)
  ) +
  # Dual Axis Configuration
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), 
      breaks = derive() 
    )) +
  
  labs(
    x = "Age (in months)", 
    y = "Similarity to adult S-A axis (Fisher z)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )


ggsave(paste(resdir_fig, 'bam_margtraj_main_age_on_rcoef_F.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


## PDS

# 1. Build prediction grid for PDS
pds_grid <- seq(
  min(covar_df_F$pds, na.rm = TRUE),
  max(covar_df_F$pds, na.rm = TRUE),
  length.out = 200
)

newdata_pds <- expand.grid(
  pds = pds_grid,
  age = mean(covar_df_F$age, na.rm = TRUE),  # hold age at mean
  tot_SA = mean(covar_df_F$tot_SA, na.rm = TRUE),  # hold totSA at mean
  site_id = covar_df_F$site_id[1],
  site.family = covar_df_F$site.family[1],
  site.family.subject = covar_df_F$site.family.subject[1]
)

# 2. Generate predictions using gratia
smooth_preds_pds <- fitted_values(
  bam_main_rcoef_F, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"
)

# 3. Plot marginal trajectories for PDS
plot_margtraj_pds <- ggplot(
  smooth_preds_pds, 
  aes(x = pds, y = .fitted)
) +
  geom_ribbon(
    aes(ymin = .lower_ci, ymax = .upper_ci), 
    alpha = 0.25, 
    fill = "indianred3", 
    color = NA
  ) +
  geom_line(linewidth = 1.2, color = "indianred3") +
  
  # Dual Axis Configuration
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), 
      breaks = derive() 
    )) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Similarity to adult S-A axis (Fisher z)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

# --- Harmonizing y-axes ---

# Find the global range for marginal trajectories
marg_ylim <- range(
  c(smooth_preds_age$.lower_ci, smooth_preds_age$.upper_ci,
    smooth_preds_pds$.lower_ci, smooth_preds_pds$.upper_ci)
)

# Update Age Plot with shared limits
plot_margtraj_age <- plot_margtraj_age + 
  scale_y_continuous(
    limits = marg_ylim, 
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)),
      breaks = derive()
    )
  )

# Update PDS Plot with shared limits
plot_margtraj_pds <- plot_margtraj_pds + 
  scale_y_continuous(
    limits = marg_ylim, 
    sec.axis = sec_axis(
      trans = ~ ., 
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)),
      breaks = derive()
    )
  )

# Save the updated PDS plot
ggsave(paste(resdir_fig, 'bam_margtraj_main_pds_on_rcoef_F.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)



# Combine plots
# Note: guides = "none" because there's no longer a Sex legend to collect
plot_margtraj_panel <- (plot_margtraj_age + patchwork::plot_spacer() + plot_margtraj_pds) + 
  patchwork::plot_layout(widths = c(1, 0.1, 1)) & 
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_panel_main_ageEpds_on_rcoef_F.svg', sep = ''), plot = plot_margtraj_panel, width = 11, height = 5, dpi = 300)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



