
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Development

# S-A axis expansion analyses 
# (S-A expansion metric computed at the PARCEL-level)
# SENSITIVITY ANALYSIS - with OUTLIERS REMOVED 
# (note: outliers determined at the node level, so different sample than main dispersion analysis at network level)

# Testing for the effects of age/pds on S-A axis expansion (leg1) and of S-A axis expansion on rcoef/grad_flip (leg2) 
# including sex interaction effects

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
resdir = '/data/p_02667/development/results/supplementary/SA_expansion_node/'  # supplementary directory
resdir_fig = '/data/p_02667/development/results/supplementary/SA_expansion_node/figures/'  # supplementary directory


#### set directory to path of current script
setwd(codedir) 




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### S-A axis expansion metrics (computed at the parcel-level)

## Baseline 
mat_SA_expansion_metrics_baseline <- readMat(paste(datadir_local, 'SA_expansion_metrics_baseline_rm_outliers_suppl_SA_exp_node.mat', sep = ''))
names(mat_SA_expansion_metrics_baseline)  # print contents of matfile

SA_expansion_node_baseline <- mat_SA_expansion_metrics_baseline$SA.expansion.node.baseline  # shape: 1, 3955
SA_expansion_node_baseline <- t(SA_expansion_node_baseline)  #transposing to get the shape: 3955, 1

# sub ID 
sub_ID_baseline <- mat_SA_expansion_metrics_baseline$sub.ID.baseline   # shape: 3955, 1

dim(SA_expansion_node_baseline)
dim(sub_ID_baseline)



## 2y follow-up
mat_SA_expansion_metrics_fu2y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu2y_rm_outliers_suppl_SA_exp_node.mat', sep = ''))
names(mat_SA_expansion_metrics_fu2y)  # print contents of matfile

# Network-level
SA_expansion_node_fu2y <- mat_SA_expansion_metrics_fu2y$SA.expansion.node.fu2y  # shape: 1, 1255
SA_expansion_node_fu2y <- t(SA_expansion_node_fu2y)  #transposing to get the shape: 1255, 1

# sub ID 
sub_ID_fu2y <- mat_SA_expansion_metrics_fu2y$sub.ID.fu2y   # shape: 1255, 1

dim(SA_expansion_node_fu2y)
dim(sub_ID_fu2y)



## 4y follow-up
mat_SA_expansion_metrics_fu4y <- readMat(paste(datadir_local, 'SA_expansion_metrics_fu4y_rm_outliers_suppl_SA_exp_node.mat', sep = ''))
names(mat_SA_expansion_metrics_fu4y)  # print contents of matfile

# Network-level
SA_expansion_node_fu4y <- mat_SA_expansion_metrics_fu4y$SA.expansion.node.fu4y  # shape: 1, 912
SA_expansion_node_fu4y <- t(SA_expansion_node_fu4y)  #transposing to get the shape: 912, 1

# sub ID 
sub_ID_fu4y <- mat_SA_expansion_metrics_fu4y$sub.ID.fu4y   # shape: 912, 1

dim(SA_expansion_node_fu4y)
dim(sub_ID_fu4y)



### Demographics with covariates (with corresponding outliers removed)

# baseline
abcd_demo_baseline = read.csv(paste(datadir_local, 'abcd_demo_baseline_clean_rm_outliers_suppl_SA_exp_node.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

dim(abcd_demo_baseline)


# 2y
abcd_demo_fu2y = read.csv(paste(datadir_local, 'abcd_demo_fu2y_clean_rm_outliers_suppl_SA_exp_node.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

dim(abcd_demo_fu2y)


# 4y
abcd_demo_fu4y = read.csv(paste(datadir_local, 'abcd_demo_fu4y_clean_rm_outliers_suppl_SA_exp_node.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')

dim(abcd_demo_fu4y)



### Brain organization metrics (with corresponding outliers removed)

## baseline
# r coefs
list_SA_axis_corr_to_ref_baseline <- mat_SA_expansion_metrics_baseline$list.SA.axis.corr.to.ref.baseline  # shape: 1, 3955
list_SA_axis_corr_to_ref_baseline <- t(list_SA_axis_corr_to_ref_baseline)  #transposing to get the shape: 3955, 1

# gradient number corresponding to S-A axis
list_SA_axis_grad_num_baseline <- mat_SA_expansion_metrics_baseline$list.SA.axis.grad.num.baseline  # shape: 1, 3955
list_SA_axis_grad_num_baseline <- t(list_SA_axis_grad_num_baseline)  #transposing to get the shape: 3955, 1


## 2y follow-up
# r coefs
list_SA_axis_corr_to_ref_fu2y <- mat_SA_expansion_metrics_fu2y$list.SA.axis.corr.to.ref.fu2y  # shape: 1, 1255
list_SA_axis_corr_to_ref_fu2y <- t(list_SA_axis_corr_to_ref_fu2y)  #transposing to get the shape: 1255, 1

# gradient number corresponding to S-A axis
list_SA_axis_grad_num_fu2y <- mat_SA_expansion_metrics_fu2y$list.SA.axis.grad.num.fu2y  # shape: 1, 1255
list_SA_axis_grad_num_fu2y <- t(list_SA_axis_grad_num_fu2y)  #transposing to get the shape: 1255, 1


## 4y follow-up
# r coefs
list_SA_axis_corr_to_ref_fu4y <- mat_SA_expansion_metrics_fu4y$list.SA.axis.corr.to.ref.fu4y  # shape: 1, 912
list_SA_axis_corr_to_ref_fu4y <- t(list_SA_axis_corr_to_ref_fu4y)  #transposing to get the shape: 912, 1

# gradient number corresponding to S-A axis
list_SA_axis_grad_num_fu4y <- mat_SA_expansion_metrics_fu4y$list.SA.axis.grad.num.fu4y  # shape: 1, 912
list_SA_axis_grad_num_fu4y <- t(list_SA_axis_grad_num_fu4y)  #transposing to get the shape: 912, 1



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
SA_expansion_node_long = c(SA_expansion_node_baseline, SA_expansion_node_fu2y, SA_expansion_node_fu4y)

str(SA_expansion_node_long)


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



##### Transforming S-A axis expansion metric due to skew? -> currently haven't done it

# Plot histogram
hist(SA_expansion_node_long,
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
# RUN BAM analyses for SA_expansion_node_long
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Leg 1

# S-A axis expansion ~ age + pds + sex + totSA + (1 | sub+fam+site)

bam_SA_exp_leg1 <- bam(
  SA_expansion_node_long ~
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

ggsave(paste(resdir_fig, 'bam_margtraj_panel_main_leg1_ageEpds_on_SAexp_rm_outliers.svg', sep = ''), plot = plot_margtraj_panel, width = 11, height = 5, dpi = 300)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Leg 2

# rcoef ~ S-A axis expansion + age + pds + sex + totSA + (1 | sub+fam+site)

bam_SA_exp_leg2_rcoef <- bam(
  list_SA_axis_corr_to_ref_long_z ~
    sex +
    s(age, k = 10, fx = FALSE) +  # default
    s(pds, k = 10, fx = FALSE) +  # default
    s(SA_expansion_node_long, k = 10, fx = FALSE) +  # default
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

df_SAexp <- smooth_estimates(bam_SA_exp_leg2_rcoef, select = "s(SA_expansion_node_long)")
df_age <- smooth_estimates(bam_SA_exp_leg2_rcoef, select = "s(age)")
df_pds <- smooth_estimates(bam_SA_exp_leg2_rcoef, select = "s(pds)")


plot_smes_SAexp = ggplot(df_SAexp, aes(x = SA_expansion_node_long, y = .estimate)) +
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

plot_smes_age = ggplot(df_age, aes(x = age, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "honeydew3", alpha = 0.25) +
  geom_line(color = "honeydew4", linewidth = 1.2) +
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

ggsave(paste(resdir_fig, 'bam_smes_leg2_SAexp_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_smes_SAexp, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_leg2_age_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_leg2_pds_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)




### Predicted marginal trajectories (using fitted_values() from gratia package)

# 1. Build prediction grid
SAexp_grid<- seq(
  min(SA_expansion_node_long, na.rm=TRUE),
  max(SA_expansion_node_long, na.rm=TRUE),
  length.out=200
)

newdata_SAexp <- expand.grid(
  SA_expansion_node_long = SAexp_grid,
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
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
  aes(x = SA_expansion_node_long, y = .fitted, color = sex, fill = sex)
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
  SA_expansion_node_long = mean(SA_expansion_node_long, na.rm = TRUE),  # hold SAexp at mean
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions
smooth_preds_age <- fitted_values(
  bam_SA_exp_leg2_rcoef, 
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
  
  # Dual Axis Configuration
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., # No transformation of the scale itself (stay in alignment)
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), # Transform the labels only (and format as string with fitted paddings)
      breaks = derive() # This ensures it uses the exact same tick positions as the Z axis
    )) +
  
  labs(x = "Age (in months)", y = "Similarity to adult S-A axis (Fisher z)", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  ) 

ggsave(paste(resdir_fig, 'bam_margtraj_leg2_age_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)


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
  SA_expansion_node_long = mean(SA_expansion_node_long, na.rm = TRUE),  # hold SAexp at mean
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold totSA at mean
  site_id = covar_df$site_id[1],
  site.family = covar_df$site.family[1],
  site.family.subject = covar_df$site.family.subject[1]
)

# 2. Generate predictions using gratia
smooth_preds_pds <- fitted_values(
  bam_SA_exp_leg2_rcoef, 
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
  
  # Dual Axis Configuration
  scale_y_continuous(
    sec.axis = sec_axis(
      trans = ~ ., # No transformation of the scale itself (stay in alignment)
      name = "Equivalent r coefficient",
      labels = function(x) sprintf("%.2f", tanh(x)), # Transform the labels only (and format as string with fitted paddings)
      breaks = derive() # This ensures it uses the exact same tick positions as the Z axis
    )) +
  
  labs(x = "Pubertal stage (PDS score)", y = "Similarity to adult S-A axis (Fisher z)", color = "Sex", fill = "Sex") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5)
  )

ggsave(paste(resdir_fig, 'bam_margtraj_leg2_pds_on_rcoef_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)



##### grad_flip ~ S-A axis expansion + pds + sex + totSA + (1 | sub+fam+site)

bam_SA_exp_leg2_gradflip <- bam(
  grad_flip ~ 
    sex +
    s(age, k = 10, fx = FALSE) +  # default (including age because -although it wasnt a significant predictor in gradflip main analysis- it was for M)
    s(pds, k = 10, fx = FALSE) +  # default
    s(SA_expansion_node_long, k = 10, fx = FALSE) +  # default
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

df_SAexp <- smooth_estimates(bam_SA_exp_leg2_gradflip, select = "s(SA_expansion_node_long)")
df_age <- smooth_estimates(bam_SA_exp_leg2_gradflip, select = "s(age)")
df_pds <- smooth_estimates(bam_SA_exp_leg2_gradflip, select = "s(pds)")


plot_smes_SAexp = ggplot(df_SAexp, aes(x = SA_expansion_node_long, y = .estimate)) +
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

plot_smes_age = ggplot(df_age, aes(x = age, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se, ymax = .estimate + 2 * .se),
              fill = "honeydew3", alpha = 0.25) +
  geom_line(color = "honeydew4", linewidth = 1.2) +
  labs(
    x = "Age (in months)",
    y = "Partial effect on gradient flip (log-odds)",
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
    y = "Partial effect on gradient flip (log-odds)",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black")
  )

ggsave(paste(resdir_fig, 'bam_smes_leg2_SAexp_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_smes_SAexp, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_leg2_age_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_smes_age, width = 5, height = 5, dpi = 300)
ggsave(paste(resdir_fig, 'bam_smes_leg2_pds_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_smes_pds, width = 5, height = 5, dpi = 300)





### Predicted marginal trajectories (using fitted_values() from gratia package)

## S-A axis expansion

# 1. Build prediction grid
SAexp_grid<- seq(
  min(SA_expansion_node_long, na.rm=TRUE),
  max(SA_expansion_node_long, na.rm=TRUE),
  length.out=200
)

newdata_SAexp <- expand.grid(
  SA_expansion_node_long = SAexp_grid,
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
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
  aes(x = SA_expansion_node_long, y = .fitted, color = sex, fill = sex)
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
  pds = mean(covar_df$pds, na.rm = TRUE),  # hold PDS at mean
  SA_expansion_node_long = mean(SA_expansion_node_long, na.rm = TRUE),  # hold SAexp at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold total SA at mean
  site_id = covar_df$site_id[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family = covar_df$site.family[1],  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
  site.family.subject = covar_df$site.family.subject[1]  # Use any existing level (e.g. the first level) - NA yields errors (levels will be ignored anyway)
)


# 2. Generate the predictions (gratia will detect the binomial link)
smooth_preds_age <- fitted_values(
  bam_SA_exp_leg2_gradflip, 
  data = newdata_age, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"  # scale = "response" is usually the default for binomial, but being explicit ensures we get probabilities and correct intervals
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
  
  labs(x = "Age (in months)", y = "Probability of gradient flip", color = "Sex", fill = "Sex") +
  #scale_y_continuous(labels = scales::label_percent()) + # Optional: converts 0.5 to 50%
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(color = "black"))


ggsave(paste(resdir_fig, 'bam_margtraj_leg2_age_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_margtraj_age, width = 6, height = 5, dpi = 300)



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
  SA_expansion_node_long = mean(SA_expansion_node_long, na.rm = TRUE),  # hold SAexp at mean
  age = mean(covar_df$age, na.rm = TRUE),  # hold age at mean
  tot_SA = mean(covar_df$tot_SA, na.rm = TRUE),  # hold totSA at mean
  site_id = covar_df$site_id[1],
  site.family = covar_df$site.family[1],
  site.family.subject = covar_df$site.family.subject[1]
)

# 2. Generate predictions using gratia
smooth_preds_pds <- fitted_values(
  bam_SA_exp_leg2_gradflip, 
  data = newdata_pds, 
  exclude = c("s(site_id)", "s(site.family)", "s(site.family.subject)"),
  scale = "response"  # scale = "response" is usually the default for binomial, but being explicit ensures we get probabilities and correct intervals
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
  
  labs(x = "Pubertal stage (PDS score)", y = "Probability of gradient flip", color = "Sex", fill = "Sex") +
  #scale_y_continuous(labels = scales::label_percent()) + # Optional: converts 0.5 to 50%
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(color = "black"))

ggsave(paste(resdir_fig, 'bam_margtraj_leg2_pds_on_gradflip_rm_outliers.svg', sep = ''), plot = plot_margtraj_pds, width = 6, height = 5, dpi = 300)
