
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  

# Project: Development

# Data cleaning for supplementary analyses: computing BMI z-scores based on CDC 2000 Growth Charts

# see https://github.com/CDC-DNPAO/CDCAnthro for package instructions and info

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

# CDC ANTHROpometry values- to generate sex- and age-standardized BMI metrics from the 2000 CDC growth charts
library(cdcanthro)  

### Clear environment
rm(list = ls())


#### set up directories
codedir = dirname(getActiveDocumentContext()$path)  # get path to current script
datadir = '/data/pt_02667/data/ABCD/'
datadir_local = '/data/p_02667/development/data/'


#### set directory to path of current script
setwd(codedir) 



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Demographics data (i.e., sex)

abcd_p_demo = read.csv(paste(datadir, 'abcd-data-release-5.1/core/abcd-general/abcd_p_demo.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')
str(abcd_p_demo)


### Anthropometrics data

abcd_y_anthro = read.csv(paste(datadir, 'abcd-data-release-5.1/core/physical-health/ph_y_anthro.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')
str(abcd_y_anthro)


### Longitudinal tracking data

abcd_y_lt = read.csv(paste(datadir, 'abcd-data-release-5.1/core/abcd-general/abcd_y_lt.csv', sep = ''), header = TRUE, fileEncoding = 'UTF-8-BOM')
str(abcd_y_lt)




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Cleaning data into one dataframe
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



### Demographics ----

# take one entry per subject (baseline only)
abcd_p_demo_baseline <- abcd_p_demo %>%
  filter(eventname == "baseline_year_1_arm_1")

# keep only subject ID
abcd_demo_clean_baseline <- abcd_p_demo_baseline %>%
  select(src_subject_id)

# clean SEX variable
# demo_sex_v2:
# 1 = Male; 2 = Female; 3 = Intersex-Male; 4 = Intersex-Female;
# 999 = Don't know; 777 = Refuse to answer
abcd_demo_clean_baseline <- abcd_demo_clean_baseline %>%
  mutate(
    sex = case_when(
      abcd_p_demo_baseline$demo_sex_v2 == 1 ~ "M",
      abcd_p_demo_baseline$demo_sex_v2 == 2 ~ "F",
      abcd_p_demo_baseline$demo_sex_v2 == 3 ~ "Intersex-M",
      abcd_p_demo_baseline$demo_sex_v2 == 4 ~ "Intersex-F",
      TRUE ~ "n/a"
    )
  )


### BMI ----

# keep only baseline, 2y, and 4y follow-ups
abcd_y_anthro_clean <- abcd_y_anthro %>%
  filter(eventname %in% c(
    "baseline_year_1_arm_1",
    "2_year_follow_up_y_arm_1",
    "4_year_follow_up_y_arm_1"
  )) %>%
  select(
    src_subject_id,
    eventname,
    anthroheightcalc,
    anthroweightcalc
  ) %>%
  # convert units: lb -> kg, inches -> meters
  mutate(
    weight_kg = anthroweightcalc * 0.45359237,
    height_m = anthroheightcalc * 0.0254,
    BMI = weight_kg / (height_m^2)
  )


### Age in months ----
# Note from cdcanthro: age: age in months specified as accurately as possible. 
# If age is given as the completed number of months (as in NHANES), add 0.5. 
# If age is given in days, divide by 30.4375.

abcd_y_lt_clean <- abcd_y_lt %>%
  filter(eventname %in% c(
    "baseline_year_1_arm_1",
    "2_year_follow_up_y_arm_1",
    "4_year_follow_up_y_arm_1"
  )) %>%
  select(src_subject_id, eventname, interview_age) %>%
  # adding 0.5 to age in mmonths (assuming that ABCD recorded completed number of months)
  mutate(
    age = interview_age + 0.5)


### Merge all dataframes ----

abcd_merged <- abcd_demo_clean_baseline %>%
  left_join(abcd_y_anthro_clean,
            by = "src_subject_id") %>%
  left_join(abcd_y_lt_clean,
            by = c("src_subject_id", "eventname"))



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Compute BMIz scores using cdcanthro

# The calculation of BMI z-scores for children without obesity is Z = (((BMI / M) ^ L) -1) / (L*S) where 
# BMI is the child’s BMI, L is Box-Cox transformation for normality for the child’s sex and age, M is median, 
# and S is coefficient of variation. Reference data are the merged LMS files at 
# https://www.cdc.gov/growthcharts/percentile_data_files.htm (Centers for Disease Control and Prevention (CDC),
# 2022). Values of sigma for children with obesity are based on formulas in the Wei et al. (2020) paper.

# Output:
# bmip and bmiz: CDC BMI percentile and z-score. These are based on the LMS method for children without obesity 
# and the 'extended' method for children with obesity. See Wei et al. (2020) for the 'extended' method based on 
# modeling high BMIs as a half-normal distribution. Extended BMIz is obtained by taking the inverse CDF of a normal 
# distribution (qnorm(value)). If 'value' is extremely close to 1, such as 1 - 1e-17, the result will be 'Inf'. 
# The function converts these values to a z-score of 8.21.

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# cdcanthro(data, age = age_in_months, wt = weight_kg, ht = height_cm, bmi = bmi, all = FALSE)
results = cdcanthro(abcd_merged, age = age, wt = weight_kg, ht = height_m, bmi = BMI)

# convert the output data.table to dataframe
results = setDF(results)

str(results)



# clean dataframe 
abcd_bmi <- results %>%
  select(src_subject_id, sex, eventname, anthroheightcalc, anthroweightcalc, 
         weight_kg, height_m,  BMI, bmiz, mod_bmiz) %>%
  dplyr::rename(BMIz = bmiz, BMIz_mod = mod_bmiz)


# export dataframe as csv
write.csv(abcd_bmi, paste(datadir_local, 'abcd_bmi.csv', sep = ''), row.names = FALSE)
