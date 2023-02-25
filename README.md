# Bidirectional Mendelian randomisation analysis of venous thromboembolism and 18 cancers

This R code accompanies the manuscript titled: ..

Set libraries and generic functions

``` r{set_libraries}
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(tidyverse)
library(metafor)

`%!in%` = Negate(`%in%`)

source('./editted_funnel_plot_function.R') # the MR funnel plot function doens't produce titles so I have editted it in this script

# create function to generate odds ratio and confidence intervals from beta/log-odds
generate_odds_ratios <- function(res) {
  res$lo_ci <- res$b - 1.96 * res$`se`
  res$up_ci <- res$b + 1.96 * res$`se`
  res$OR <- exp(res$b)
  res$OR_lci95 <- exp(res$lo_ci)
  res$OR_uci95 <- exp(res$up_ci)
  return(res)

```

## MR analysis of genetic risk of VTE (exposure) and 18 cancers (outcomes)

### Format data
European ancestry VTE exposure summary data was obtained from the VTE GWAS meta-analysis by (Thibord _et al_, 2022) [https://pubmed.ncbi.nlm.nih.gov/36154123/]
Independent European ancestry VTE risk loci are found in supplementary table 2 of this paper

Cancer outcome data for each cancer was obtained from a variety of sources (see data availability supplement)

``` r{format_data}

# load VTE exposure data
VTE_exp_dat <- fread('path/to/VTE/GWAS_summ_stats') %>%
# rename all columns to correspond with the column names required for TwoSampleMR package
# i.e. SNP, chr, position, effect_allele, other_allele, eaf, beta, se, pval, ncase, ncontrol, samplesize, consortium, date, pmid, Phenotype
# check that all SNPs are associated with VTE at GWAS significant pvalue of p=5e-8 and arrange by pval
filter(pval<=5e-8) %>% arrange(.,(pval)) %>%
# format as exposure data using TwoSampleMR package
format_data(VTE_exp_dat, type="exposure") 
# clump the data using stringent MR thresholds of r2 = 0.001 in 10,000kb window to ensure all SNPs are independent
# note this function uses EUR ancestry reference panels accessible via: https://mrcieu.github.io/gwasvcf/index.html
# not all of the VTE SNPs are present in the reference panel. Those SNPs which were absent were excluded
clump_data() %>%
# save the file
fwrite(., 'VTE_exp_dat.csv', row.names =F)

# load cancer outcome data (EXAMPLE ONLY)
cancer1 <- fread('path/to/cancer1/GWAS_summ_stats') %>%
# rename all columns to correspond with the column names required for TwoSampleMR package
# i.e. SNP, chr, position, effect_allele, other_allele, eaf, beta, se, pval, ncase, ncontrol, samplesize, consortium, date, pmid, Phenotype
# format as outcome data using the TwoSampleMR package
format_data(cancer1, type="outcome", snps = VTE_exp_dat$SNP) %>%
# save the file
fwrite(., 'cancer1_outcome_dat.csv', row.names=F)
```

### Harmonise VTE exposure data and cancer outcome data

``` r{harmonise_data}

# read in the exposure and outcome data

VTE_exp_dat <- fread('VTE_exp_dat.csv')
# to read in all the cancer_outcome_dat files and unify into a single dataframe
cancer_outcomes <- list.files(pattern = "*outcome_dat.csv") %>%
lapply(., fread) %>% reduce(full_join)

# Glioma and Oesophageal cancer have missing eafs which makes harmonisation of palindromic SNPs more difficult

glioma <- cancer_outcomes %>% filter(outcome == 'Glioma') %>%
# for glioma the coding strand is not consistent for all SNPs therefore all palindromic SNPs will be excluded (i.e. action = 3 with the harmonise-data function)
harmonise_data(exposure_dat = VTE_exp_dat, outcome_dat = glioma, action = 3)

oes <- cancer_outcomes %>% filter(outcome == 'Oesophageal cancer') %>%
# for oesophageal cancer I confirmed with study authors that all SNPs are on the 5' strand (which is the same for the VTE data), therefore palindromic SNPs can be retained (i.e. action = 1 with the harmonise_data function)
harmonise_data(exposure_dat = VTE_exp_dat, outcome_dat = ., action = 1)

# for all other cancers harmonise data using the default (action = 2) which uses effect allele frequencies to resolve ambiguous/palindromic SNPs 
dat <- cancer_outcomes %>% filter(outcome != 'Glioma' & outcome != 'Oesophageal cancer') %>% 
harmonise_data(exposure_dat = VTE_exp_dat, outcome_dat = ., action =2) %>%
# create a single dataframe including the harmonised oesophageal cancer and glioma data
rbind(., glioma, oes) %>%
# save
fwrite(., 'harmonised_dat_VTE_to_cancer.csv', row.names = F)

```
### Steiger-filtering

Here I have used Steiger-filtering to exclude SNPs which explain more variance in the outcome (r2.outcome) than the exposure (r2.exposure) ans these SNPs are likely to be invalid instruments (which either act through horizontal pleiotropy or proxy a reverse causal pathway from outcome to exposure)

This function requires prevalence estimates for each outcome and exposure in order to estimate the r2 for each SNP 
(or alternatively it defaults to a prevalence of 0.1 - the same results were obtained in a sensitivity analysis where the default prevalence setting was used)

European prevalence data for each cancer was obtained from the (International Agency for Reseach on Cancer website) [https://gco.iarc.fr/today/online-analysis-table?v=2020&mode=cancer&mode_population=continents&population=900&populations=908&key=asr&sex=0&cancer=39&type=0&statistic=5&prevalence=0&population_group=0&ages_group%5B%5D=0&ages_group%5B%5D=17&group_cancer=0&include_nmsc=0&include_nmsc_other=1]

``` r{steiger_filtering}

