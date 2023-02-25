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
# clump the data using stringent MR thresholds of r2 = 0.001 in 10,000kb window to ensure all SNPs are independent
# note this function uses EUR ancestry reference panels accessible via: https://mrcieu.github.io/gwasvcf/index.html
# not all of the VTE SNPs are present in the reference panel. Those SNPs which were absent were excluded
clump_data() %>%
# save the file
fwrite(., 'VTE_exp_dat.csv', row.names =F)

```
