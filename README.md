# Bidirectional Mendelian randomisation analysis of venous thromboembolism and 18 cancers

This R code accompanies the manuscript titled: Causal relationships between risk of venous thromboembolism and 18 cancers: a bidirectional Mendelian randomisation analysis

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
European ancestry VTE exposure summary data was obtained from the VTE GWAS meta-analysis by (Thibord _et al_, 2022)[https://pubmed.ncbi.nlm.nih.gov/36154123/]
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
      format_data(., type="exposure") %>% 
# clump the data using stringent MR thresholds of r2 = 0.001 in 10,000kb window to ensure all SNPs are independent
# note this function uses EUR ancestry reference panels accessible via: https://mrcieu.github.io/gwasvcf/index.html
# not all of the VTE SNPs are present in the reference panel. Those SNPs which were absent were excluded
      clump_data()

# save the file
fwrite(VTE_exp_dat, 'VTE_exp_dat.csv', row.names =F)

# load cancer outcome data (EXAMPLE ONLY)
cancer1_outcome_dat <- fread('path/to/cancer1/GWAS_summ_stats') %>%
# rename all columns to correspond with the column names required for TwoSampleMR package
# i.e. SNP, chr, position, effect_allele, other_allele, eaf, beta, se, pval, ncase, ncontrol, samplesize, consortium, date, pmid, Phenotype
# format as outcome data using the TwoSampleMR package
  format_data(., type="outcome", snps = VTE_exp_dat$SNP)

# save the file
fwrite(cancer1_outcome_dat, 'cancer1_outcome_dat.csv', row.names=F)
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
        harmonise_data(exposure_dat = VTE_exp_dat, outcome_dat = ., action = 3)

oes <- cancer_outcomes %>% filter(outcome == 'Oesophageal cancer') %>%
# for oesophageal cancer I confirmed with study authors that all SNPs are on the 5' strand, therefore palindromic SNPs can be retained (i.e. action = 1 with the harmonise_data function)
        harmonise_data(exposure_dat = VTE_exp_dat, outcome_dat = ., action = 1)

# for all other cancers harmonise data using the default (action = 2) which uses effect allele frequencies to resolve ambiguous/palindromic SNPs 
harmonised_dat <- cancer_outcomes %>% filter(outcome != 'Glioma' & outcome != 'Oesophageal cancer') %>% 
    harmonise_data(exposure_dat = VTE_exp_dat, outcome_dat = ., action =2) %>%
# create a single dataframe including the harmonised oesophageal cancer and glioma data
    rbind(., glioma, oes)

# save
fwrite(harmonised_dat, 'harmonised_dat_VTE_to_cancer.csv', row.names = F)

```
### Steiger-filtering

Here I have used Steiger-filtering to exclude SNPs which explain more variance in the outcome (r2.outcome) than the exposure (r2.exposure) as these SNPs are likely to be invalid instruments (which either act through horizontal pleiotropy or proxy a reverse causal pathway from outcome to exposure)

This function requires prevalence estimates for each outcome and exposure in order to estimate the r2 for each SNP 
(or alternatively it defaults to a prevalence of 0.1 if no prevalence is supplied)

European prevalence data for each cancer was obtained from the (International Agency for Reseach on Cancer website)[https://gco.iarc.fr/today/online-analysis-table?v=2020&mode=cancer&mode_population=continents&population=900&populations=908&key=asr&sex=0&cancer=39&type=0&statistic=5&prevalence=0&population_group=0&ages_group%5B%5D=0&ages_group%5B%5D=17&group_cancer=0&include_nmsc=0&include_nmsc_other=1] which shows cumulative incidence statistics for each cancer.

``` r{prevalence_estimates}

IARC_dat <- fread('./IARC_cancer_prevalence_EUR.csv')

# Some of the cancer names do not match the cancer names in my harmonised data so need to change these. Note here the prevalence of Glioma, CLL and DLBCL is estimated from the stated cumulative incidence of 'Brain/CNS tumours', 'Leukaemia' and 'Non-Hodgkin lymphoma' respectively. 
IARC_dat$Cancer <- ifelse(IARC_dat$Cancer == 'Oropharynx', 'Oropharyngeal cancer',
                          ifelse(IARC_dat$Cancer == 'Lip,oralcavity', 'Oral cancer',
                                 ifelse(IARC_dat$Cancer == 'Leukaemia', 'Chronic lymphocytic leukaemia',
                                        ifelse(IARC_dat$Cancer == 'Brain,centralnervoussystem', 'Glioma',
                                               ifelse(IARC_dat$Cancer == 'Non-Hodgkinlymphoma',
                                                      'Diffuse large B cell lymphoma',
                                        IARC_dat$Cancer)))))

prev <- IARC_dat %>% dplyr::select(Cancer, `Cum. risk`) %>% mutate(prevalence.outcome = `Cum. risk`/100)

harmonised_dat <- left_join(harmonised_dat, prev, by = c('outcome' = 'Cancer')) %>%
# follicular lymphoma and MZL do not have a prevalence estimate since they are not in the IARC table. Approximate using the prevalence of DLBCL (0.0183)
  mutate(prevalence.outcome = ifelse(is.na(prevalence.outcome), 0.0183, prevalence.outcome)) %>%
# Since VTE is an acute event, the estimated incidence of 2 per 1000 person years has been used in this formula. 
  mutate(prevalence.exposure = 0.002)


# first set units columns so the Steiger-filtering function will work
harmonised_dat_steigered <- harmonised_dat %>% 
                      mutate(units.exposure = 'log odds',
                      units.outcome = 'log odds') %>%

# the function also requires effect allele frequencies 
#Glioma and Oesophageal cancer datasets did not have effect allele frequencies available (so missing eaf.outcome. Since the data has already been harmonised and palindromic SNPs discarded as appropriate, I will approximate the eafs for these cancers using the VTE (outcome effect allele frequencies - eaf.exposure)
                    mutate(eaf.outcome = ifelse(is.na(eaf.outcome), eaf.exposure, 
                           eaf.outcome) %>%
# perform the Steiger filtering
                    steiger_filtering %>%
# for the Steiger 'FALSE' SNPs (where r2.outcome > r2.exposure), change mr_keep to FALSE to exclude from the analysis
                    mutate(mr_keep = ifelse(steiger_dir == F,
                                      FALSE, mr_keep)
# save the file
fwrite(harmonised_dat_steigered, 'harmonised_dat_steigered_VTE_to_cancer.csv'

```

### Summary of harmonised data

To create Table 2: showing how many SNPs used for each VTE-cancer analysis

``` r{table2}

# create a vector of the individual cancer outcomes
outcomes <- cancer_outcomes %>% arrange(desc(samplesize.outcome)) %>% select(outcome) %>% distinct %>% unlist

# 1) How many VTE SNPs were missing from the cancer summary data for each cancer
snps_unavailable <- list()
for (i in outcomes) {
  eachcancer <- cancer_outcomes %>% filter(outcome == i)
  # for each cancer how many VTE SNPs are not available
  n_snp <- VTE_exp_dat %>% filter(SNP %!in% eachcancer$SNP) %>% count
  outcome <- i
  snps_unavailable[[i]] <- cbind(outcome, n_snp)
  }
snps_unavailable <- do.call(rbind, snps_unavailable) %>% rename('snps_unavailable' = 'n')

# 2) How many VTE SNPs are excluded from each VTE-cancer analysis after harmonisation due to either an ambiguous coding strand (palindromic SNPs which could not be resolved) or Steiger-filtering
snps_excluded <- list()
for (i in outcomes) {
  n_snp <- harmonised_dat_steigered %>% filter(outcome == i) %>% filter(mr_keep == F) %>% count
  outcome <- i
  snps_excluded[[i]] <- cbind(outcome, n_snp)
  }
snps_excluded <- do.call(rbind, snps_excluded) %>% 
  rename('total_snps_excluded' = 'n')
  
  # 3) How many VTE SNPs are used in each VTE-cancer analysis
  snps_used <- list()
for (i in outcomes) {
  n_snp <- harmonised_dat_steigered %>% filter(outcome == i) %>% filter(mr_keep == T) %>% count
  outcome <- i
  snps_used[[i]] <- cbind(outcome, n_snp)}
snps_used <- do.call(rbind, snps_used) %>% 
  rename('snps_used' = 'n')
  
  # create table 2
  list(snps_unavailable, snps_excluded, snps_used) %>% reduce(full_join, by = 'outcome') %>% fwrite(., 'table2.csv')
  
```
To create supplementary table 1

``` r{supptable1}
supp_table1 <- harmonised_dat_steigered %>%
  # delete redundant columns and reorder columns
  select("SNP", "chr", "position", "exposure", "outcome", "effect_allele.exposure", "other_allele.exposure","effect_allele.outcome","other_allele.outcome", "eaf.exposure","eaf.outcome", "beta.exposure", "se.exposure", "pval.exposure", "ncase.exposure",         "ncontrol.exposure", "samplesize.exposure",  "beta.outcome", "se.outcome",    "pval.outcome", "ncase.outcome", "ncontrol.outcome", "samplesize.outcome", "units.exposure", "prevalence.exposure" ,"rsq.exposure", "units.outcome", "prevalence.outcome", "rsq.outcome", "palindromic","ambiguous",           
"steiger_dir", "steiger_pval", "mr_keep") %>%
  # round numeric columns
  mutate_at(vars(rsq.outcome, steiger_pval, pval.exposure), ~formatC(., format = "e", digits = 2)) %>%
  mutate_at(vars( 'rsq.exposure', 'beta.exposure', 'se.exposure', 'eaf.exposure'), ~(round(., digits=4))) %>%
  arrange(., desc(samplesize.outcome), chr, position)
``` 

### calculate F-statistic

The mr_wrapper function returns the F statistic (under `info`)

``` r{Fstatistic}
# run it on the non-filtered Steiger dataframe
wrapper <- mr_wrapper(harmonised_dat, parameters = default_parameters())
# extract the info results results to a new list 
temp <- list()
for (i in 1:length(wrapper)) {
  temp[[i]] <- wrapper[[i]]$info
  }
# merge into a single df
info <- do.call(rbind, temp)
# assign the exposure names and outcome names to the info results
outcome_names <- dat %>% dplyr::select(id.exposure, exposure, id.outcome, outcome) %>% distinct(id.exposure, id.outcome, exposure, outcome) 

info <- left_join(info, outcome_names)
```

### Power calculations
Calculate the power using the Rcode derived from the supplement in: Burgess 2014: https://academic.oup.com/ije/article/43/3/922/761826

``` r{Power}
# first split the harmonised data by exposure (using only the SNPs that are mr_keep: SNPs that will be used in MR)
list <- harmonised_dat_steigered %>% filter(mr_keep == T) %>% group_by(outcome) %>% group_split()

# Now I want to apply a function to this list to calculate the following:
sum.rsq.exposure <- lapply(list, function(x) {sum(x$rsq.exposure)})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('sum.rsq.exposure' = 'V1')
mean.ncase.exposure <- lapply(list, function(x) {round(mean(x$ncase.exposure))})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('mean.ncase.exposure' = 'V1')
mean.ncontrol.exposure <- lapply(list, function(x) {round(mean(x$ncontrol.exposure))})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('mean.ncontrol.exposure' = 'V1')
mean.ncase.outcome <- lapply(list, function(x) {round(mean(x$ncase.outcome))})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('mean.ncase.outcome' = 'V1')
mean.ncontrol.outcome<- lapply(list, function(x) {round(mean(x$ncontrol.outcome))})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('mean.ncontrol.outcome' = 'V1')
exposure <- lapply(list, function(x) {x$exposure %>% head(1)})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('exposure' = 'V1')
outcome <- lapply(list, function(x) {x$outcome %>% head(1)})%>% 
  do.call(rbind, .) %>% as.data.frame %>% rename('outcome' = 'V1')

# merge these variables into a single df
power <- cbind(exposure, mean.ncase.exposure, mean.ncontrol.exposure, outcome,  mean.ncase.outcome, mean.ncontrol.outcome, sum.rsq.exposure) %>%
   mutate(ratio.case.control.exposure = round((mean.ncontrol.exposure/mean.ncase.exposure), digits=2),
         samplesize.exposure = mean.ncontrol.exposure+mean.ncase.exposure,
         effective.samplesize.exposure = effective_n(mean.ncase.exposure, mean.ncontrol.exposure),
         ratio.case.control.outcome = round((mean.ncontrol.outcome/mean.ncase.outcome), digits=2),
         samplesize.outcome= mean.ncontrol.outcome+mean.ncase.outcome)

The OUTCOME GWAS sample size and case control ratio is used
ratio = power$ratio.case.control.outcome # ratio of cases:controls = 1:ratio
n = power$samplesize.outcome # sample size

# calculate power to detect an odds ratio of 1.5 (for cancer for every SD increase in risk of VTE)
b1 = log(1.5) # or log of OR per SD of causal effect that I need to detect
power$power_alt_OR.1.5 <- pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))

```

### Run the MR-IVW and sensitivity analyses

``` r{results_VTE_to_cancer}

results <- mr(harmonised_dat_steigered, method_list = c('mr_ivw', 'mr_egger_regression', 'mr_weighted_median', 'mr_weighted_mode')) %>%
# convert logodds to OR and confidence intervals using the previously created function (see set_up chunk)
  generate_odds_ratios %>%
# add FDR pvalues  
  mutate(FDR_pval = p.adjust(pval, method = 'fdr')) %>% 
  mutate_if(is.numeric, ~round(., 4))

# look at pleiotropy
mr_pleiotropy <- mr_pleiotropy_test(harmonised_dat_steigered) %>% dplyr::select(exposure, outcome, egger_intercept, se, pval)

# look at heterogeneity
mr_heterogeneity <- mr_heterogeneity(
  harmonised_dat_steigered,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj
)

```
Combined results from these three analyses are shown in supplementary table 2

### Graphs

#### Forest plot of MR-IVW analysis

``` r {forest_plot}
# merge the MR results with the heterogeneity stats in order to forest plot them:
res2<- left_join(results, mr_heterogeneity) %>%
dplyr::select(outcome, method, OR, OR_lci95, OR_uci95, pval, FDR_pval, Q, Q_pval) %>% filter(method == 'Inverse variance weighted') %>% arrange(., pval) %>% 
  as.list # not essential to convert to a list but allows the next line of code to work
res2$length <- length(res2) # this allows me to plot the ilab headings on the graph

forest(x=res2$OR, ci.lb = res2$OR_lci95, ci.ub = res2$OR_uci95, 
       slab = res2$outcome, # label on left side
       cex = 1.0, # text size
       refline = 1.0, # line of null effect
       main = 'MR-IVW analysis of effect of VTE (exposure) on cancer risk (outcome)', # title for the overall graph
       header= c('Outcome', 'OR [95% CI]'),
       textpos=c(0.40, 1.6), # position of annotations
       ilab = cbind(format(round(res2$pval, digits=2), nsmall=2), format(round(res2$FDR_pval, digits=2), nsmall=2), formatC(res2$Q_pval, format = "e", digits = 1)), # extra annotations
       ilab.xpos = c(1.65, 1.75, 1.9), # position for extra annotations
       pch = 15, # shape of point
       psize = 1.0, # size to plot the points
       xlab = 'OR [95% CI] for cancer per log-odds increase in genetic risk of VTE',
       xlim = c(0.38, 1.95))
# to add labels for the ilabs
text(c(1.65, 1.75, 1.9), res2$length+10,
     c('P', 'FDR-P', 'het-P'), cex=1.0)
par(cex.lab = 3)
```

#### Scatter plots (not shown in manuscript)

``` r{scatter_plots}
scatter_plots <- mr_scatter_plot(results, harmonised_dat_steigered)
setwd('path/to/results')
pdf("./scatterplots.pdf")
for (i in 1:length(scatter_plots)) {print(scatter_plots[[i]])}
dev.off()
```

#### Single SNP, leave one out and funnel plots (supplementary figures)

``` r{sensitivity_analyses_plots}

res_single <- mr_singlesnp(harmonised_dat_steigered, all_method = 'mr_ivw')

# single SNP

single_snp_plots <- mr_forest_plot(res_single)
setwd('path/to/results')
pdf("./singleSNPplots.pdf")
for (i in 1:length(single_snp_plots)) {print(single_snp_plots[[i]])}
dev.off()

# funnel plots

funnel_plots <- mr_funnel_plot(res_single)
setwd('path/to/results')
pdf("./funnel_plots.pdf")
for (i in 1:length(funnel_plots)) {print(funnel_plots[[i]])}
dev.off()

# leave one out

res_loo <- mr_leaveoneout(harmonised_dat_steigered)
plot_loo <- mr_forest_plot(res_loo)
setwd('path/to/results')
pdf("./leave_one_out_plots.pdf")
for (i in 1:length(plot_loo)) {print(plot_loo[[i]])}
dev.off()
```

#### Sensitivity analyses with and without rs687289
The leave one out plots indicate an outlier SNP rs687289 which seems to be distorting the analysis for pancreatic, ovarian and endometrial cancer. Re-run the analysis leaving out rs687289
Draw forest plots including MR-IVW, MR-Egger, Weighted mean/mode for analyses with and without rs687289

``` r{sensitivity_analyses_plots}

res_without_rs687289 <- dat_steigered %>% filter(outcome == 'Pancreatic cancer'|outcome == 'Ovarian cancer'|outcome == 'Endometrial cancer'|outcome == 'Oral cancer') %>% filter(SNP != 'rs687289') %>% mr %>% generate_odds_ratios %>% mutate_if(is.numeric, ~round(., 4))

# Forest sensitivity with and without rs687289 using metafor package

sensitivity <- MRresults %>% 
# select the relevant columns from the overall results
  dplyr::select(outcome, method, OR, OR_lci95, OR_uci95, pval) %>% filter(outcome == 'Pancreatic cancer'|outcome == 'Oral cancer'|outcome == 'Ovarian     cancer'|outcome == 'Endometrial cancer') %>% filter(method == 'Weighted mode'|method == 'Weighted median'|method == 'MR Egger'|method == 'Inverse       variance weighted') %>%
  mutate(method = factor(method, levels = c("Weighted mode", "Weighted median", "MR Egger", "Inverse variance weighted")),
         outcome = factor(outcome, levels = c("Oral cancer", "Endometrial cancer", "Ovarian cancer","Pancreatic cancer")))%>% 
  arrange(outcome, method)

sensitivity_without_rs687289 <- res_without_rs687289 %>% 
# select the relevant columns 
dplyr::select(outcome, method, OR, OR_lci95, OR_uci95, pval) %>% filter(method == 'Weighted mode'|method == 'Weighted median'|method == 'MR Egger'|method == 'Inverse variance weighted') %>%
mutate(method = factor(method, levels = c("Weighted mode", "Weighted median", "MR Egger", "Inverse variance weighted")),
       outcome = factor(outcome, levels = c("Oral cancer", "Endometrial cancer", "Ovarian cancer","Pancreatic cancer"))) %>%
arrange(outcome, method)

# Draw forest plot
mycols = rep(c('darkred', 'green4', 'blue', 'black'), times=4)
par(mfrow = c(1,2)) # to draw two plots side by side in same window
forest(x=sensitivity$OR, ci.lb = sensitivity$OR_lci95, ci.ub = sensitivity$OR_uci95,
slab = sensitivity$method, # label on left side
cex = 1.0, # text size
refline = 1.0, # line of null effect
header= c('MR method', 'OR [95% CI]'),
col = mycols,
pch = 15, # shape of point
psize = 1.0, # size to plot the points
order = c(sensitivity$outcome), # order by outcome
rows=c(1:4, 8:11, 15:18, 22:25), # say what rows of the graph to plot each group
xlab = 'OR [95% CI] for cancer per log-odds increase in genetic risk of VTE',
main = 'Sensitivity analyses including all SNPs',
ylim = c(0, 30)) # specify ylimits otherwise labels don't fit
text(-0.3, c(5.5,12.5,19.5, 26.5), c("Oral cancer",
"Endometrial cancer",
"Ovarian cancer",
"Pancreatic cancer"), cex = 1.0, pos=4, col = 'black')
rect(-6, 0, 3.0, 6.0, col=adjustcolor("grey60", 0.1), border=NA)
rect(-6, 7, 3.0, 13.0, col=adjustcolor("grey60", 0.1), border=NA)
rect(-6, 14, 3.0, 20.0, col=adjustcolor("grey60", 0.1), border=NA)
rect(-6, 21, 3.0, 27, col=adjustcolor("grey60", 0.1), border=NA)
forest(x=sensitivity_without_rs687289$OR, ci.lb = sensitivity_without_rs687289$OR_lci95, ci.ub = sensitivity_without_rs687289$OR_uci95,
slab = sensitivity_without_rs687289$method, # label on left side
cex = 1.0, # text size
refline = 1.0, # line of null effect
header= c('MR method', 'OR [95% CI]'),
col = mycols,
pch = 15, # shape of point
psize = 1.0, # size to plot the points
order = c(sensitivity_without_rs687289$outcome), # order by outcome
rows=c(1:4, 8:11, 15:18, 22:25), # say what rows of the graph to plot each group
xlab = 'OR [95% CI] for cancer per log-odds increase in genetic risk of VTE',
main = 'Sensitivity analyses excluding rs687289',
ylim = c(0, 30)) # specify ylimits otherwise labels don't fit
text(0.15, c(5.5,12.5,19.5, 26.5), c("Oral cancer",
"Endometrial cancer",
"Ovarian cancer",
"Pancreatic cancer"), cex = 1.0, pos=4, col = 'black')
rect(0, 0, 3.0, 6.0, col=adjustcolor("grey60", 0.1), border=NA)
rect(0, 7, 3.0, 13.0, col=adjustcolor("grey60", 0.1), border=NA)
rect(0, 14, 3.0, 20.0, col=adjustcolor("grey60", 0.1), border=NA)
rect(0, 21, 3.0, 27, col=adjustcolor("grey60", 0.1), border=NA)

```

#### Sensitivity analysis with no Steiger-filtering and replicated VTE SNPs

A sensitivity analysis was performed using the same code above on the non-Steiger filtered harmonised dataframe

```r{sensitivity_non_Steigered_dat}
dat <- fread('./harmonised_dat_VTE_to_cancer.csv')

MRresults <- mr(dat) %>%
# ensure the generate odds ratios function has been created (top of the page)
  generate_odds_ratios %>%
# add FDR 
  mutate(FDR_pval = p.adjust(pval, method = 'fdr')) %>% 
  mutate_if(is.numeric, ~round(., 4))
MR_heterogeneity <- mr_heterogeneity(
  dat,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj
)
res2<- left_join(MRresults, Mr_heterogeneity) %>% 
  dplyr::select(outcome, method, OR, OR_lci95, OR_uci95, pval, FDR_pval, bonferroni_pval, Q, Q_pval) %>% filter(method == 'Inverse variance weighted') %>% arrange(., pval) %>% 
  as.list # not essential to convert to a list but allows the next line of code to work
res2$length <- length(res2) # this allows me to plot the ilab headings on the graph

forest(x=res2$OR, ci.lb = res2$OR_lci95, ci.ub = res2$OR_uci95, 
       slab = res2$outcome, # label on left side
       cex = 1.0, # text size
       refline = 1.0, # line of null effect
       main = 'MR-IVW analysis of effect of VTE (exposure) on cancer risk (outcome) - non-Steiger filtered results', # title for the overall graph
       header= c('Outcome', 'OR [95% CI]'),
       textpos=c(0.40, 1.6), # position of annotations
       ilab = cbind(format(round(res2$pval, digits=2), nsmall=2), format(round(res2$FDR_pval, digits=2), nsmall=2), formatC(res2$Q_pval, format = "e", digits = 1)), # extra annotations
       ilab.xpos = c(1.65, 1.75, 1.9), # position for extra annotations
       pch = 15, # shape of point
       psize = 1.0, # size to plot the points, otherwise they are a function of the weighting which looks odd
       xlab = 'OR [95% CI] for cancer per log-odds increase in genetic risk of VTE',
       xlim = c(0.38, 1.95))
# to add labels for the ilabs it is a bit of a faff
text(c(1.65, 1.75, 1.9), res2$length+10, 
     c('P', 'FDR-P', 'het-P'), cex=1.0)
par(cex.lab = 3)

```

A further sensitvity analysis was performed using SNPs which were identified as 'known' (i.e. well replicated) in supplementary table 2 of the publication by (Thibord _et al_, 2022)[https://pubmed.ncbi.nlm.nih.gov/36154123/]

```r{sensitivity_known_VTE_SNPs_only}
knownSNPs <- fread('./Thibord_published_supp2.csv', skip=2) %>% 
filter(`Novel/known` == 'Known') # 39 SNPs

# Look for these SNPs in the VTE_exp_data
VTE_exp_known_only <- VTE_exp_dat %>% filter(SNP %in% knownSNPs$SNP) %>% clump_data

# identify these SNPs in the Steiger filtered dataframe
dat_steigered_known <- harmonised_dat_steigered %>% filter(SNP %in% VTE_known_SNPs$SNP)

# these are shown in supplementary table 3 of the paper:
supp_table3 <- dat_steigered_known %>%
  # delete redundant columns and reorder columns
  select("SNP", "chr", "position", "exposure", "outcome", "effect_allele.exposure", "other_allele.exposure","effect_allele.outcome","other_allele.outcome", "eaf.exposure","eaf.outcome", "beta.exposure", "se.exposure", "pval.exposure", "ncase.exposure",         "ncontrol.exposure", "samplesize.exposure",  "beta.outcome", "se.outcome",    "pval.outcome", "ncase.outcome", "ncontrol.outcome", "samplesize.outcome", "units.exposure", "prevalence.exposure" ,"rsq.exposure", "units.outcome", "prevalence.outcome", "rsq.outcome", "palindromic","ambiguous",           
"steiger_dir", "steiger_pval", "mr_keep") %>%
  # round numeric columns
  mutate_at(vars(rsq.outcome, steiger_pval, pval.exposure), ~formatC(., format = "e", digits = 2)) %>%
  mutate_at(vars( 'rsq.exposure', 'beta.exposure', 'se.exposure', 'eaf.exposure'), ~(round(., digits=4))) %>%
  arrange(., desc(samplesize.outcome), chr, position)

# get the results
MRresults <- mr(dat_steigered_known) %>% generate_odds_ratios %>% 
  mutate(fdr_pval = p.adjust(pval, method = 'fdr'))  %>%
  mutate_if(is.numeric, ~round(., 4))
  
MR_heterogeneity <- mr_heterogeneity(
 dat_steigered_known,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj
) %>% arrange (., Q_pval)

mr_pleiotropy <- mr_pleiotropy_test(dat_steigered_known) %>% dplyr::select(exposure, outcome, egger_intercept, se, pval)

# forest plot
res2<- left_join(MRresults, mr_heterogeneity) %>% 
  dplyr::select(outcome, method, OR, OR_lci95, OR_uci95, pval, fdr_pval, Q, Q_pval) %>% filter(method == 'Inverse variance weighted') %>% arrange(., pval) %>% 
  as.list # not essential to convert to a list but allows the next line of code to work
res2$length <- length(res2) # this allows me to plot the ilab headings on the graph


forest(x=res2$OR, ci.lb = res2$OR_lci95, ci.ub = res2$OR_uci95, 
       slab = res2$outcome, # label on left side
       cex = 1.0, # text size
       refline = 1.0, # line of null effect
       main = 'MR-IVW analysis of effect of VTE (proxied using replicated instruments only) \n on cancer risk (outcome)', # title for the overall graph
       header= c('Outcome', 'OR [95% CI]'),
       textpos=c(0.40, 1.7), # position of annotations
       ilab = cbind(format(round(res2$pval, digits=2), nsmall=2), format(round(res2$fdr_pval, digits=2), nsmall=2), formatC(res2$Q_pval, format = "e", digits = 1)), # extra annotations
       ilab.xpos = c(1.75, 1.85, 1.95), # position for extra annotations
       pch = 15, # shape of point
       psize = 1.0, # size to plot the points, otherwise they are a function of the weighting which looks odd
       xlab = 'OR [95% CI] for cancer per log-odds increase in genetic risk of VTE',
       xlim = c(0.38, 2.0))
# to add labels for the ilabs it is a bit of a faff
text(c(1.75, 1.85, 1.95), res2$length+11,
     c('P', 'FDR-P', 'het-P'), cex=1.0)
par(cex.lab = 3)

```
#### MR Wald ratios for the association between VTE risk proxied by either Factor V Leiden (rs6025) or Prothrombin G20210A (rs1799963), and cancer


```r{FVL_PT_Waldratios}

# rs1799963 wasn't included in the main analysis as it was omitted from the LD_ref panel due to a low MAF. So I need to reharmonise the data just for these SNPs.

# identify FVL and PT variant in the VTE summary data
VTE_PT_FVL <- fread('path/to/VTE/GWAS_summ_stats') %>%
# rename all columns to correspond with the column names required for TwoSampleMR package
# i.e. SNP, chr, position, effect_allele, other_allele, eaf, beta, se, pval, ncase, ncontrol, samplesize, consortium, date, pmid, Phenotype
# select just the relevant SNPs
      filter(SNP == 'rs6025'|SNP == 'rs1799963') %>%
# format as exposure data using TwoSampleMR package
      format_data(., type="exposure")

# read in the cancer outcome data
setwd('./OUTCOME_data')
cancer_outcomes <- list.files(pattern = "*dat.csv")       # Get all file names. # use this prefix to avoid accidentally reading in the VTE file...
cancer_outcomes <- lapply(cancer_outcomes, fread) %>% reduce(full_join)

# Harmonise as previously described
# Glioma and Oesophageal cancer have missing eafs which makes harmonisation of palindromic SNPs more difficult
glioma <- cancer_outcomes %>% filter(outcome == 'Glioma') %>%
# not clear if all SNPs on same coding strand use action 3 = drop all palindromic SNPs
      harmonise_data(exposure_dat = VTE_exp_with_PT, outcome_dat = ., action = 3) 

oes <- cancer_outcomes %>% filter(outcome == 'Oesophageal cancer') %>% 
# all SNPs on forward strand - same for VTE therefore palindromic SNPs retained.
    harmonise_data(exposure_dat = VTE_PT_FVL, outcome_dat = ., action = 1)

harmonised_dat_FVL_PT <- cancer_outcomes %>% filter(outcome != 'Glioma' & outcome != 'Oesophageal cancer') %>%
# for all other cancer use action = 2: the function will try to infer ambiguous/palindromix SNP issues using allele frequencies.
      harmonise_data(exposure_dat = VTE_exp_with_PT, outcome_dat = ., action =2) %>%
      # merge with the oesophageal and glioma harmonised df
      rbind(., oes, glioma)

# get the results
resFVL <- harmonised_dat_FVL_PT %>% filter(SNP=='rs6025') %>% mr() %>% generate_odds_ratios() %>%
  mutate('FDR_pval' = p.adjust(pval, method='fdr')) %>% mutate(SNP = 'rs6025')
resPT <- harmonised_dat_FVL_PT %>% filter(SNP=='rs1799963') %>% r() %>% generate_odds_ratios() %>%
  mutate('FDR_pval' = p.adjust(pval, method='fdr')) %>% mutate(SNP = 'rs1799963')

# draw Forest plots
# insert the cancers where there was no data available for the PT G20210A mutations
missing_cancers <- resFVL  %>% filter(outcome %!in% resPT$outcome) %>% select(outcome)
resPT <- resPT %>% full_join(., missing_cancers) %>%
  # arrange in the same order as the FVL plot:
   mutate(outcome = factor(outcome, levels = c( "Colorectal cancer", "Follicular lymphoma", "Melanoma", "Oral cancer", "Chronic lymphocytic leukaemia", "Oropharyngeal cancer", "Lung Cancer", "Pancreatic cancer", "Prostate cancer", "Diffuse large B cell lymphoma", "Ovarian cancer",   "Marginal zone lymphoma", "Bladder cancer", "Breast cancer", "Endometrial cancer", "Oesophageal cancer",   "Glioma", "Kidney Cancer"))) %>% arrange(., outcome)
 
# FV Leiden Wald ratio forest plot 
   forest(x=resFVL$OR, ci.lb = resFVL$OR_lci95, ci.ub = resFVL$OR_uci95, 
       slab = resFVL$outcome, # label on left side
       cex = 1.0, # text size
       refline = 1.0, # line of null effect
       main = '[A] MR Wald ratios for effect of VTE proxied by Factor 5 Leiden (exposure) \non cancer risk (outcome)', # title for the overall graph
       header= c('Outcome', 'OR [95% CI]'),
       textpos=c(0.20, 1.80), # position of annotations
       ilab = cbind(format(round(resFVL$pval, digits=2), nsmall=2), format(round(resFVL$FDR_pval, digits=2), nsmall=2)), # extra annotations
       ilab.xpos = c(1.85, 2.0), # position for extra annotations
       pch = 15, # shape of point
       psize = 1.0, # size to plot the points
       xlab = 'OR [95% CI] for cancer per log-odds increase in risk of VTE proxied by Factor V Leiden genotype',
       xlim = c(0.1, 2.1))
# to add labels for the ilabs 
text(c(1.85, 2.0), 20, 
     c('P', 'FDR-P'), cex=1.0)

# Prothrombin variant Wald ratio forest plot 
options(na.action = "na.pass") # temporarily set na.action to pass so that it plots these
forest(x=resPT$OR, ci.lb = resPT$OR_lci95, ci.ub = resPT$OR_uci95,
       slab = resPT$outcome, # label on left side
       cex = 1.0, # text size
       refline = 1.0, # line of null effect
       main = '[B] MR Wald ratios for effect of VTE proxied by Prothrombin G20210A (exposure) \non cancer risk (outcome)', # title for the overall graph
       header= c('Outcome', 'OR [95% CI]'),
       textpos=c(-7, 14), # position of annotations
       ilab = cbind(format(round(resPT$pval, digits=2), nsmall=2), format(round(resPT$FDR_pval, digits=2), nsmall=2)), # extra annotations
       ilab.xpos = c(15, 17), # position for extra annotations
       pch = 15, # shape of point
       psize = 1.0, # size to plot the points
       xlab = 'OR [95% CI] for cancer per log-odds increase in risk of VTE proxied by Prothrombin G20210A genotype',
       xlim = c(-7, 18))
# to add labels for the ilabs it is a bit of a faff
text(c(15, 17), 20,
     c('P', 'FDR-P'), cex=1.0)
# re-set the na.action
options(na.action = "na.omit")
```
## MR analysis of genetic risk of 18 cancers (exposures) and VTE (outcome)

Exposure data for each cancer and VTE outcome data was formatted as described above. 


### Harmonise Cancer exposure data and VTE outcome data

```r{harmonise_cancer_VTE}
# read in exposure and outcome data
# set correct wds
cancer_exp_dat <- list.files(pattern = "*clumped.csv") %>% # Get all file names
    # read in files and merge into single df
    lapply(., fread) %>%  rbindlist(., fill=T)
VTE_outcome_dat <- fread('./OUTCOME_data/outcome_dat_VTE.csv')

# Glioma and Oesophageal cancer have missing eafs which makes harmonisation of palindromic SNPs more difficult

glioma <- cancer_exp_dat %>% filter(exposure == 'Glioma') %>%
# for glioma the coding strand is not consistent for all SNPs therefore all palindromic SNPs will be excluded (i.e. action = 3 with the harmonise-data function)
        harmonise_data(exposure_dat = ., outcome_dat = VTE_outcome_dat, action = 3)

oes <- ccancer_exp_dat %>% filter(exposure == 'Oesophageal cancer') %>%
# for oesophageal cancer I confirmed with study authors that all SNPs are on the 5' strand, therefore palindromic SNPs can be retained (i.e. action = 1 with the harmonise_data function)
       harmonise_data(exposure_dat = ., outcome_dat = VTE_outcome_dat, action = 1)

# for all other cancers harmonise data using the default (action = 2) which uses effect allele frequencies to resolve ambiguous/palindromic SNPs 
harmonised_dat <- cancer_exp_dat %>% filter(exposure != 'Glioma' & exposure != 'Oesophageal cancer') %>% 
    harmonise_data(exposure_dat = ., outcome_dat = VTE_outcome_dat, action =2) %>%
# create a single dataframe including the harmonised oesophageal cancer and glioma data
    rbind(., glioma, oes)
```

### Steiger-filtering

Prevalence estimates obtained from IARC data as described above. 

```r{Steiger_cancer_VTE}

# insert prevalence estimates
prev <- IARC_dat %>% dplyr::select(Cancer, `Cum. risk`) %>% mutate(prevalence.exposure = `Cum. risk`/100)

harmonised_dat_steigered <- left_join(harmonised_dat, prev, by = c('exposure' = 'Cancer')) %>%
  # approximate prevalence of Follicular lymphoma and Marginal zone lymphoma using DLBCL data
  mutate(prevalence.exposure = ifelse(exposure == 'Follicular lymphoma'|exposure == 'Marginal zone lymphoma', 0.0183, prevalence.exposure)) %>%
  # insert VTE prevalence
  mutate(prevalence.outcome = 0.002) %>%
  # add units 
  mutate(units.exposure = 'log odds', units.outcome = 'log odds') %>%
  # Glioma and Oesophageal cancer datasets did not have effect allele frequencies available (so missing eaf.exposure. Since the data has already been harmonised and palindromic SNPs discarded as appropriate, I will approximate the eafs for these cancers using the VTE eafs
  mutate(eaf.exposure <- ifelse(is.na(eaf.exposure), eaf.outcome, eaf.exposure)) %>%
   # Perform Steiger filtering 
    steiger_filtering()
    mutate(mr_keep = ifelse(steiger_dir == F, FALSE, mr_keep)
```

Supplementary table 6

```{Supp_table6}
supp_table6 <- left_join(harmonised_dat_steigered, coordinates, by = 'SNP') %>%
  # delete redundant columns and reorder columns
  select("SNP", "chr", "position", "exposure", "outcome", "effect_allele.exposure", "other_allele.exposure","effect_allele.outcome","other_allele.outcome", "eaf.exposure","eaf.outcome", "beta.exposure", "se.exposure", "pval.exposure", "ncase.exposure",         "ncontrol.exposure", "samplesize.exposure",  "beta.outcome", "se.outcome",    "pval.outcome", "ncase.outcome", "ncontrol.outcome", "samplesize.outcome", "units.exposure", "prevalence.exposure" ,"rsq.exposure", "units.outcome", "prevalence.outcome", "rsq.outcome", "palindromic","ambiguous",           
"steiger_dir", "steiger_pval", "mr_keep") %>%
  # round numeric columns
  mutate_at(vars(rsq.outcome, steiger_pval, pval.exposure), ~formatC(., format = "e", digits = 2)) %>%
  mutate_at(vars( 'rsq.exposure', 'beta.exposure', 'se.exposure', 'eaf.exposure'), ~(round(., digits=4))) %>%
  arrange(., desc(samplesize.exposure), chr, position)

```
