---
# Age Prediction Modeling Helper
---


## Overview
<!-- badges: start -->

<!-- badges: end -->

This package provides a set of convenience functions designed to support the development of ageing-related predictive models using DNA methylation values inferred from SNP genotypes.

Direct DNA methylation profiling can be costly or unavailable in many studies. As an alternative, this package helps researchers evaluate whether SNP-predicted methylation levels at CpG sites can serve as effective substitutes for measured methylation data in building epigenetic clocks and other ageing-related models.


## Installation

You can install the development version of somepackage from <https://github.com/Kemda-Imsid/MethAMSNP.git> with:

``` r
# Pre-requisites
OSMi required these packages to be available prior insatllation
 
if (!require("BiocManager", quietly = TRUE)) ##Version>3.19
 
    install.packages("BiocManager")
 
BiocManager::install("minfi")
 
BiocManager::install("impute")
 
BiocManager::install("methylclockData")
 
BiocManager::install("methylclock")
 
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
 
devtools::install_github("pdilena/methyLImp")
``` 

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/Kemda-Imsid/MethAMSNP.git")##recommanded

install_github("https://github.com/Kemda-Imsid/MethAMSNP.git")

###You can download the .zip data from github  

#Local installation with base

install.packages("/Yourpath/MethAMSNP.zip",  repos = NULL, type = 'source')

# Local installation with devtools

library(devtools)
devtools::install("/Yourpath/MethAMSNP.zip")

```

## Main functions

``` r
library("MethAMSNP")



#read the documentation
 
?MethAMSNP::impute_missing_CpGs
```
``` r
## Missing values imputation
## basic example code

#find_missing_CpGs: A function to identify missing CpG sites required by a specific epigenetic clock
#This function compares the CpG sites present in a methylation dataset against those required by a specific epigenetic clock model (e.g., Horvath, Hannum, PhenoAge) and returns the names of CpGs that are missing. It is useful for quality control or preprocessing before applying clock models that rely on a defined set of CpG predictors.

missing_CpGs<-find_missing_CpGs(imputed_methylation_snp_beta_s,"Horvath")

#add_missing_CpGs: A function to add placeholder columns for missing CpG sites in a methylation DataFrame
#This function ensures that a methylation DataFrame includes a complete set of CpG sites by adding empty (e.g., NA-filled) columns for CpGs that are missing from the dataset.

data<-add_missing_CpGs(imputed_methylation_snp_beta_s,missing_CpGs[[1]])

#find_nearest_CpGs A function to identify the nearest neighboring CpG site (based on genomic position) for a missing CpG in a given "clock" (a set of CpG markers),

data<-find_nearest_CpGs(data, imputed_methylation_snp_beta_s,missing_CpGs[[1]])

#replace_missings_CpG: A function to impute missing methylation values using neighboring CpG positions
#Specifically, it replaces missing values using the methylation value of the nearest neighboring CpG site on the same chromosome and strand, based on genomic position. This preserves local methylation patterns and is particularly useful in analyses where spatial continuity of methylation is important, such as regional methylation profiling or epigenetic landscape reconstruction.

data<-replace_missings_CpGs(data,missing_CpGs)
```
``` r
## Modeling 

#Task 1: Horvath age estimation

# DNA_age estimate DNA methylation age using a methylation clock.

dat <- DNA_age(dat, c("Horvath"))

#The DNA_age_metrics function computes metrics to evaluate the accuracy of DNA methylation age estimates

DNA_age_metrics<-DNA_age_metrics(dat)

#The DNA_age_plot function generates a scatter plot comparing chronological age with DNA methylation-estimated biological age. 

DNA_age_plot<-DNA_age_plot(dat,"Chroage")

###Save the plot
png(filename = "snp_horvath_age.png", width = 90, height = 90,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot[[1]],ncol =1,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


##Task 2 A: Model building using lasso on a data set with measured methylation values

dat_test_prep<-data_prep(dat_test,"Chroage")

dat_test_age<-E_model(dat_test_prep[["x_train_set"]], dat_test_prep[["y_train_set"]],dat_test_prep[["x_test_set"]],
              dat_test_prep[["y_test_set"]],dat_test_prep[["tmp"]],1)

png(filename = "dat_test_72k_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(dat_test_age[[5]])
dev.off()

#retrieve the metrics for lambda se

DNA_age_metrics_dat_test<-DNA_age_metrics(dat_test_age[[2]])

# Scatter plot Chronological vs Biological age

DNA_age_plot_dat_test<-DNA_age_plot(dat_test_age[[2]],"Chroage")

# Test the model on snp data, data set dat

y_predicted <- predict(dat_test_age[[4]],newx =  x_test_set)

dat_age<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))

DNA_age_plot_dat_dat_test<-DNA_age_plot(dat_age,"Chroage")

##Task 2 B: Model building  using preselected CpGs identified through a generic algorithm and followed by fitting a simple linear regression model.




#Task 3: Model building was performed using a dataset containing measured blood methylation values. The dataset was reduced to include only CpG sites common with the predicted values dataset, filtered for those with an absolute correlation coefficient (|r|) ≥ 0.8.

dat_test_cor80_age<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1,tmp)

##retrieve the metrics for lamda se

DNA_age_metrics_dat_test_cor80_age<-DNA_age_metrics(dat_test_cor80_agee[[2]])

```







