---
# Ageing Model Helper
---


## Overview
<!-- badges: start -->

<!-- badges: end -->

This package provides a set of convenience functions designed to support the development of ageing-related predictive models using DNA methylation values inferred from SNP genotypes.

Direct DNA methylation profiling can be costly or unavailable in many studies. As an alternative, this package helps researchers evaluate whether SNP-predicted methylation levels at CpG sites can serve as effective substitutes for measured methylation data in building epigenetic clocks and other ageing-related models.


## Installation

You can install the development version of somepackage from <https://github.com/Kemda-Imsid/MethAMSNP.git> with:

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
 
### basic example code
#find_missing_CpGs: A function to identify missing CpG sites required by a specific epigenetic clock
#This function compares the CpG sites present in a methylation dataset against those required by a specific epigenetic clock model (e.g., Horvath, Hannum, PhenoAge) and returns the names of CpGs that are missing. It is useful for quality control or preprocessing before applying clock models that rely on a defined set of CpG predictors.

missing_CpGs<-find_missing_CpGs(imputed_methylation_snp_beta_s,"Horvath")

#add_missing_CpGs: A function to add placeholder columns for missing CpG sites in a methylation DataFrame
#This function ensures that a methylation DataFrame includes a complete set of CpG sites by adding empty (e.g., NA-filled) columns for CpGs that are missing from the dataset.

data<-add_missing_CpGs(imputed_methylation_snp_beta_s,missing_CpGs[[1]])

#replace_missings_CpG: A function to impute missing methylation values using neighboring CpG positions
#Specifically, it replaces missing values using the methylation value of the nearest neighboring CpG site on the same chromosome and strand, based on genomic position. This preserves local methylation patterns and is particularly useful in analyses where spatial continuity of methylation is important, such as regional methylation profiling or epigenetic landscape reconstruction.

data<-replace_missings_CpGs(data,missing_CpGs)

```



