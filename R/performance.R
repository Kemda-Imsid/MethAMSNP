
###This is an usage example of the functions implemented in this package
### We are using a data set of ...

##Make data ready to impute age

#dat<-imputed_methylation_snp_horvath_beta ###CpGs in rows and sample in col
#meta_dat<-read_xlsx("yourMetadata")
#meta_dat<-read_xlsx("Yourpath.xlsx")
sum(colnames(dat)%in%meta_dat$SIC)
colnames(dat)<-gsub("X","",colnames(dat))
meta_dat<-meta_dat[meta_dat$SIC%in%colnames(dat),]
dat<-dat[,meta_dat$SIC]
dat["Chroage",]<-meta_dat$ADULT_PROB_AGE

##Task 1

Horvath_age <- DNA_age(dat, "Horvath")

# #sum(age_snp_h_r$id %in%rownames(PV0697_datajoin))
#
# PV0697_datajoin<-PV0697_datajoin[age_snp_h_r$id,]
#
# age_snp_h_r<-data.frame(cbind("Chroage"=PV0697_datajoin[,"ADULT_PROB_AGE"],"Bioage"=age_snp_h_r$Horvath))
# age_snp_h_r[]<-lapply(age_snp_h_r,as.numeric)
#
# #
# #
