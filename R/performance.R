
###This is an usage example of the functions implemented in this package
### We are using a data set of ...

##Make data ready to impute age

#dat<-imputed_methylation_snp_horvath_beta ###CpGs in rows and sample in col
#meta_dat<-read_xlsx("yourMetadata")

sum(colnames(dat)%in%meta_dat$SIC)
colnames(dat)<-gsub("X","",colnames(dat))
meta_dat<-meta_dat[meta_dat$SIC%in%colnames(dat),]
dat<-dat[,meta_dat$SIC]
dat["Chroage",]<-meta_dat$ADULT_PROB_AGE

##Task 1

dat <- DNA_age(dat, c("Horvath"))

DNA_age_metrics<-DNA_age_metrics(dat)

DNA_age_plot<-DNA_age_plot(dat,"Chroage")

png(filename = "snp_horvath_age.png", width = 90, height = 90,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot[[1]],ncol =1,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()

