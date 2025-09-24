
###This is an usage example of the functions implemented in this package
### The data set contains predicted methylation values from snps genotypes

##Make data ready to impute age

#dat: imputed_methylation_snp_horvath_beta with CpGs in rows and sample in col
#dat_test: measured methylation values with CpGs in columns and sample in rows
#meta_dat: meta data of dat set dat containing chronological age informations
#meta_dat_test: meta data of dat_test set dat containing chronological age informations
#meta_dat<-read_xlsx("yourMetadata")

##prepare dat data set

sum(colnames(dat)%in%meta_dat$SIC)
colnames(dat)<-gsub("X","",colnames(dat))
meta_dat<-meta_dat[meta_dat$SIC%in%colnames(dat),]
dat<-dat[,meta_dat$SIC]
dat["Chroage",]<-meta_dat$ADULT_PROB_AGE

#prepare dat_test data set

meta_dat_test<-meta_dat_test[!is.na(meta_dat_test$age),]
dat_test<-dat_test[meta_dat_test$sample_id,]
dat_test["Chroage",]<-meta_dat_test$age
dat_test[]<-lapply(dat_test,as.numeric)

##Task 1: Horvath age estimation

dat <- DNA_age(dat, c("Horvath"))

DNA_age_metrics<-DNA_age_metrics(dat)

DNA_age_plot<-DNA_age_plot(dat,"Chroage")

png(filename = "snp_horvath_age.png", width = 90, height = 90,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot[[1]],ncol =1,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


##Task 2: Model building using a data set with measured methylation values

set.seed(23060830)

x<-rep(1:round(ncol(dat_test)*0.8))
nfold<-createFolds(1:round(ncol(dat_test)*0.8), k = 10, list = TRUE, returnTrain = FALSE)
tmp<-1:round(ncol(dat_test)*0.8)
for(k in 1:10){tmp[nfold[[k]]]<-k}
dat_test<-dat_test[,rownames(meta_dat_test)]
dat_test<-data.frame(t(dat_test[c("age",namescpgsnpblood),])) ##common names to both, dat and dat_test


set.seed(23060830)
train.names<-sample(rownames(dat_test), 1498)
train_set<-dat_test[train.names,]
test_set<-dat_test[!rownames(dat_test)%in%train.names,]

y_train_set<-train_set[,colnames(train_set)=="age"]
x_train_set<-as.matrix(train_set[,!colnames(train_set)=="age"])


y_test_set<-test_set[,colnames(test_set)=="age"]
x_test_set<-as.matrix(test_set[,!colnames(test_set)=="age"])

res_72k<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1)
png(filename = "dat_test_72k_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(res_72k[[5]])
dev.off()



