
# This is an usage example of the functions implemented in this package
# The data set contains predicted methylation values from snps genotypes

# Make data ready to impute age

# dat: predicted blood methylation values structured with CpG sites as columns and samples as rows. Imputation driven by the life_adlut biobank
# dat_test: measured blood methylation values structured with CpG sites as columns and samples as rows. Data set was retrieved from the open source Ewas datahub (https://ngdc.cncb.ac.cn/ewas/datahub/download)
# meta_dat: meta data of dat set dat containing chronological age informations
# meta_dat_test: meta data of dat_test set dat containing chronological age informations
# dat_predict_cor:  data frame containing the correlation coefficients between predicted and measured CpG methylation values.

# Prepare dat data set

sum(rownames(dat)%in%meta_dat$SIC)
rownames(dat)<-gsub("X","",rownames(dat))
meta_dat<-meta_dat[meta_dat$SIC%in%rownames(dat),]
dat<-dat[meta_dat$SIC,]
dat[,"Chroage"]<-meta_dat$ADULT_PROB_AGE

# Prepare dat_test data set

meta_dat_test<-meta_dat_test[!is.na(meta_dat_test$age),]
dat_test<-data.frame(t(dat_test[rownames(dat_test)%in%colnames(dat),meta_dat_test$sample_id])) ##common names to both, dat and dat_test
dat_test[]<-lapply(dat_test,as.numeric)
dat_test<- dat_test[, colMeans(is.na(dat_test)) <= 0.8]
dat_test<-as.matrix(dat_test)
dat_test<-impute.knn(dat_test)$data
dat_test<-dat_test[meta_dat_test$sample_id,]
dat_test<-cbind(dat_test,"Chroage"=meta_dat_test$age)


# Task 1: Horvath age estimation

dat<-t(dat)
dat <- DNA_age(dat, c("Horvath"))

DNA_age_metrics<-DNA_age_metrics(dat)

DNA_age_plot_horvath<-DNA_age_plot(dat,"Chroage")

png(filename = "snp_horvath_age.png", width = 90, height = 90,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot[[1]],ncol =1,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


# Task 2: Model building using a data set with measured methylation values

set.seed(23060830)

x<-rep(1:round(nrow(dat_test)*0.8))
nfold<-createFolds(1:round(nrow(dat_test)*0.8), k = 10, list = TRUE, returnTrain = FALSE)
tmp<-1:round(nrow(dat_test)*0.8)
for(k in 1:10){tmp[nfold[[k]]]<-k}
train.names<-sample(rownames(dat_test), round(nrow(dat_test)*0.8))
train_set<-dat_test[train.names,]
test_set<-dat_test[!rownames(dat_test)%in%train.names,]

y_train_set<-train_set[,colnames(train_set)=="Chroage"]
x_train_set<-as.matrix(train_set[,!colnames(train_set)=="Chroage"])


y_test_set<-test_set[,colnames(test_set)=="Chroage"]
x_test_set<-as.matrix(test_set[,!colnames(test_set)=="Chroage"])

dat_test_age<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1,tmp)


png(filename = "dat_test_72k_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(dat_test_age[[5]])
dev.off()


# Retrieve the metrics for lambda se

DNA_age_metrics_dat_test<-DNA_age_metrics(dat_test_age[[2]])

# Scatter plot

dat_list<-list("dat_test_age"=dat_test_age[[2]])

DNA_age_plot_dat_test<-DNA_age_plot(dat_list,"Chroage",c(paste( "LASSO","*\"  (EWAS)\"")))


# Apply the model on dat data set

dat_common_CpGs_with_dat_test<-data.frame(dat[,colnames(dat_test)])
y_test_set<-dat_common_CpGs_with_dat_test[,colnames(dat_common_CpGs_with_dat_test)=="Chroage"]
x_test_set<-as.matrix(dat_common_CpGs_with_dat_test[,!colnames(dat_common_CpGs_with_dat_test)=="Chroage"])

y_predicted <- predict(dat_test_age[[4]],newx =  x_test_set)
dat_age<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))

dat_list<-list("dat_test_age"=dat_test_age[[2]],"dat_age"=dat_age)

DNA_age_plot_dat<-DNA_age_plot(dat_list,"Chroage",
                  c(paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\"")))




#Task 3: Model building was performed using a dataset containing measured blood methylation values.
#The dataset was reduced to include only CpG sites common with the predicted values dataset, filtered for those with an absolute correlation coefficient (|r|) ≥ 0.8.

dat_predict_cor<-subset(dat_predict_cor, probeID %in% colnames(dat_test) & r >= 0.8)
dat_test_cor80<-dat_test[,c("Chroage",dat_predict_cor$probeID)]

set.seed(23060830)

x<-rep(1:round(nrow(dat_test_cor80)*0.8))
nfold<-createFolds(1:round(nrow(dat_test_cor80)*0.8), k = 10, list = TRUE, returnTrain = FALSE)
tmp<-1:round(nrow(dat_test_cor80)*0.8)
for(k in 1:10){tmp[nfold[[k]]]<-k}
train.names<-sample(rownames(dat_test_cor80), round(nrow(dat_test_cor80)*0.8))
train_set<-dat_test_cor80[train.names,]
test_set<-dat_test_cor80[!rownames(dat_test_cor80)%in%train.names,]

y_train_set<-train_set[,colnames(train_set)=="Chroage"]
x_train_set<-as.matrix(train_set[,!colnames(train_set)=="Chroage"])


y_test_set<-test_set[,colnames(test_set)=="Chroage"]
x_test_set<-as.matrix(test_set[,!colnames(test_set)=="Chroage"])

dat_test_cor80_age<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1,tmp)

##retrieve the metrics for lambda se

DNA_age_metrics_dat_test_cor80_age<-DNA_age_metrics(dat_test_cor80_age[[2]])
dat_list<-list("dat_test_age"=dat_test_age[[2]],"dat_age"=dat_age[[2]],"dat_test_cor80_age"=dat_test_cor80_age[[2]])

DNA_age_plot_dat_dat_test_cor80_age<-DNA_age_plot(dat_list,"Chroage",
                                    c(paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\""),paste( "LASSO","*\"  (LA)\"")))
# Apply the model on dat

dat_cor80<-dat[,c("Chroage",dat_predict_cor$probeID)]

y_test_set<-dat_cor80[,colnames(dat_cor80)=="Chroage"]
x_test_set<-as.matrix(dat_cor80[,!colnames(dat_cor80)=="Chroage"])

y_predicted <- predict(dat_test_age[[4]],newx =  x_test_set)
dat_cor80<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))

dat_list<-list("dat_test_age"=dat_test_age[[2]],"dat_age"=dat_age[[2]],
               "dat_test_cor80_age"=dat_test_cor80_age[[2]],"dat_cor80_age"=dat_cor80_age[[2]])

p<-DNA_age_plot(dat_list,"Chroage")



