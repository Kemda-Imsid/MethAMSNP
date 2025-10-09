
# This is an usage example of the functions implemented in this package
# The data set contains predicted methylation values from snps genotypes

# Make data ready to impute age

# meth_pred: data frame containg predicted blood methylation values structured with CpG sites as columns and samples as rows. Imputation driven by the life_adult biobank
# meth_obs: data frame containing measured blood methylation values structured with CpG sites as columns and samples as rows. data set was retrieved from the open source Ewas datahub (https://ngdc.cncb.ac.cn/ewas/datahub/download)
# meth_pred_cor:  data frame containing the correlation coefficients between predicted and measured CpG methylation values.
# meta_meth_pred: meta data of data set meth_pred containing chronological age informations
# meta_meth_obs: meta data of data set meth_pred containing chronological age informations

# Prepare meth_pred data set

sum(rownames(meth_pred)%in%meta_meth_pred$SIC)
rownames(meth_pred)<-gsub("X","",rownames(meth_pred))
meta_meth_pred<-meta_meth_pred[meta_meth_pred$SIC%in%rownames(meth_pred),]
meth_pred<-meth_pred[meta_meth_pred$SIC,]
meth_pred[,"Chroage"]<-meta_meth_pred$ADULT_PROB_AGE

# Prepare meth_obs data set

meta_meth_obs<-meta_meth_obs[!is.na(meta_meth_obs$age),]
meth_obs<-data.frame(meth_obs[meta_meth_obs$sample_id,colnames(meth_obs)%in%colnames(meth_pred)])##common names to both, meth_pred and meth_obs
meth_obs[]<-lapply(meth_obs,as.numeric)
meth_obs<- meth_obs[, colMeans(is.na(meth_obs)) <= 0.8]
meth_obs<-as.matrix(meth_obs)
meth_obs<-data.frame(impute.knn(meth_obs)$data)
meth_obs<-meth_obs[meta_meth_obs$sample_id,]
meth_obs<-cbind(meth_obs,"Chroage"=meta_meth_obs$age)


# Task 1: Horvath age estimation

ref450<- read.csv("Data/ref450.csv")
ref850<- read.csv("Data/ref850.csv")
df_annotation<-data.frame(rbind(cbind("chr"=ref850$chr,"pos"=ref850$pos,"strand"=ref850$strand,"Name"=ref850$Name),
                                cbind("chr"=ref450$chr,"pos"=ref450$pos,"strand"=ref450$strand,"Name"=ref450$Name)))##a little bit more than the 850
df_annotation<-df_annotation[!duplicated(df_annotation$Name),]
meth_pred<-t(meth_pred)
meth_pred<-replace_missings_CpGs(meth_pred,df_annotation, "Horvath")
meth_pred <- DNA_age(meth_pred, c("Horvath"))

DNA_age_metrics<-DNA_age_metrics(meth_pred)

DNA_age_plot_horvath<-DNA_age_plot(meth_pred,"Chroage")

png(filename = "snp_horvath_age.png", width = 90, height = 90,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot[[1]],ncol =1,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


# Task 2 A: Model building using a data set with measured methylation values

set.seed(23060830)

DNA_age_plot_meth_pred_meth_obs_GA<-Lasso_GA_lm_model(meth_obs,meth_pred,"Chroage",maxIter=3000,nPop=1000)


png(filename = "meth_obs_72k_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(meth_obs_age[[7]])
dev.off()

png(filename = "res_blood_snp_CpGs_Blood_72k_116CpGs.png", width = 170, height = 170,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot_meth_pred_meth_obs_GA[[6]][[1]],DNA_age_plot_meth_pred_meth_obs_GA[[6]][[2]],
                   DNA_age_plot_meth_pred_meth_obs_GA[[6]][[3]],DNA_age_plot_meth_pred_meth_obs_GA[[6]][[4]],ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()

#Task 3 : Model building was performed using a dataset containing measured blood methylation values.
#The dataset was reduced to include only CpG sites common with the predicted values dataset, filtered for those with an absolute correlation coefficient (|r|) ≥ 0.8.

meth_pred_cor<-subset(meth_pred_cor, probeID %in% colnames(meth_obs) & r >= 0.8)
meth_obs_cor80<-meth_obs[,c("Chroage",meth_pred_cor$probeID)]

set.seed(23060830)


DNA_age_plot_meth_pred_meth_obs_GA<-Lasso_GA_lm_model(meth_obs_cor80,meth_pred,"Chroage",maxIter=300,nPop=100)


png(filename = "meth_obs_750_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(meth_obs_age[[7]])
dev.off()

png(filename = "res_blood_snp_CpGs_Blood_750_53CpGs.png", width = 170, height = 170,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot_meth_pred_meth_obs_GA[[6]][[1]],DNA_age_plot_meth_pred_meth_obs_GA[[6]][[2]],
                   DNA_age_plot_meth_pred_meth_obs_GA[[6]][[3]],DNA_age_plot_meth_pred_meth_obs_GA[[6]][[4]],ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()





