
# This is an usage example of the functions implemented in this package
# The data set contains predicted methylation values from snps genotypes

# Make data ready to impute age

# meth_pred: data frame containg predicted blood methylation values structured with CpG sites as columns and samples as rows. Imputation driven by the life_adult biobank
# meth_obs: data frame containing measured blood methylation values structured with CpG sites as columns and samples as rows. data set was retrieved from the open source Ewas datahub (https://ngdc.cncb.ac.cn/ewas/datahub/download)
# meth_pred_cor_CpG: A data frame with one correlation coefficient per CpG site
# meth_pred_cor_sample: A data frame with one correlation coefficient per sample
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

##Validation of DNA-methylation prediction

meth_pred_cor_sample<-data.frame(vgl[["cor_ind"]])
colnames(meth_pred_cor_sample)<-"r"
meth_pred_cor_sample$r<-as.numeric(meth_pred_cor_sample$r)
meth_pred_cor_CpG<-data.frame(vgl[["cor_transkript"]])
colnames(meth_pred_cor_CpG)<-"r"
meth_pred_cor_CpG$r<-as.numeric(meth_pred_cor_CpG$r)
mean_value <-mean(meth_pred_cor_CpG$r,na.rm = TRUE)
median_value <-median(meth_pred_cor_CpG$r,na.rm = TRUE)



# Open a PNG file with width = 150mm (convert mm to inches: 1 inch = 25.4 mm)
png(filename = "snp_imputation97_hist.png", width = 170, height = 90,units = "mm", res=300)

# Layout: 2 rows, 2 columns
par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))  # margins: bottom, left, top, right

# Histogram 1
hist_data_cpgs <- hist(meth_pred_cor_CpG$r, breaks = 7000,main = "Pearson correlation methyl-probes",
                  xlab = "r",
                  density = NULL, angle = 45, border = NULL,panel.first = grid(nx = NA, ny = NULL, col = "white", lty = "dashed"),
                  font.lab = 2,
                  cex.main=0.83, #change font size of title
                  cex.sub=0.83, #change font size of subtitle
                  cex.lab=0.83, #change font size of axis labels
                  cex.axis=0.83) #change font size of axis text )

# Draw a rectangle behind the bars for panel background
usr <- par("usr")  # get plot limits
rect(usr[1], usr[3], usr[2], usr[4], col = "#ebebeb", border = NA)
grid(nx = NULL, ny = NULL, col = "white", lty = "solid", lwd = 1)
# Redraw the histogram bars on top
plot(hist_data, col = "gray70", border = "black", add = TRUE)


# boxplot(meth_pred_cor_CpG$r, breaks = 7000,main = "",
#         density = NULL, angle = 45, border = NULL,horizontal = TRUE,xlab = "r",
#
#         cex.main=0.83, #change font size of title
#         cex.sub=0.83, #change font size of subtitle
#         cex.lab=0.83, #change font size of axis labels
#         cex.axis=0.83) #change font size of axis text )


# Histogram 2


hist_data <- hist(meth_pred_cor_sample$r, plot = FALSE)

# Draw a colored background first
plot(
  hist_data,
  col = "gray70",      # bar fill color
  border = "black",    # bar borders
  main = "Pearson correlation individuals",
  xlab = "r",
  font.lab = 2,
  cex.main = 0.83,
  cex.lab = 0.83,
  cex.axis = 0.83,
  axes = TRUE,
  frame.plot = TRUE
)

# Draw a rectangle behind the bars for panel background
usr <- par("usr")  # get plot limits
rect(usr[1], usr[3], usr[2], usr[4], col = "#ebebeb", border = NA)
grid(nx = NULL, ny = NULL, col = "white", lty = "solid", lwd = 1)
# Redraw the histogram bars on top
plot(hist_data, col = "gray70", border = "black", add = TRUE)

# Boxplot 2
# boxplot(meth_pred_cor_sample$r, main = "",
#         density = NULL, angle = 45, border = NULL,horizontal = TRUE,xlab = "r",
#
#         cex.main=0.83, #change font size of title
#         cex.sub=0.83, #change font size of subtitle
#         cex.lab=0.83, #change font size of axis labels
#         cex.axis=0.83) #change font size of axis text )

dev.off()


# Task 2 Horvath age estimation

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


# Task 3: Model building using a data set with measured methylation values

set.seed(23060830)

DNA_age_plot_meth_pred_meth_obs_GA<-Lasso_GA_lm_model(meth_obs,meth_pred,"Chroage",maxIter=3000,nPop=1000)


DNA_age_plot_meth_pred_meth_data<-DNA_age_plot(DNA_age_plot_meth_pred_meth_obs_GA[["meth_pred_list"]],"Chroage",
                                               c(paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\""),
                                                 paste( "GA","*\"  (EWAS)\""),paste( "GA","*\"  (LA)\"")))

png(filename = "meth_obs_72k_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(meth_obs_age[[7]])
dev.off()

png(filename = "res_blood_snp_CpGs_Blood_72k_116CpGs.png", width = 170, height = 170,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot_meth_pred_meth_data[[1]],DNA_age_plot_meth_pred_meth_data[[2]],
                   DNA_age_plot_meth_pred_meth_data[[3]],DNA_age_plot_meth_pred_meth_data[[4]],ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


#Task 4: Model building was performed using a dataset containing measured blood methylation values.
#The dataset was reduced to include only CpG sites common with the predicted values dataset, filtered for those with an absolute correlation coefficient (|r|) ≥ 0.8.

meth_pred_cor<-subset(meth_pred_cor, probeID %in% colnames(meth_obs) & r >= 0.8)
meth_obs_cor80<-meth_obs[,c("Chroage",meth_pred_cor$probeID)]

set.seed(23060830)


DNA_age_plot_meth_pred_meth_obs_GA<-Lasso_GA_lm_model(meth_obs_cor80,meth_pred,"Chroage",maxIter=3000,nPop=1000)

DNA_age_plot_meth_pred_meth_data<-DNA_age_plot(DNA_age_plot_meth_pred_meth_obs_GA[["meth_pred_list"]],"Chroage",
                                               c(paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\""),
                                                 paste( "GA","*\"  (EWAS)\""),paste( "GA","*\"  (LA)\"")))
png(filename = "meth_obs_750_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(meth_obs_age[[7]])
dev.off()

png(filename = "res_blood_snp_CpGs_Blood_750_53CpGs.png", width = 170, height = 170,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot_meth_pred_meth_data[[1]],DNA_age_plot_meth_pred_meth_data[[2]],
                   DNA_age_plot_meth_pred_meth_data[[3]],DNA_age_plot_meth_pred_meth_data[[4]],ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


#Task 5: Model building was performed using a dataset containing from snp predicted blood methylation values.


set.seed(23060830)


DNA_age_plot_meth_pred_meth_pred<-Lasso_GA_lm_model(meth_pred,meth_pred,"Chroage",maxIter=3000,nPop=1000)

DNA_age_plot_meth_pred_meth_pred_data<-DNA_age_plot(DNA_age_plot_meth_pred_meth_pred[["meth_pred_list"]],"Chroage",
                                               c(paste( "LASSO","*\"  (LA)\""),paste( "LASSO","*\"  (LA)\""),
                                                 paste( "GA","*\"  (LA)\""),paste( "GA","*\"  (LA)\"")))
png(filename = "DNA_age_plot_meth_pred_meth_predmodelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
plot(meth_obs_age[[7]])
dev.off()

png(filename = "DNA_age_plot_meth_pred_meth_pred.png", width = 170, height = 170,units = "mm", res=300)
cowplot::plot_grid(DNA_age_plot_meth_pred_meth_pred_data[[1]],DNA_age_plot_meth_pred_meth_pred_data[[2]],
                   DNA_age_plot_meth_pred_meth_pred_data[[3]],DNA_age_plot_meth_pred_meth_pred_data[[4]],ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()


