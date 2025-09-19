
#' DNAm age estimation giving a DNA methylation clock.
#'
#' @param dat data.frame with samples in columns and CpGs in rows. A row called "Chroage" which contains the chronological age infos.
#' @param clocks the methods used for estimating DNAmAge. Currrently "Horvath", "Hannum",
#' "Levine", "BNN", "skinHorvath", "PedBE", "Wu", "TL", "BLUP", "EN" and "all" are available.
#' @returns a data.frame with the chronological and estimated age in colums, samples in rows.
#' @export

DNA_age<-function(dat,clocks){

  Chroage<-data.frame(t(dat["Chroage",]))
  dat <- dat[rownames(dat) != "Chroage", ]
  dat<-cbind("ProbeID"=rownames(dat),dat)
  dat[,-1]<-lapply(dat[,-1],as.numeric)
  dat_age<-data.frame(DNAmAge(dat,clocks = clocks))
  rownames(dat_age)<-dat_age$id
  dat_age<-dat_age[rownames(Chroage),]
  dat_age<-data.frame(cbind("Chroage"=as.numeric(Chroage[,1]),"Bioage"=as.numeric(dat_age[,-1])))
  return(dat_age)
}

DNA_age_metrics<-function(dat){

  dat$diff<-dat$Bioage-dat$Chroage
  dat<- data.frame("min_diff"=min(dat$diff),"median_diff"=median(dat$diff),"mean_diff"=mean(dat$diff),
                  "max_diff"=max(dat$diff), "sd_diff"=sd(dat$diff), "rmse"=rmse(dat$Chroage,dat$Bioage),
                   "mae"=mae(dat$Chroage,dat$Bioage),"mad"=mad(dat$Chroage,dat$Bioage))
 return(dat)
}


DNA_age_plot <- function(dat, target_col) {
  plot_list <- list()

  # Check if target column exists
  if (!(target_col %in% names(dat))) {
    stop("Target column not found in data frame.")
  }

  # Only proceed if target is numeric
  if (!is.numeric(dat[[target_col]])) {
    stop("Target column must be numeric for scatter plots.")
  }

  # Loop through other columns
  for (col_name in names(dat)) {
    # Skip the target column itself
    if (col_name == target_col) next

    # Only plot if the other column is numeric
    if (is.numeric(dat[[col_name]])) {
      p <-  ggplot(dat, aes(x = .data[[target_col]], y = .data[[col_name]])) +
        geom_point(size=0.01) +
        geom_abline(linewidth = 0.2)+
        #geom_smooth(method=lm, linewidth= 0.2,,se=FALSE,
        #  color="blue3")+
        theme(text = element_text(size = 10, face = "bold"),legend.position = "none",legend.text=element_text(size=10)
        )+
        xlim(0,100)+
        ylim(0,100)+
        xlab("Chronological age (years)")+
        ylab("Biological age (years)")

      # Add to list
      plot_list[[col_name]] <- p
    }
  }

  return(plot_list)
}








# #####
# set.seed(23060830)
# x<-rep(1:1498)
# nfold<-createFolds(1:1498, k = 10, list = TRUE, returnTrain = FALSE)
# tmp<-rep(1:1498)
# for(k in 1:10){tmp[nfold[[k]]]<-k}
# blood_download<-blood_download[,rownames(sample_blood)]
# blood_download_72<-data.frame(t(blood_download[namescpgsnpblood,]))##blood_download_72<-data.frame(blood_download[,namescpgsnpblood])
#
# ##error in x[!which, , drop = FALSE] :
# ##(subscript) logical subscript too long => adjust the number of samples
#
# ######
#
E_model<-function(x_train_set, y_train_set,x_test_set, y_test_set,alpha){

  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x_train_set, y_train_set, alpha = alpha,foldid=tmp)

  #find optimal lambda value that minimizes test MSE
  best_lambda.se <- cv_model$lambda.1se
  best_lambda.min <- cv_model$lambda.min
  best_model.min <- glmnet(x_train_set, y_train_set,alpha = alpha, lambda = best_lambda.min)
  best_model.se <- glmnet(x_train_set, y_train_set,alpha = alpha, lambda = best_lambda.se)
  best_model.min.coefs<-as.matrix(coef(best_model.min))
  best_model.se.coefs<-as.matrix(coef(best_model.se))

  best_model.se.coefs<-names(best_model.se.coefs[best_model.se.coefs[,1]!=0,])
  best_model.min.coefs<-names(best_model.min.coefs[best_model.min.coefs[,1]!=0,])

  y_predicted.min <- predict(best_model.min, s = best_lambda.min, newx =  x_test_set)
  y_predicted.se <- predict(best_model.se, s = best_lambda.se, newx =  x_test_set)
  res.min<-data.frame(cbind("Bioage"=as.numeric(y_predicted.min),"Chroage"=as.numeric(y_test_set)))
  res.se<-data.frame(cbind("Bioage"=as.numeric(y_predicted.se),"Chroage"=as.numeric(y_test_set)))
  return(list(res.min,res.se,best_model.min,best_model.se,cv_model,best_model.min.coefs,best_model.se.coefs))
}


# #blood_download_72<-data.frame(t(blood_download_72)) ## CpGs should be in columns
#
# blood_download_72<-cbind("age"=as.numeric(sample_blood$age),blood_download_72)###1872x72219
# set.seed(23060830)
# train.names<-sample(rownames(blood_download_72), 1498)
# train_set<-blood_download_72[train.names,]
# test_set<-blood_download_72[!rownames(blood_download_72)%in%train.names,]
#
# y_train_set<-train_set[,colnames(train_set)=="age"]
# x_train_set<-as.matrix(train_set[,!colnames(train_set)=="age"])
#
#
# y_test_set<-test_set[,colnames(test_set)=="age"]
# x_test_set<-as.matrix(test_set[,!colnames(test_set)=="age"])
#
# res_72k<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1)
# png(filename = "blood_download_72k_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
# plot(res_72k[[5]])
# dev.off()
#
#
#
# # Measure: Mean-Squared Error
# #
# # Lambda Index Measure    SE Nonzero
# # min 0.1572   100   34.58 2.353     573
# # 1se 0.3631    82   36.91 1.899     212
# #
#
# ###### test the model with snp data
#
# imputed_methylation_snp_beta_small<-data.frame(imputed_methylation_snp_beta[,c("age",namescpgsnpblood)])
# y_test_set<-imputed_methylation_snp_beta_small[,colnames(imputed_methylation_snp_beta_small)=="age"]
# x_test_set<-as.matrix(imputed_methylation_snp_beta_small[,!colnames(imputed_methylation_snp_beta_small)=="age"])
#
# E_model_plot<-function(res, res_model,file_name, plt_name1,plt_name2){
#
#   y_predicted <- predict(res_model,newx =  x_test_set)
#   res_snp<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))
#
#
#   p1<- ggplot(res, aes(Chroage,Bioage))+
#     geom_point(size=0.01) +
#     geom_abline(linewidth = 0.2)+
#     xlim(0,100)+
#     ylim(0,100)+
#     #geom_smooth(method=lm, linewidth= 0.2,,se=FALSE,
#     # color="blue3")+
#     theme(text = element_text(size = 10, face = "bold"),legend.position = "none",legend.text=element_text(size=10)
#     )+
#     stat_regline_equation(label.y = 2,label.x =100,aes(label = plt_name1),hjust = 1,parse = TRUE)+
#
#     xlab("")+
#     ylab("")
#
#   p2<- ggplot(res_snp, aes(Chroage,Bioage))+
#     geom_point(size=0.01) +
#     geom_abline(linewidth = 0.2)+
#     xlim(0,100)+
#     ylim(0,100)+
#     # geom_smooth(method=lm, linewidth= 0.2,,se=FALSE,
#     #  color="blue3")+
#     theme(text = element_text(size = 10, face = "bold"),legend.position = "none",legend.text=element_text(size=10)
#     )+
#     stat_regline_equation(label.y = 2,label.x =100,aes(label = plt_name2),hjust = 1,parse = TRUE)+
#
#     xlab("")+
#     ylab("")
#   png(filename = file_name, width = 170, height = 80,units = "mm", res=300)
#   cowplot::plot_grid(p1,p2,ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
#   dev.off()
#
#   return(list(p1,p2,res,res_snp))
# }
#
# #res_72k.min<-  E_model_plot( res_72k[[1]], res_72k[[3]],"blood_download_72K_cpgs.snp.png","LASSO (EWAS)","LASSO (LA)")
# res_72k.se<-  E_model_plot( res_72k[[2]], res_72k[[4]],"blood_download_72K_cpgs1se.snp.png"
#                             ,paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\""))
#
# png(filename = "blood_download_72K_cpgs1se.snp.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_72k.se[[1]],res_72k.se[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
# E_model_summary<-function(res){
#
#   res$diff<-res$Bioage-res$Chroage
#
#   res_summary<-c("min"=min(res$diff),"median"=median(res$diff),"mean"=mean(res$diff),
#                  "max"=max(res$diff),"sd"=sd(res$diff),"rmse"=rmse(res$Bioage,res$Chroage),
#                  "mae"=mae(res$Bioage,res$Chroage),"mad"=mad(res$Bioage,res$Chroage))
#   return(res_summary)
# }
#
#
# res_72k.se_summary<-E_model_summary(res_72k.se[[3]])
#
# # min      median        mean         max          sd        rmse
# # -22.2943685  -0.4670545  -0.2701903  35.9178891   6.6333143   6.6299481
# # mae         mad
# # 4.5538265   4.9508761
#
#
# res_snp_72k.se_summary<-E_model_summary(res_72k.se[[4]])
#
# # min   median     mean      max       sd     rmse      mae      mad
# # -7.34220 13.18514 14.73918 53.82221 12.68188 19.44355 15.41007 19.54829
#
# ### cor>80 load("/mnt/qnap/researchData/cngueda/Impuls_Own/snp_data/blood_download_snp_cor80.RData")
# #blood_download_snp_cor80<-blood_download_snp_cor80[,!colnames(blood_download_snp_cor80)=="age"]
# #blood_download_snp_cor80<-data.frame(t(blood_download[colnames(blood_download_snp_cor80),]))
# blood_download_snp_cor80<-blood_download_snp_cor80[rownames(sample_blood),]
# blood_download_snp_cor80<-data.frame(cbind("age"=as.numeric(sample_blood$age),blood_download_snp_cor80))###1872x743
# blood_download_snp_cor80[]<-lapply(blood_download_snp_cor80,as.numeric)
# set.seed(23060830)
# train.names<-sample(rownames(blood_download_snp_cor80), 1498)
# train_set<-blood_download_snp_cor80[train.names,]
# test_set<-blood_download_snp_cor80[!rownames(blood_download_snp_cor80)%in%train.names,]
#
# y_train_set<-train_set[,colnames(train_set)=="age"]
# x_train_set<-as.matrix(train_set[,!colnames(train_set)=="age"])
#
#
# y_test_set<-test_set[,colnames(test_set)=="age"]
# x_test_set<-as.matrix(test_set[,!colnames(test_set)=="age"])
#
# res_742<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1)
# png(filename = "blood_download_742_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
# plot(res_742[[5]])
# dev.off()
#
# # Measure: Mean-Squared Error
# #
# # Lambda Index Measure    SE Nonzero
# # min 0.1513    45   114.5 4.163     306
# # 1se 0.2644    39   117.8 3.836     220
#
#
# ###### test the model with snp data
#
# imputed_methylation_snp_beta_small<-data.frame(imputed_methylation_snp_beta[,colnames(blood_download_snp_cor80)])##7473x743###imputed_methylation_snp_beta_small<-data.frame(t(imputed_methylation_snp_beta[colnames(blood_download_snp_cor80),]))
# #sum(rownames(imputed_methylation_snp_beta_small)%in%PV0697_datajoin$SIC)
# #imputed_methylation_snp_beta_small<-cbind("age"=PV0697_datajoin$ADULT_PROB_AGE,imputed_methylation_snp_beta_small)
# PV0697_datajoin<-PV0697_datajoin[rownames(imputed_methylation_snp_beta_small),]
# y_test_set<-imputed_methylation_snp_beta_small[,colnames(imputed_methylation_snp_beta_small)=="age"]
# x_test_set<-as.matrix(imputed_methylation_snp_beta_small[,!colnames(imputed_methylation_snp_beta_small)=="age"])
#
# res_742.se<-  E_model_plot( res_742[[2]], res_742[[4]],"blood_download_742cpgs1se.snp.png",
#                             paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\""))
#
#
# png(filename = "blood_download_742_cpgs1se.snp.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_742.se[[1]],res_742.se[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
#
# # E_model_summary(res_742.se[[3]])
# # min      median        mean         max          sd        rmse
# # -35.8002800   0.1464212  -0.2268346  43.0350893  11.0258444  11.0134303
# # mae         mad
# # 8.6152057  10.3512341
#
# # E_model_summary(res_742.se[[4]])
# # min     median       mean        max         sd       rmse        mae
# # -26.943177  -5.123898  -3.636307  35.229518  12.705073  13.214386  11.198450
# # mad
# # 15.468440
#
#
#
#
# ### With CpGs from Lutz: lutz pGs were selected with Generic proceedure. and the model was build with simple lineare regression
#
# #load("GA_CpG_subsets.RData") ### "CpGs_Blood_72k", "CpGs_Blood_742" ,"CpGs_SNP_179k"
#
# # for the common CpGs 72k was selected 116 CpGs
#
# set.seed(23060830)
# CpGs_Blood_72k_116<-CpGs_Blood_72k
# x<-rep(1:1498)
# nfold<-createFolds(1:1498, k = 10, list = TRUE, returnTrain = FALSE)
# tmp<-rep(1:1498)
# for(k in 1:10){tmp[nfold[[k]]]<-k}
#
# set.seed(23060830)
# train.names<-sample(rownames(CpGs_Blood_72k_116), 1498)
# train_set<-data.frame(CpGs_Blood_72k_116[train.names,])
# test_set<-data.frame(CpGs_Blood_72k_116[!rownames(CpGs_Blood_72k_116)%in%train.names,])
#
# y_train_set<-train_set[,colnames(train_set)=="age"]
# x_train_set<-train_set[,!colnames(train_set)=="age"]
# y_test_set<-test_set[,colnames(test_set)=="age"]
# x_test_set<-test_set[,!colnames(test_set)=="age"]
#
# res_CpGs_Blood_72k_116_model<-lm(age~.,data =train_set)
# length(summary(res_CpGs_Blood_72k_116_model)$coefficients[-1,1]!=0)
# y_predicted <- predict(res_CpGs_Blood_72k_116_model,  data.frame(x_test_set))
# res_CpGs_Blood_72k_116<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))
#
# ###test on snp data
#
# imputed_methylation_snp_beta_small<-data.frame(imputed_methylation_snp_beta[,colnames(CpGs_Blood_72k_116)])
# sum(rownames(imputed_methylation_snp_beta_small)%in%PV0697_datajoin$SIC)
# imputed_methylation_snp_beta_small<-imputed_methylation_snp_beta_small[rownames(PV0697_datajoin),]
# imputed_methylation_snp_beta_small[]<-lapply(imputed_methylation_snp_beta_small,as.numeric)
#
# y_test_set<-imputed_methylation_snp_beta_small[,colnames(imputed_methylation_snp_beta_small)=="age"]
# x_test_set<-imputed_methylation_snp_beta_small[,!colnames(imputed_methylation_snp_beta_small)=="age"]
# y_predicted <- predict(res_CpGs_Blood_72k_116_model,  data.frame(x_test_set))
# res_snp_CpGs_Blood_72k_116<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))
#
#
# ##plots
# E_model_plotlm<-function(res, res_snp,file_name, plt_name1,plt_name2){
#
#   p1<- ggplot(res, aes(Chroage,Bioage))+
#     geom_point(size=0.01) +
#     geom_abline(linewidth = 0.2)+
#     # geom_smooth(method=lm, linewidth= 0.2,,se=FALSE,
#     #  color="blue3")+
#     theme(text = element_text(size = 10, face = "bold"),legend.position = "none",legend.text=element_text(size=10)
#     )+
#     stat_regline_equation(label.y = 2,label.x =100,aes(label = plt_name1),hjust = 1,parse = TRUE)+
#     xlim(0,100)+
#     ylim(0,100)+
#     xlab("")+
#     ylab("")
#
#   p2<- ggplot(res_snp, aes(Chroage,Bioage))+
#     geom_point(size=0.01) +
#     geom_abline(linewidth = 0.2)+
#     # geom_smooth(method=lm, linewidth= 0.2,,se=FALSE,
#     # color="blue3")+
#     theme(text = element_text(size = 10, face = "bold"),legend.position = "none",legend.text=element_text(size=10)
#     )+
#     stat_regline_equation(label.y = 2,label.x =100,aes(label = plt_name2),hjust = 1,parse = TRUE)+
#     xlim(0,100)+
#     ylim(0,100)+
#     xlab("")+
#     ylab("")
#   png(filename = file_name, width = 170, height = 90,units = "mm", res=300)
#   cowplot::plot_grid(p1,p2,ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
#   dev.off()
#
#   return(list(p1,p2,res,res_snp))
# }
#
# res_72k_116.plt<-  E_model_plotlm( res_CpGs_Blood_72k_116, res_snp_CpGs_Blood_72k_116,"blood_download_72K_cpgs.snp.png",
#                                    paste( "GA","*\"  (EWAS)\""),paste( "GA","*\"  (LA)\""))
#
# png(filename = "res_blood_snp_CpGs_Blood_72k_116.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_72k_116.plt[[1]],res_72k_116.plt[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
#
# res_CpGs_Blood_72k_116_summary<-E_model_summary(res_CpGs_Blood_72k_116)
# #
# # min       median         mean          max           sd         rmse
# # -14.95575499  -0.36989420  -0.05209278  27.09965559   4.60670613   4.60083823
# # mae          mad
# # 3.36440522   4.04375973
#
#
# res_snp_CpGs_Blood_72k_116_summary<-E_model_summary(res_snp_CpGs_Blood_72k_116)
#
#
# # min    median      mean       max        sd      rmse       mae       mad
# # -33.46893 -13.00501 -11.42999  27.69799  12.68477  17.07415  14.25480  21.31923
#
# # for the common CpGs 742 was selected 53 CpGs
#
# CpGs_Blood_742_53<-blood_download_72[,colnames(CpGs_Blood_742)]
# sum(rownames(imputed_methylation_snp_beta_small)%in%PV0697_datajoin$SIC)
# imputed_methylation_snp_beta_small<-imputed_methylation_snp_beta_small[rownames(PV0697_datajoin),]
# set.seed(23060830)
# train.names<-sample(rownames(CpGs_Blood_742_53), 1498)
# train_set<-CpGs_Blood_742_53[train.names,]
# test_set<-CpGs_Blood_742_53[!rownames(CpGs_Blood_742_53)%in%train.names,]
#
# y_train_set<-train_set[,colnames(train_set)=="age"]
# x_train_set<-as.matrix(train_set[,!colnames(train_set)=="age"])
# y_test_set<-test_set[,colnames(test_set)=="age"]
# x_test_set<-as.matrix(test_set[,!colnames(test_set)=="age"])
#
# res_CpGs_Blood_742_53_model<-lm(age~.,data =train_set)
# length(summary(res_CpGs_Blood_742_53_model)$coefficients[-1,1]!=0)
# y_predicted <- predict(res_CpGs_Blood_742_53_model,  data.frame(x_test_set))
# res_CpGs_Blood_742_53<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))
#
# ###test on snp data
#
# imputed_methylation_snp_beta_small<-data.frame(imputed_methylation_snp_beta[,colnames(CpGs_Blood_742)])
#
#
# imputed_methylation_snp_beta_small[]<-lapply(imputed_methylation_snp_beta_small,as.numeric)
# y_test_set<-imputed_methylation_snp_beta_small[,colnames(imputed_methylation_snp_beta_small)=="age"]
# x_test_set<-imputed_methylation_snp_beta_small[,!colnames(imputed_methylation_snp_beta_small)=="age"]
# y_predicted <- predict(res_CpGs_Blood_742_53_model,  data.frame(x_test_set))
# res_snp_CpGs_Blood_742_53<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))
#
#
# ##plots
#
# res_742_53.plt<-  E_model_plotlm( res_CpGs_Blood_742_53, res_snp_CpGs_Blood_742_53,"blood_download_72K_cpgs.snp.png",
#                                   paste( "GA","*\"  (EWAS)\""),paste( "GA","*\"  (LA)\""))
#
#
# png(filename = "res_blood_snp_CpGs_Blood_742_53CpGs.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_742_53.plt[[1]],res_742_53.plt[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
#
# res_CpGs_Blood_742_53_summary<-E_model_summary(res_CpGs_Blood_742_53)
#
# # min       median         mean          max           sd         rmse
# # -30.31097794  -0.01432100  -0.06597923  39.62553896  10.94876241  10.93431430
# # mae          mad
# # 8.55785532  10.56764828
#
#
# res_snp_CpGs_Blood_742_53_summary<-E_model_summary(res_snp_CpGs_Blood_742_53)
#
#
# # min    median      mean       max        sd      rmse       mae       mad
# # -39.73892 -17.70745 -16.13586  23.09645  12.71594  20.54360  17.57941  26.73403
# #
#
# ### with all snp CpG
#
# #imputed_methylation_snp_beta
# set.seed(23060830)
# x<-rep(1:5978)
# nfold<-createFolds(1:5978, k = 10, list = TRUE, returnTrain = FALSE)
# tmp<-rep(1:5978)
# for(k in 1:10){tmp[nfold[[k]]]<-k}
#
# set.seed(23060830)
# train.names<-sample(rownames(imputed_methylation_snp_beta), 5978)
# train_set<-imputed_methylation_snp_beta[train.names,]
# test_set<-imputed_methylation_snp_beta[!rownames(imputed_methylation_snp_beta)%in%train.names,]
#
#
# y_train_set<-train_set[,colnames(train_set)=="age"]
# x_train_set<-as.matrix(train_set[,!colnames(train_set)=="age"])
#
# dim(x_train_set)
# y_test_set<-test_set[,colnames(test_set)=="age"]
# x_test_set<-as.matrix(test_set[,!colnames(test_set)=="age"])
#
# #
# res_snp_all<-E_model(x_train_set, y_train_set,x_test_set, y_test_set,1)
#
# png(filename = "snp_all_modelplotlambdacpgs.snp.png", width = 170, height = 90,units = "mm", res=300)
# plot(res_snp_all[[5]])
# dev.off()
#
# # Call:  cv.glmnet(x = x_train_set, y = y_train_set, foldid = tmp, alpha = alpha)
# #
# # Measure: Mean-Squared Error
# #
# # Lambda Index Measure    SE Nonzero
# # min 0.7246     1   162.4 2.632       0
# # 1se 0.7246     1   162.4 2.632       0
#
# #ici
#
# res_snp_all.se<-  E_model_plot( res_snp_all[[2]], res_snp_all[[4]],"snp_all_modelplotlambdacpgs1se.snp.png",
#                                 paste( "LASSO","*\"  (EWAS)\""),paste( "LASSO","*\"  (LA)\""))
#
# png(filename = "res_snp_all_cpgs1se.snp.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_snp_all.se[[1]],res_snp_all.se[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
#
#
#
#
# res_snp_all.se_summary<-E_model_summary(res_snp_all.se[[4]])
#
# # min      median        mean         max          sd        rmse
# # -21.8917548  -1.4817548  -0.2998885  37.6282452  12.4480411  12.4474903
# # mae         mad
# # 10.5792763  15.4957716
#
#
#
# res_72k_116.plt<-  E_model_plotlm( res_CpGs_Blood_72k_116, res_snp_CpGs_Blood_72k_116,"blood_download_72K_cpgs.snp.png",
#                                    paste( "GA","*\"  (EWAS)\""),paste( "GA","*\"  (LA)\""))
#
#
# png(filename = "res_blood_snp_CpGs_Blood_72k_116.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_72k_116.plt[[1]],res_72k_116.plt[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
# # for the all snp CpGs 3 was selected  CpGs
# imputed_methylation_snp_beta_small<-data.frame(imputed_methylation_snp_beta[,CpGs_Snp_snp_3])###CpGs_Snp_snp_3<-c("age","cg14795672","cg10731507","cg11994425")
#
# set.seed(23060830)
# train.names<-sample(rownames(imputed_methylation_snp_beta_small), 5978)
# train_set<-imputed_methylation_snp_beta_small[train.names,]
# test_set<-imputed_methylation_snp_beta_small[!rownames(imputed_methylation_snp_beta_small)%in%train.names,]
#
# y_train_set<-train_set[,colnames(train_set)=="age"]
# x_train_set<-as.matrix(train_set[,!colnames(train_set)=="age"])
# y_test_set<-test_set[,colnames(test_set)=="age"]
# x_test_set<-as.matrix(test_set[,!colnames(test_set)=="age"])
#
# res_CpGs_snp_3_model<-lm(age~.,data =train_set)
# length(summary(res_CpGs_snp_3_model)$coefficients[-1,1]!=0)
# y_predicted <- predict(res_CpGs_snp_3_model,  data.frame(x_test_set))
# rmse(as.numeric( y_predicted),as.numeric(y_test_set))
# res_CpGs_snp_3<-data.frame(cbind("Bioage"=as.numeric( y_predicted),"Chroage"=as.numeric(y_test_set)))
#
#
#
# res_snp_all.3<-  E_model_plotlm( res_CpGs_snp_3, res_CpGs_snp_3,"snp_all_modelplotlambdacpgs1se.snp.png",
#                                  paste( "GA","*\"  (EWAS)\""),paste( "GA","*\"  (LA)\""))
#
#
# png(filename = "res_CpGs_snp_3.png", width = 170, height = 90,units = "mm", res=300)
# cowplot::plot_grid(res_snp_all.3[[1]],res_snp_all.3[[2]],ncol =2,nrow=1,label_colour ="Blue3", label_x = '0', label_y = '1')
# dev.off()
#
#
# res_CpGs_snp_3_summary<- E_model_summary(res_CpGs_snp_3)
#
# # sd        rmse
# # -23.0955440  -1.5497850  -0.3447418  38.0972109  12.3841641  12.3848206
# # mae         mad
# # 10.5262697  15.3840901
#
#
# ##7473 179095
# png(filename = "blood_download_72K_cpgs.snplasso_116GA2.png", width = 170, height = 150,units = "mm", res=300)
# cowplot::plot_grid(res_72k.se[[1]]+theme(#axis.text.y = element_blank(),
#   #axis.ticks.y = element_blank(),
#   axis.title.y = element_blank(),
#   #axis.text.x = element_blank(),
#   #axis.ticks.x = element_blank(),
#   axis.title.x = element_blank() ),
#   res_72k.se[[2]]+theme(
#     axis.title.y = element_blank(),
#
#     axis.title.x = element_blank() ),
#   res_72k_116.plt[[1]]+theme(
#     axis.title.y = element_blank(),
#
#     axis.title.x = element_blank() ),
#   res_72k_116.plt[[2]]+theme(
#     axis.title.y = element_blank(),
#
#     axis.title.x = element_blank() ),
#   labels = "AUTO",vjust=3.3,hjust=-5.5,scale=0.90,label_size = 10
# )+
#   draw_label("Chronological age (years)", x=0.5, y=  0, vjust= -0.3,size = 12,fontface = "bold", angle= 0)+
#   draw_label("Biological age (years)", x=  0, y=0.5, vjust= 1.1,size = 12,fontface = "bold", angle=90)
# dev.off()
#
#
#
#
#
#
# png(filename = "blood_download_snp_cor80_742scpgs.snplasso_53GA.png", width = 170, height = 150,units = "mm", res=300)
# cowplot::plot_grid(res_742.se[[1]]+theme(#axis.text.y = element_blank(),
#   #axis.ticks.y = element_blank(),
#   axis.title.y = element_blank(),
#   #axis.text.x = element_blank(),
#   #axis.ticks.x = element_blank(),
#   axis.title.x = element_blank() ),
#   res_742.se[[2]]+theme(
#     axis.title.y = element_blank(),
#
#     axis.title.x = element_blank() ),
#   res_742_53.plt[[1]]+theme(
#     axis.title.y = element_blank(),
#
#     axis.title.x = element_blank() ),
#   res_742_53.plt[[2]]+theme(
#     axis.title.y = element_blank(),
#
#     axis.title.x = element_blank() ),
#   labels = "AUTO",vjust=3.3,hjust=-5.5,scale=0.90,label_size = 10)+
#   draw_label("Chronological age (years)", x=0.5, y=  0, vjust= -0.3,size = 12,fontface = "bold", angle= 0)+
#   draw_label("Biological age (years)", x=  0, y=0.5, vjust= 1.1,size = 12,fontface = "bold", angle=90)
# dev.off()
#
#
#
# png(filename = "snp_allcpgs.snp_snp_3.png", width = 170, height = 150,units = "mm", res=300)
# cowplot::plot_grid(res_snp_all.se[[2]]+theme(#axis.text.y = element_blank(),
#   #axis.ticks.y = element_blank(),
#   axis.title.y = element_blank(),
#   #axis.text.x = element_blank(),
#   #axis.ticks.x = element_blank(),
#   axis.title.x = element_blank() ),
#   res_snp_all.3[[2]]+theme(
#     axis.title.y = element_blank(),
#     axis.title.x = element_blank() ),
#   res_snp_all.se[[2]]+theme(
#     axis.title.y = element_blank(),
#     axis.title.x = element_blank() ),
#   res_snp_all.3[[2]]+theme(
#     axis.title.y = element_blank(),
#     axis.title.x = element_blank() ),
#   labels = "AUTO",vjust=3.3,hjust=-5.5,scale=0.90,label_size = 10)+
#   draw_label("Chronological age (years)", x=0.5, y=  0, vjust= -0.3,size = 12,fontface = "bold", angle= 0)+
#   draw_label("Biological age (years)", x=  0, y=0.5, vjust= 1.1,size = 12,fontface = "bold", angle=90)
# dev.off()
#
# #####convert M to beta values
#
# x<-x[,-1]
# x[]<-lapply(x, as.numeric)
#
# x_b<-2^x/(1+2^x)
#
# #BiocManager::install("lumi")
# library(lumi)
#
# x_beta<-as.matrix(x)
# x_beta<-m2beta(x_beta)
#
#
# #### use of methylimp to impute missing horvath CpG
# ##we use annotation illumina reference file 850 K and 72K CpGs for chrom and Strand information
# load("../data450_850.RData")##same for 27.Rdata
#
# imputed_methylation_snp_beta_hknn<-imputed_methylation_snp_beta[,-1]
#
# imputed_methylation_snp_beta_hknn[,snp_horvath_missing_CpGs]<- as.numeric(NA)
#
# data27<-data27[,-c(2,4)]###same for 850
# data27<-data.frame(data27)###same for 850
# x<-colnames(imputed_methylation_snp_beta_hknn)[!colnames(imputed_methylation_snp_beta_hknn)%in%rownames(data850)]
# data850<-rbind(data27[x,],data850)
# data850<-data850[colnames(imputed_methylation_snp_beta_hknn),]
# data850<-t(data850)
# #data8500<-data850[!(data850$chr %in% "chr8" & data850$strand %in% "+"),]
#
#
# imputed_methylation_snp_beta_hknn<-data.frame(rbind(data850,imputed_methylation_snp_beta_hknn))
#
# imputed_methylation_snp_beta_hknn <- split(imputed_methylation_snp_beta_hknn, f=list(imputed_methylation_snp_beta_hknn[,1],
#                                                                                      imputed_methylation_snp_beta_hknn[,2]))
#
# for(i in 1:length(imputed_methylation_snp_beta_hknn)){
#   if(nrow(imputed_methylation_snp_beta_hknn[[i]])==0){
#     print(names(imputed_methylation_snp_beta_hknn)[i])
#     v<-c(v,names(imputed_methylation_snp_beta_hknn)[i])
#   }
# }
# imputed_methylation_snp_beta_hknn[[v]]<-NULL
#
#
# for(i in 1:length(imputed_methylation_snp_beta_hknn)){
#
#   imputed_methylation_snp_beta_hknn[[i]]<-imputed_methylation_snp_beta_hknn[[i]][,-1:-2]
#   imputed_methylation_snp_beta_hknn[[i]][]<-lapply(imputed_methylation_snp_beta_hknn[[i]],as.numeric)
#   imputed_methylation_snp_beta_hknn[[i]]<-t(imputed_methylation_snp_beta_hknn[[i]])
# }
#
# for(i in 1:length(imputed_methylation_snp_beta_hknn)){
#   imputed_methylation_snp_beta_hknn[[i]]<-methyLImp2_internal(imputed_methylation_snp_beta_hknn[[i]], min=0, max=1)
# }
#
#
# for(i in 1:length(imputed_methylation_snp_beta_hknn_c)){
#   print(sum(is.na(imputed_methylation_snp_beta_hknn_c[[i]])))
# }
#
#
#




# library(AnnotationHub)
# library(minfi)
# library(Biobase)
# #library(tibble)
# library(ggpmisc)
# library(GEOquery)
# library(Metrics)
# library(impute)
# library(ggplot2)
# library(dplyr)
# library(ggpubr)
# library(glmnet)
# library(caret)
# library(readxl)
# library(cowplot)

