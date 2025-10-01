
#' DNA_age help in prepossessing of the data and estimates DNAm age giving a DNA methylation clock.
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

#' DNA_age_metrics provides some metrics for the estimation of the age imputation
#'
#' @param dat data frame containing in columns Biological and chronological age
#'
#' @returns a table with the metrics
#' @export

DNA_age_metrics<-function(dat){

  dat$diff<-dat$Bioage-dat$Chroage
  dat<- data.frame("min_diff"=min(dat$diff),"median_diff"=median(dat$diff),"mean_diff"=mean(dat$diff),
                  "max_diff"=max(dat$diff), "sd_diff"=sd(dat$diff), "rmse"=rmse(dat$Chroage,dat$Bioage),
                   "mae"=mae(dat$Chroage,dat$Bioage),"mad"=mad(dat$Chroage,dat$Bioage))
 return(dat)
}


#' DNA_age_plot provides a plot based on the computed data
#'
#' @param dat_list list of data frames containing in column Chronological age and estiomated Biological age
#' @param target_col String column name of the Chronological age
#'
#' @returns a list of resulted plots
#' @export

DNA_age_plot <- function(dat_list, target_col) {
  plot_list <- list()

  # Loop through the data frame list
  for (dat in 1:length(dat_list)){
    lab<<-names(dat_list[dat])
    # Check if target column exists
    if (!(target_col %in% colnames(dat_list[[dat]]))) {
      stop("Target column not found in data frame.")
    }

    # Only proceed if target is numeric
    if (!is.numeric(dat_list[[dat]][,target_col])) {
      stop("Target column must be numeric for scatter plots.")
    }
  # Loop through other columns
  # for (col_name in colnames(dat_list[[dat]])) {
  #   # Skip the target column itself
  #   if (col_name == target_col) next

    col_name <-colnames(dat_list[[dat]])[!colnames(dat_list[[dat]]) == target_col]
    print(col_name)
    # Only plot if the other column is numeric
    if (is.numeric(dat_list[[dat]][,col_name])) {
      print(lab)
      p <-  ggplot(dat_list[[dat]], aes(x = .data[[target_col]], y = .data[[col_name]],label = lab)) +
        geom_point(size=0.01) +
        geom_abline(linewidth = 0.2)+
        #geom_smooth(method=lm, linewidth= 0.2,,se=FALSE,
        #  color="blue3")+
        theme(text = element_text(size = 10, face = "bold"),legend.position = "none",legend.text=element_text(size=10)
        )+
        stat_regline_equation(label.y = 2,label.x =125,hjust = 1,parse = TRUE)+
        xlim(0,125)+
        ylim(0,125)+
        xlab("Chronological age (years)")+
        ylab("Biological age (years)")
      # Add to list
      plot_list[[lab]] <- p
    }
  #
  # }


  }

  return(plot_list)
}



######

#' E_model uses elastic net to build models
#'
#' @param x_train_set a matrix, train data set
#' @param y_train_set a numeric vector containing the chronological age  for the train set
#' @param x_test_set a matrix, test data set
#' @param y_test_set a numeric vector containing the chronological age  for the test set
#' @param alpha a number between 0 and 1 in elastic net. 0 corrsponds to ridge and to 1 to lasso lasso regression
#' @param foldid a vector of integers indicating which fold each observation belongs to
#' @return A list containing the following elements:
#' \describe{
#'   \item{res.min}{Results (e.g., predictions or metrics) using lambda.min}
#'   \item{res.se}{Results using lambda.1se}
#'   \item{best_model.min}{Model fitted with lambda.min (the lambda with minimum CV error)}
#'   \item{best_model.se}{Model fitted with lambda.1se (more regularized model)}
#'   \item{cv_model}{The full cross-validation glmnet model object}
#'   \item{best_model.min.coefs}{Coefficients from the model at lambda.min}
#'   \item{best_model.se.coefs}{Coefficients from the model at lambda.1se}
#' }
#' @export

E_model<-function(x_train_set, y_train_set,x_test_set, y_test_set,alpha,foldid){

  #perform k-fold cross-validation to find optimal lambda value

  cv_model <- cv.glmnet(x_train_set, y_train_set, alpha = alpha,foldid=foldid)

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





