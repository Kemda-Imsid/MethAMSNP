#' find_missing_CpGs
#' Find the missing CpGs names for a given data set for a given clock.
#'
#' @param data_miss a dataframe containing the CpGs in columns and Sample in rows. the dataframe in wich the missing CpGs be imputed.
#' @param clock a clock names amount these, "coefBLUP", "coefEN", "coefHannum", "coefHorvath", "coefLevine", "coefPedBE", "coefSkin", "coefTL", "coefWu")
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{missing}{CpGs required for the clock that are not found in the input data frame.}
#'   \item{present}{CpGs required for the specified clock that are present in the input data frame.}
#' }
#' @export

find_missing_CpGs<-function(data_miss,clock){
  data_miss<-data.frame(t(data_miss[1:3,]))
  data_miss<-data.frame(cbind("ProbeID"=rownames(data_miss),data_miss))
  data_miss[,-1]<-lapply(data_miss[,-1],as.numeric)
  missing_CpGs <- methylclock::checkClocks(data_miss,clocks=clock,localHub=TRUE)
  if(length(missing_CpGs[[clock]]!=0)){
    missing_CpGs<-missing_CpGs[[clock]]
  }else{
    missing_CpGs<-missing_CpGs[[clock]]
  }
  coefs<-get(paste0("coef",clock))
  present_CpGs<-coefs$CpGmarker[-1][!coefs$CpGmarker[-1]%in%missing_CpGs]
  return(list(missing_CpGs,present_CpGs))
}
#load("Data/imputed_methylation_snp_beta_s.RData")
#missing_CpGs<-find_missing_CpGs(imputed_methylation_snp_beta_s,"Horvath")

#' add_missing_CpGs
#' the function crates some empty value columns to the dataframe for the missings CpGs
#' @param data a dataframe containing the CpGs in rows and Sample in columns. the dataframe in wich the missing CpGs be imputed.
#' @param x missing CpGs names
#'
#' @returns a data frame containing empty columns for the missing CpGs names
#' @export

add_missing_CpGs<-function(data, x){
    # Create a list of empty columns (NA values)
    x <- data.frame(setNames(
      replicate(length(x), rep(NA, nrow(data)), simplify = FALSE), x))
    # Bind the empty columns to the original data frame
    data <- cbind(x,data)
    return(data)
  }

#data<-add_missing_CpGs(imputed_methylation_snp_beta_s,missing_CpGs[[1]])



#' Title
#'
#' @param data a dataframe containing the CpGs in rows and Sample in columns. the dataframe in wich the missing CpGs be imputed.
#' @param x missing CpGs names
#'
#' @returns a data frame containing CpGs and they nearest neighbours
#' @export
#'
find_nearest_CpGs<- function(data, x) {
  # Ensure required columns exist
  if (!all(c("chr", "pos", "strand", "Name") %in% colnames(data))) {
    stop("Data must have columns: chr, pos, strand, Name")
  }
  data$pos<-as.numeric(data$pos)
  # Remove rows with NA positions
  data <- data[!is.na(data$pos), ]
  data<-data[order(data$pos),]
  # Split data into available and missing CpGs
  available <- data[!data$Name %in% x, ]
  missing   <- data[data$Name %in% x, ]

  # Output structure
  result <- data.frame(
    chr = character(),
    pos = numeric(),
    strand = character(),
    Name = character(),
    cpgstoreplace = character(),
    dist = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(missing))) {
    m <- missing[i, ]

    # Only compare to same chr and strand
    subset <- available[available$chr == m$chr & available$strand == m$strand, ]
    if (nrow(subset) == 0) next  # skip if no match

    # Find closest by genomic position
    distances <- abs(as.numeric(subset$pos) - as.numeric(m$pos))
    closest_cpg <- subset[which.min(distances),]

    result <- rbind(result, data.frame(
      chr = m$chr,
      pos = m$pos,
      strand = m$strand,
      Name = m$Name,
      cpgstoreplace = closest_cpg$Name,
      dist = distances[closest_idx],
      stringsAsFactors = FALSE
    ))
  }

  return(result)
}


#' replace_missings_CpGs
#' replace missing CpGs values whithin single samples with their nearest neighbor CpGs values availble in the dataframe
#' @param data a dataframe in which the missing CpGs values for a given clock should be replaced, containing the CpGs in rows and Sample in columns.
#' @param x a vector containing the missing CpGs to replace the missing ones for a given clock
#' @returns a data frame containing CpGs
#' @export


replace_missings_CpGs<-function(data,x){

  ref450<- read.csv("Data/ref450.csv")
  ref850<- read.csv("Data/ref850.csv")
  df_annotation<-data.frame(rbind(cbind("chr"=ref850$chr,"pos"=ref850$pos,"strand"=ref850$strand,"Name"=ref850$Name),
                                  cbind("chr"=ref450$chr,"pos"=ref450$pos,"strand"=ref450$strand,"Name"=ref450$Name)))##a little bit more than the 850
  df_annotation<-df_annotation[!duplicated(df_annotation$Name),]
  data_split<-data.frame(cbind("chr" =  df_annotation$chr[match(colnames(data), df_annotation$Name)],
                         "pos" = as.numeric(df_annotation$pos[match(colnames(data), df_annotation$Name)]),
                         "strand" =  df_annotation$strand[match(colnames(data), df_annotation$Name)],"Name"=colnames(data)))
  # data<-data.frame(data[order(data[,"chr"],data[,"strand"],data[,"pos"]),])
  data_split <- split(data_split, f=list(data_split$chr, data_split$strand))
  v<-c()
  for(i in 1:length(data_split)){
    if(nrow(data_split[[i]])==0){
      v<-c(v,names(data_split)[i])
    }
  }
  if(length(v)!=0){
    data_split[[v]]<-NULL
  }
  # find nearest neighbors for each element of the list
  nearest_CpGs <- mclapply(data_split, function(item) {
    find_nearest_CpGs(item, x)
  }, mc.cores = 4)
  # Combine results into a single data frame
  nearest_CpGs <- do.call(rbind, nearest_CpGs)
  for(i in 1:nrow(nearest_CpGs)){
    data[,nearest_CpGs[i,"Name"]]<-data[,nearest_CpGs[i,"cpgstoreplace"]]
  }
  return(list( data,nearest_CpGs))
}

data<-replace_missings_CpGs(data,x)


