LinRegGA<-function(trainingSet,age_col_name,maxIter=3000,nPop=1000){

require(glmnet)
require(GA)
avg_CpG_start_num <-200 ### internal constant
# to determine the age column
p <- ncol(trainingSet) - 1
n <- nrow(trainingSet)
ageCol <- which(colnames(trainingSet) == age_col_name)

# the fitness function for GA
decodingFunc <- function(chromosome, trainingSet, nSamples, ageCol, cpgNames) {

  # to init
  nPredictors <- sum(chromosome)
  if ((nPredictors > 0) & (nPredictors < nSamples)) {
    # to conduct linear regressions
    predNames <- cpgNames[chromosome == 1]
    dfNames <- c(age_col_name, predNames)
    linModel <- lm(as.formula(paste(age_col_name, "~ .")), data=trainingSet[, dfNames])
    return(-BIC(linModel))
  }
  else {
    return(-1e8)
  }
}

# selection function with exponential ranking for GA
selectionFunc <- function(object) {

  ranks <- rank(object@fitness, ties.method="random")

  # high selection pressure
  # tau <- log(0.1) / ((1 - object@popSize)) # P(selection of the worst individual) / P(selection of the best individual) = 1 / 10
  # moderate selection pressure
  tau <- log(0.2) / ((1 - object@popSize)) # P(selection of the worst individual) / P(selection of the best individual) = 1 / 5
  # low selection pressure
  # tau <- log(0.5) / ((1 - object@popSize)) # P(selection of the worst individual) / P(selection of the best individual) = 1 / 2
  weights <- exp(tau * (ranks - object@popSize))
  weights <- weights / sum(weights)

  sel <- sample(1:object@popSize, size=object@popSize, prob=weights, replace=TRUE)
  out <- list(population=object@population[sel, , drop=FALSE], fitness=object@fitness[sel])

  return(out)
}

# balanced flip mutation function
balancedFlipMutationFunc <- function(object, parent) {

  mutRate <- 1 / object@nBits
  mutate <- as.vector(object@population[parent,])

  # single flip mutation
  if (runif(1) < 0.5) {
    idx <- which(mutate == 0)
  } else {
    idx <- which(mutate == 1)
  }
  if(length(idx) > 0) {
    idx <- sample(idx, 1)
    mutate[idx] <- 1 - mutate[idx]
  }

  # balanced flip mutation
  nmbOfFlips <- rbinom(1, object@nBits, mutRate)
  if (nmbOfFlips > 0) {
    zeroPos <- which(mutate == 0)
    onePos <- which(mutate == 1)
    numPairs <- min(nmbOfFlips, length(zeroPos), length(onePos))

    if (numPairs > 0) {
      zeroInd <- sample(zeroPos, numPairs, replace=FALSE)
      oneInd <- sample(onePos, numPairs, replace=FALSE)
      mutate[zeroInd] <- 1
      mutate[oneInd] <- 0
    }

    # add on for SNP data
    if (sum(mutate) < 2) {
      mutate[sample(length(mutate), 4, replace=FALSE)] <- 1
    }
  }

  return(mutate)
}

# Genetic Algorithm
set.seed(as.numeric(Sys.time()))
nBits <- ncol(trainingSet) - 1
cpgNames <- colnames(trainingSet[,-ageCol])
GA <- ga(type = "binary",
    crossover = GA:::gabin_uCrossover,
    mutation = balancedFlipMutationFunc,
    fitness = decodingFunc,
    selection = selectionFunc,
    nBits = nBits,
    popSize = nPop,
    maxiter = maxIter,
    pmutation = 1, # chromosome based
    pcrossover = 0.90,
    elitism = 1,
    suggestions = 1 * matrix(runif(nPop * nBits), nrow=nPop) < avg_CpG_start_num / nBits,
    trainingSet = trainingSet,
    nSamples = nrow(trainingSet),
    ageCol = ageCol,
    cpgNames =cpgNames)

 solution<-cpgNames[GA@solution==1]
 return(solution)
}
