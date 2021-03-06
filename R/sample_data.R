#' Sample the prepped data of the boys of the Berkeley Growth Study Data
#'
#' This returns an oversampling of "boys", the returned object from "prep_data()".
#' @param data_in an object returned from prep_data()
#' @param tot_subj total number of subjects sampled with replacement -- numeric, multiple of 100.
#' @param seed an integer for set.seed(seed) to assist reproducible simulations
#' @keywords oversampling
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' overSampMat<-sample_data(d,1000)
sample_data <- function(data_in=NULL, tot_subj=1000, seed=101, timepoints=as.numeric(names(data_in[,-1]))){
  ## @knitr sampleBoys
  ## we set a seed for reproducibility and randomly sample with replacement
  ## nSampId subjects from 'boys' using a row number indicator, not ids.
  ## we then sort on the 32nd column, which is final recorded height
  10*seed
  ##set.seed(seed)
  nSampId <- tot_subj # needs to be multiple of 100
  overSampRow <- sample(nrow(data_in), nSampId, replace=T)
  overSampMat <- data_in[overSampRow,]
  ## with perturbation...
  overSampMat.pert <- overSampMat[,-1] + runif(nrow(overSampMat),-2,2)
  overSampMat <- cbind(overSampMat$id, overSampMat.pert)
  colnames(overSampMat)[1] <- "id"
  overSampMat <- overSampMat[order(overSampMat[,32],decreasing=TRUE),]
}