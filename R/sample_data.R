#' Sample the prepped data of the boys of the Berkeley Growth Study Data
#'
#' This returns an oversampling of "boys", the returned object from "prep_data()".
#' @param data_in an object returned from prep_data()
#' @param tot_subj total number of subjects sampled with replacement -- numeric, multiple of 100.
#' @keywords oversampling
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' overSampMat<-sample_data(d,1000)
sample_data <- function(data_in=NULL, tot_subj=1000){
  ## @knitr sampleBoys
  ## we set a seed for reproducibility and randomly sample with replacement
  ## nSampId subjects from 'boys' using a row number indicator, not ids.
  ## we then sort on the 32nd column, which is final recorded height
  set.seed(101)
  nSampId <- tot_subj # needs to be multiple of 100
  overSampRow <- sample(nrow(data_in), nSampId, replace=T)
  overSampMat <- data_in[overSampRow,]
  overSampMat <- overSampMat[order(overSampMat[,32],decreasing=TRUE),]
}