#' Sample the prepped data of the boys of the Berkeley Growth Study Data
#' based on the fpc as opposed to the oversampling in "sample_data()"
#'
#' This returns an oversampling of "boys", the returned object from "prep_data()".
#' @param data_in an object returned from prep_data()
#' @param tot_subj total number of subjects sampled with replacement -- numeric, multiple of 100.
#' @param seed an integer for set.seed(seed) to assist reproducible simulations
#' @param timepoints the grid that the columns of data_in represent
#' @keywords oversampling
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' overSampMat<-sample_data(d,1000)
sample_data_fpc <- function(data_in=NULL, tot_subj=1000, seed=101, timepoints=age_vec){
  ff<-fpca.face(Y=data_in[,-1], argvals=timepoints, knots=(length(timepoints)-2)) 
  ## empirical covariance matrix
  covmat = cov(ff$scores)
  ## explicitly put eigenvalues as variances on the diagonal
  diag(covmat) <- ff$eigenvalues
  ## independent:
  ##covmat <- diag(ff$eigenvalues)
  library(mvtnorm)
  sim_scores<-rmvnorm(tot_subj,
                      rep(0,length(ff$eigenvalues)),
                      covmat) 
  
  sim_curves<-ff$eigenvectors%*%t(sim_scores) + ff$mu
  ## with some meas error (upward bias) -- to see if lme fails less often
  sim_curves<-sim_curves +  
    matrix(runif(prod(dim(sim_curves)), -.125, 3), nrow=nrow(sim_curves))
  ##range(sim_curves[31,]) ## range seems smaller than what's observed...
  
  sim_data <- data.frame(cbind(1:tot_subj, t(sim_curves)))
  colnames(sim_data) <- colnames(data_in)
  sim_data <- sim_data[order(sim_data[,32],decreasing=TRUE),]
  sim_data$id <- 1:nrow(sim_data)
  sim_data
}