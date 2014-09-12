#' Sample the prepped data of the boys of the Berkeley Growth Study Data
#'
#' Using the oversampling of "boys", the returned object from "sample_data()";
#' imputes a made-up linear relationship between fabricated SES variable (think household income)
#' and the final height in column 32 of data_in.  In future I will build in
#' params for the 100 and the 12 and the rnorm in teh SES assignment.
#' @param data_in an object returned from sample_data()
#' @keywords oversampling
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' long_ses <- calculate_ses(over_samp_mat)
#' 
calculate_ses <- function(data_in=NULL){
  ## @knitr SESHeight
  ## induce linear relationship between SES and Height,
  ## that ranges between 12 and 112.  We can suppose that this SES is household income
  ## in 1000 USD or something like that.  Then add to overSampMat (now it is 33rd column)
  SES <- 100*(data_in[,32]-min(data_in[,32]))/
             (max(data_in[,32])-min(data_in[,32])) + 12 +
             rnorm(dim(data_in)[1])
  data_in$ses <- SES
}