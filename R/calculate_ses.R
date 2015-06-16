#' Make up SES (household income) relationship with final height for oversampled data
#'
#' Using the oversampling of "boys", the returned object from "sample_data()";
#' imputes a made-up linear relationship between fabricated SES variable (think household income)
#' and the final height in column 32 of data_in.  In future I will build in
#' params for the 100 and the 12 and the rnorm in teh SES assignment.
#' @param data_in an object returned from sample_data()
#' @param slope a number that determines the linear relationship of final height and SES.
#' For instance, if slope=100 and intercept=12, then the SES for individual i is 
#' SES(i) = 100*( height(i) - minimum of all heights )/(max of all heights - min of all heights) + 12 + rnorm(1)
#' @param intercept a number, see explanation of slope
#' @keywords oversampling
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' with_ses <- calculate_ses(over_samp_mat)
#' 
calculate_ses <- function(data_in=NULL, slope=100, intercept=12){
  ## @knitr SESHeight
  ## induce linear relationship between SES and Height,
  ## that ranges between 12 and 112.  We can suppose that this SES is household income
  ## in 1000 USD or something like that.  Then add to overSampMat (now it is 33rd column)
  nc <- ncol(data_in)
  SES <- slope*(data_in[,nc]-min(data_in[,nc]))/
               (max(data_in[,nc])-min(data_in[,nc])) + intercept +
               rnorm(dim(data_in)[1])
  data_in$ses <- SES
  data_in
}