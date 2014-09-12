#' Make up SES (household income) relationship with final height for oversampled data
#'
#' Takes wide format with individual covariates and makes it long while
#' redoing the ids (original, oversampled ids retained as 'id', 'new.id' is unique) for a sensical melt.  Names are hard coded, so stick
#' to script in the examples.
#' @param data_in an object returned from calculate_ses()
#' @keywords oversampling
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' with_ses <- calculate_ses(over_samp_mat)
#' long <-make_long(with_ses)
#' head(long)
make_long <- function(data_in=NULL){
  ## @knitr moreDataPrep
  ## To get overSampMat to where it is now, we've done some sorting and some sampling but will
  ## need newid for sensical melting; so make one that is akin to row number.
  ## Then melt() the data to a long, ggplot friendly format, from the matrix, one subject per row format in 'boys'
  ## After melt()ing we rename some variables and convert the age (which was in the columns of 'boys') to a numeric (otherwise it is a factor)
  data_in$newid <- 1:nrow(data_in)
  long <- melt(data=data_in,
               id.vars=c("newid","id","ses"), variable.name="age",value.name="inches")
  long <- long[order(long$newid, long$age),]
  names(long)[names(long)=="value"] <- "inches"
  long$age <- (as.numeric(levels(long$age))[long$age] )
  long
}