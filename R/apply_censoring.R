#' Apply censoring scheme that is related to SES
#'
#' Takes output from make_long() and applies censoring
#' by creating a made-up, linear, hard-coded relationship between the probability
#' of being censored at each observation time and SES  
#' Names are hard coded, so stick
#' to script in the examples.
#' @param data_in an object returned from calculate_ses()
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' with_ses <- calculate_ses(over_samp_mat)
#' long <-make_long(with_ses)
#' head(long)
#' censored <- apply_censoring(long)
#' head(censored)
apply_censoring <- function(data_in=NULL){
  data_in$prob.cens <- .1*(1-(data_in$ses-min(data_in$ses))/(max(data_in$ses)-min(data_in$ses))) + 
                    .2*((data_in$age-min(data_in$age))/(max(data_in$age)-min(data_in$age)))
  #####################
  
  ##data_in$prob.cens[data_in$age==1] <- 0 ## everyone has at least baseline value, set prob.cens to 0
  ## BJS EDIT:  for fpca.face to work, we need at least four values. So,
  ## I'm updating accordingly.  The first four corresponds to 1,1.25,1.5,1.75 years ofage.
  data_in$prob.cens[data_in$age %in% c(1,1.25,1.5,1.75,2.0,3)] <- 0
  ## everyone has at least 4 baseline values
  ## set prob.cens to 0
  #plot(data_in$age,data_in$prob.cens)
  #plot(data_in$ses,data_in$prob.cens)
  #qplot(x=ses,y=prob.cens,colour=age, data=data_in)
  
  ## @knitr instudy
  ## We randomly draw drop-out based bernoulli variables based on prob.cens,
  ## independent for every measurement.  This means each subject has a vector of 0s and 1s
  ## across age.  To enforece "once dropped out stay dropped out", we then take a cumulative
  ## product ascending through age where once a 0 occurs the remaining elements will be 0.
  ## (a 0 means dropped out, 1 means instudy)
  data_in$instudy.sim <- 1-rbinom(nrow(data_in),1,data_in$prob.cens)
  stay <- ddply(data_in, .(newid), function(w) cumprod(w$instudy.sim))
  stay.melt <- melt(stay, id.vars="newid")
  stay.melt <- stay.melt[order(stay.melt$newid,stay.melt$variable),]
  summary(stay.melt)
  data_in$instudy <- stay.melt[,3]
  #qplot(x=age,y=instudy,colour=ses, data=data_in)
  data_in
}