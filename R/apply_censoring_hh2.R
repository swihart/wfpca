#' Apply censoring scheme that is related to SES
#'
#' Takes output from make_long() and applies censoring
#' by predicting/fitting from a simple model 
#' 
#'     
#' Names are hard coded, so stick
#' to script in the examples.
#' @param data_in an object returned from calculate_ses()
#' @param protected a vector corresponding to age / colnames of age_vec that are protected from censoring so that fpca.face() wont bonk
#' @export
#' @return prob.cens as a column added to data_in representing the probability of being censored at that time.
#' @return instudy.sim as a column added to data_in representing whether that particular observation is instudy (1) or not (0)
#' @return instudy as a column added to data_in builds upon instudy.sim in that once instudy.sim is 0 instudy is 0 for all remaining observation times
#' @examples
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' with_ses <- calculate_ses(over_samp_mat)
#' long <-make_long(with_ses)
#' head(long)
#' censored <- apply_censoring(long)
#' head(censored,18)
apply_censoring_hh2 <- function(data_in=NULL, simple_model=simplog, protected=c(1,1.25,1.5,1.75,2.0,3)){  
  #library(plyr)
  data_in$prob.cens <- predict.glm(simplog, newdata=data_in, type="response")
  #data_in[,prob.cens:=predict.glm(simplog, newdata=data_in, type="response")]
  
  ##data_in$prob.cens[data_in$age==1] <- 0 ## everyone has at least baseline value, set prob.cens to 0
  ## BJS EDIT:  for fpca.face to work, we need at least four values. So,
  ## I'm updating accordingly.  The first four corresponds to 1,1.25,1.5,1.75 years ofage.
  #data_in$prob.cens[data_in$age %in% protected] <- 0
  
  ## errors...
#   data_in[,
#           prob.cens=c(0*protected, prob.cens[-1*protected]),
#           by="newid"]  
  ## still errors so i put these lines as `whamcast` in the insample_sim_hh2.R script
  data_in[, subject_time_index:={ c(1:.N)}, by="newid"]
  data_in[subject_time_index %in% protected, prob.cens:=0]

  
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
#stay <- ddply(as.data.frame(data_in), .(newid), function(w) cumprod(w$instudy.sim))

data_in[, instudy:=cumprod(instudy.sim), by="newid"]

#stay.melt <- melt(stay, id.vars="newid")
#stay.melt <- stay.melt[order(stay.melt$newid,stay.melt$variable),]
#summary(stay.melt)
#data_in$instudy <- stay.melt[,3]
#qplot(x=age,y=instudy,colour=ses, data=data_in)
data_in
}