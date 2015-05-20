#' Calculate STIPW based on propensity for being study
#'
#' Names are hard coded, so stick
#' to script in the examples.
#' @param data_in an object returned from apply_censoring()
#' @param na.action "keep" means the dataset returned will be same number of rows as data_in and "omit" discards all data past censoring observation
#' @export
#' @return a dataset similar to data_in with stipw in columns; and potentially fewer rows or NA-filled rows
#' for induced censoring.  
#' @examples
#' \dontrun{
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' with_ses <- calculate_ses(over_samp_mat)
#' long <-make_long(with_ses)
#' head(long)
#' censored <- apply_censoring(long)
#' head(censored,18)
#' observed_with_stipw <- calculate_stipw(censored,"keep")
#' }
calculate_stipw <- function(data_in=NULL, na.action="omit"){

  ######################################################
  # BL edits
  ######################################################
  # only do the propensity among observations upto dropout.
  # This is because the pooled logistic is a time to event propensity score.
  # Therefore keeping in censored observations
  #  is like the individual is censored over and over again
  ######################################################
  # BJS notes on BL edits
  ######################################################
  # indcen.age is just switching 1's and 0's of instudy.sim
  # see this with a
  # table(data_in$indcen.age,data_in$instudy.sim,useNA="ifany")
  # cumsum.ind, and cumsum.ind2 are shown in the following two code-lines:
  # tail(data_in[,c("newid","indcen.age", "cumsum.ind", "cumsum.ind2")])
  # (data_in[data_in$newid==1000,c("newid","indcen.age", "cumsum.ind", "cumsum.ind2")])
  # cumsum.ind  adds up 1's and 0's, which are independent
  # cumsum.ind2 adds up cumsum.ind.  By taking data_in.sub to be <=1, we
  # only keep the data up to the first indcen.age==1.  cumsum.ind2 is necessary
  # b/c we see multiple rows of cumsum.ind can equal 1.
  data_in$indcen.age[data_in$instudy.sim==0] <- c(1)
  data_in$indcen.age[is.na(data_in$indcen.age)] <- 0
  data_in$cumsum.ind <- ave(data_in$indcen.age,data_in$newid,FUN=cumsum)
  data_in$cumsum.ind2 <- ave(data_in$cumsum.ind,data_in$newid,FUN=cumsum)
  
  data_in.sub <- data_in[data_in$cumsum.ind2<=1,]
  
  
  bl.logit.nocov <- glm(instudy~bs(age,intercept=TRUE,df=7)-1    , family="binomial", data=data_in.sub)
  bl.logit.ses   <- glm(instudy~bs(age,intercept=TRUE,df=7)-1+ses, family="binomial", data=data_in.sub)
  
  #qplot(bl.logit.nocov$fitted);range(bl.logit.nocov$fitted)
  #qplot(bl.logit.ses$fitted);range(bl.logit.ses$fitted)
  
  ## add the fitted probabilities to the dataframe
  #data_in$fitted.nocov <- (logit.nocov$fitted)
  #data_in$fitted.ses   <- (logit.ses$fitted)
  data_in.sub$bl.fitted.nocov <- (bl.logit.nocov$fitted)
  data_in.sub$bl.fitted.ses   <- (bl.logit.ses$fitted)
  
  
  ## now sub again so that we only have observed for weight calculation...
  data_in.sub2 <- data_in.sub[data_in.sub$cumsum.ind2==0,]
  
  data_in.sub2$cumprod.nocov <- ave(data_in.sub2$bl.fitted.nocov,data_in.sub2$newid,FUN=cumprod)
  data_in.sub2$cumprod.ses <- ave(data_in.sub2$bl.fitted.ses,data_in.sub2$newid,FUN=cumprod)
  data_in.sub2$stipw <- data_in.sub2$cumprod.nocov/data_in.sub2$cumprod.ses
  
  ## instead of setting this to 0, once weights are calculated for everyone else,
  ## remove the censored line for each individual...
  ## data_in.sub$stipw[!data_in.sub$instudy] <- 0
  
  
  
  
  #qplot(cumprod.nocov,cumprod.ses,color=age,size=ses,data=data_in.sub)
  #qplot(stipw,cumprod.ses,color=age,size=ses,data=data_in.sub)
  #qplot(stipw,cumprod.nocov,color=age,size=ses,data=data_in.sub)
  
  ##sub.wtd_avg <- ddply(data_in.sub, .(age), function(w) sum(w$stipw*w$inches)/sum(w$stipw))
  ##names(sub.wtd_avg)[2] <- "wtd_avg"
  
  
 if(na.action=="omit"){
   return_obj=data_in.sub2
 }
 if(na.action=="keep"){
   back_in <- merge(data_in, data_in.sub2,all.x=TRUE,sort=FALSE)
   return_obj=back_in
 } 
  
  
  return_obj
}