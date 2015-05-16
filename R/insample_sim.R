#' Run a simulation and then perform prediction.  Based on readMe_sim_prediction(), which was never used..
#'
#' This function comes after many months of running readme_sim and talks with Bryan Lau on 
#' what we need to demonstrate for the method.
#'
#' @param sim_seed passed to set.seed()
#' @param sample_size (defaults to 1000) and must be a multiple of 100
#' @param sim_slope see slope in calculate_ses()
#' @param sim_intercept see intercept in calculate_ses()
#' @param sim_ses_coef see ses_coef in apply_censoring()
#' @param sim_age_coef see age_coef in apply_censoring()
#' @export
#' @return results a data frame with rmse for each approach for that simulated dataset
#' @examples
#' ---
insample_sim <- function(sim_seed=101, sample_size=1000, sim_slope=100,
                                  sim_intercept=12, sim_ses_coef=.01, sim_age_coef=.01){
  ## quick start with the default values as variables for testing.  Loads packages.
  #test_prep_script()
  ## get the data prepped, that is simulated and censored and calculate missing. 
  #simulate_censor_summarize()

  ## just for testing; comment out before github commits
  sim_seed=101; sample_size=1000; sim_slope=100; sim_intercept=12; sim_ses_coef=.01; sim_age_coef=.01;
  
  
  ## d will be a 39 x 32 matrix based on the males from the Berkeley Growth Study in `fda` package
  ## we oversample d based on the fpc as well extract out the ages of measurement
  d<-prep_data()
  age_vec <- c(as.numeric(colnames(d[-1])))
  over_samp_mat<-sample_data_fpc(d, sample_size, seed=sim_seed, timepoints=age_vec)
  
  ## we calculate ses on the oversampled dataset and turn the matrix into a long dataset.
  with_ses <- calculate_ses(over_samp_mat, slope=sim_slope, intercept=sim_intercept)
  long <-make_long(with_ses)

  ## In this chunk we apply censoring.
  censored <- apply_censoring(long, ses_coef=sim_ses_coef, age_coef=sim_age_coef)
  
  ## Measurement error: within 1/8 inch.  Can comment out or change.
  censored$inches <- censored$inches + runif(length(censored$inches), -1/8,1/8)
  
  ## observed_with_stipw has the standardized inverse probability weights (stipw)
  ## wtd_trajectories has same info as observed_with_stipw but has inches_wtd as well
  ## we calculate all_with_stipw with NAs
       all_with_stipw <- calculate_stipw(censored,"keep")
  observed_with_stipw <- calculate_stipw(censored,"omit")
  wtd_trajectories    <- calculate_wtd_trajectories(observed_with_stipw)
  
  
  
  
  ## use data.table where possible speeds up future rbinds post simulation
#setDT(long)
#setkey(long, id, age)

setDT(all_with_stipw)
setkey(all_with_stipw, newid, age, inches, ses)

setDT(wtd_trajectories)
setkey(wtd_trajectories, newid, age, inches, ses)

interim <- wtd_trajectories[all_with_stipw]
all_with_stipw<-subset(interim,
                       select=c(names(all_with_stipw),
                                "inches_wtd",
                                "remaining",
                                "stipw",
                                "stipw01",
                                "stipw01.n"))



all_with_stipw[            , inches_ltfu:=inches ]
all_with_stipw[is.na(stipw), inches_ltfu:=NA ]

## standardize the names selected across all approaches and their resultant  data sets
selected<-c("newid","age", "ses", 
            
            "inches", "inches_wtd", "inches_ltfu", "inches_predicted",
            
            "remaining",
            "stipw",
            "stipw01",
            "stipw01.n",
            
            "approach")

  ## Note: we can calulate this post-simulation
  ## truth, all data:
  ##true_avg <- ddply(long, .(age), function(w)  mean(w$inches))
  ##true_avg$approach<- "true_avg"

## note: we're simplifying to (w)fpca and (w)lme
## compare the bias of mean trajectories of the following approaches applied to the
## censored/observed data, not the truth:
##a) naive-nonparametric: that is, just na.rm=TRUE and turn the mean() crank
# naive_non_parm_avg <- ddply(observed_with_stipw, .(age), function(w)  mean(w$inches))
# naive_non_parm_avg$approach <- "naive_non_parm_avg"
# ## combine with long for the prediction:
# setDT(naive_non_parm_avg)
# setnames(naive_non_parm_avg, "V1", "inches_predicted")
# setkey(naive_non_parm_avg, age, inches_predicted)
# setkey(long, age) ## no id in this method.
# dim(naive_non_parm_avg)
# dim(long)
# naive_non_parm_avg <- naive_non_parm_avg[long]
# setkey(naive_non_parm_avg, newid, age, inches_predicted)
# naive_non_parm_avg <- subset(naive_non_parm_avg, select=c("newid","age", "inches", "inches_predicted", "approach"))


## note: we're simplifying to (w)fpca and (w)lme
##b) weighted-nonparametric: a little more sophisticated, na.rm=TRUE and turn the weighted mean() crank
# wtd_non_parm_avg <- ddply(observed_with_stipw, .(age), function(w) sum(w$stipw*w$inches)/sum(w$stipw))
# wtd_non_parm_avg$approach <- "wtd_non_parm_avg"
# ## combine with long for the prediction:
# setDT(wtd_non_parm_avg)
# setnames(wtd_non_parm_avg, "V1", "inches_predicted")
# setkey(wtd_non_parm_avg, age, inches_predicted)
# setkey(long, age) ## no id in this method.
# dim(wtd_non_parm_avg)
# dim(long)
# wtd_non_parm_avg <- wtd_non_parm_avg[long]
# setkey(wtd_non_parm_avg, newid, age, inches_predicted)
# wtd_non_parm_avg <- subset(wtd_non_parm_avg, select=c("newid","age", "inches", "inches_predicted", "approach"))


##c)  naive-FPC,
unwtd_fncs      <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches")
fpca_unwtd_fncs <- fpca.face(Y=as.matrix(unwtd_fncs)[,-1], argvals=age_vec, knots=26)
naive_fpc       <- data.frame(age=age_vec, V1=fpca_unwtd_fncs$mu, approach="naive_fpc")
## combine with long for the prediction:
## a little different than previous examples; now we do have individual level curves
## need to extract them (Yhat) and rename them and data.table them
naive_fpc_indiv <- as.data.frame(cbind(1:nrow(fpca_unwtd_fncs$Yhat),fpca_unwtd_fncs$Yhat))
colnames(naive_fpc_indiv) <- colnames(unwtd_fncs)
setDT(naive_fpc_indiv)
naive_fpc <- melt(naive_fpc_indiv,      
                  id.vars=c("newid"),  
                  variable.name = "age", 
                  variable.factor=FALSE,
                  value.name="inches_predicted")
naive_fpc[,approach:="naive_fpc",]
naive_fpc[,age:=as.numeric(age)]
setDT(naive_fpc)
setkey(naive_fpc, newid, age)
naive_fpc  <- naive_fpc[all_with_stipw]
setkey(naive_fpc, newid, age)
naive_fpc <- subset(naive_fpc, select=selected)

## visual checks
# ggplot(data=naive_fpc[newid %in% c(1,2,5)], aes(x=age, y=inches_predicted, color=factor(newid)))+geom_point()+
#   geom_point(aes(y=inches), color='black')






## note: we're simplifying to (w)fpca and (w)lme
##d) naive-FPC-post-adjusted-by-weights, and
# pabw <- dcast(data=wtd_trajectories, formula= newid~age, value.var="stipw01.n")
# pabw[pabw==0]<-NA
# pabw_fpca_unwtd_fncs <- fpca_unwtd_fncs$Yhat*pabw[,-1] ## element-wise
# naive_fpc_pabw_avg <- colMeans(pabw_fpca_unwtd_fncs,na.rm=TRUE)
# naive_fpc_pabw <- data.frame(age=age_vec, V1=naive_fpc_pabw_avg, approach="naive_fpc_pabw")
# ## combine with long for the prediction:
# ## a little different than previous examples; now we do have individual level curves
# ## need to extract them (Yhat) and rename them and data.table them
# naive_fpc_pabw_indiv <- as.data.frame(cbind(1:nrow(pabw_fpca_unwtd_fncs),pabw_fpca_unwtd_fncs))
# colnames(naive_fpc_pabw_indiv) <- colnames(unwtd_fncs)
# setDT(naive_fpc_pabw_indiv)
# naive_fpc_pabw <- melt(naive_fpc_pabw_indiv,      
#                   id.vars=c("newid"),  
#                   variable.name = "age", 
#                   variable.factor=FALSE,
#                   value.name="inches_predicted")
# naive_fpc_pabw[,approach:="naive_fpc_pabw",]
# naive_fpc_pabw[,age:=as.numeric(age)]
# setDT(naive_fpc_pabw)
# setkey(naive_fpc_pabw, newid, age, inches_predicted)
# naive_fpc_pabw <- naive_fpc_pabw[all_with_stipw]
# setkey(naive_fpc_pabw, newid, age, inches_predicted)
# naive_fpc_pabw <- subset(naive_fpc_pabw, select=c("newid","age", "inches", "inches_ltfu", "inches_predicted", "approach"))


##e) weighted-FPC.
wtd_fncs <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches_wtd")
fpca_wtd_fncs <- fpca.face(Y=as.matrix(wtd_fncs)[,-1], argvals=age_vec, knots=26)
#weighted_fpc <- data.frame(age=age_vec, V1=fpca_wtd_fncs$mu, approach="weighted_fpc")
## combine with long for the prediction:
## a little different than previous examples; now we do have individual level curves
## need to extract them (Yhat) and rename them and data.table them
weighted_fpc_indiv <- as.data.frame(cbind(1:nrow(fpca_wtd_fncs$Yhat),fpca_wtd_fncs$Yhat))
colnames(weighted_fpc_indiv) <- colnames(wtd_fncs)
setDT(weighted_fpc_indiv)
weighted_fpc <- melt(weighted_fpc_indiv,      
                  id.vars=c("newid"),  
                  variable.name = "age", 
                  variable.factor=FALSE,
                  value.name="inches_predicted")
weighted_fpc[,approach:="weighted_fpc",]
weighted_fpc[,age:=as.numeric(age)]
setDT(weighted_fpc)
setkey(weighted_fpc, newid, age)
weighted_fpc <- weighted_fpc[all_with_stipw]
## begin extra steps: (weighted inputs average out for mean, but need to be de-weighted for individual)
setDT(wtd_trajectories)
wts <- subset(wtd_trajectories, select=c("newid","age","stipw01.n"))
setkey(wts, newid, age)

## can't do curve completion with these three knuckleheads
#weighted_fpc.test<-weighted_fpc[wts]
#weighted_fpc.test[, inches_predicted_weighted:= inches_predicted]
#weighted_fpc.test[, inches_predicted_deweighted:= inches_predicted/stipw01.n]

weighted_fpc.test<-wts[weighted_fpc]
## see how many 0's
weighted_fpc.test[,table(round(stipw01.n,2))]
## change 0's to NA
weighted_fpc.test[stipw01.n==0, stipw01.n:=NA]
weighted_fpc.test[, inches_predicted_weighted:= inches_predicted]
weighted_fpc.test[!is.na(stipw01.n), inches_predicted_deweighted:= inches_predicted/stipw01.n]
## get the wtd_population_mean in there:
weighted_fpc.test[,wtd_pop_mean:=fpca_wtd_fncs$mu,by=newid]
##for now, if don't have observed data there, we just imputed the weighted mean
## kinda lame, think on it.
weighted_fpc.test[is.na(stipw01.n), inches_predicted_deweighted:= wtd_pop_mean, by=newid]
## end extra steps:
setkey(weighted_fpc.test, newid, age)#, inches_predicted_deweighted)
setnames(weighted_fpc.test, "inches_predicted", "inches_predicted_old")
setnames(weighted_fpc.test, "inches_predicted_deweighted", "inches_predicted")
weighted_fpc <- subset(weighted_fpc.test, select=selected)


# ggplot(data=weighted_fpc[newid %in% c(1,2,5)], aes(x=age, y=inches_predicted, color=factor(newid)))+geom_point()+
#   geom_point(aes(y=inches_ltfu), color='black')



##f) lme
naive_lme<-tryCatch(
{
  #naive_lme_model<-lme(inches ~ bs(age, df=15), random=~1|newid, data=observed_with_stipw);
  #naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), level=0), approach="naive_lme")
  
  naive_lme_model<-lmer(inches ~ bs(age, df=15) + (1|newid), data=observed_with_stipw);
  ## re.form=~0
  #naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="naive_lme")
  ## re.form=~1
  naive_lme <- data.table(newid            = all_with_stipw$newid,
                          age              = all_with_stipw$age,
                          ses              = all_with_stipw$ses,
                          inches           = all_with_stipw$inches,
                          inches_wtd       = all_with_stipw$inches_wtd,
                          inches_ltfu      = all_with_stipw$inches_ltfu,
                          inches_predicted = predict(naive_lme_model, 
                                                     newdata=all_with_stipw),
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="naive_lme")
},
warning =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  ##naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  ;
  naive_lme <- data.table(newid            = all_with_stipw$newid,
                          age              = all_with_stipw$age,
                          ses              = all_with_stipw$ses,
                          inches           = all_with_stipw$inches,
                          inches_wtd       = all_with_stipw$inches_wtd,
                          inches_ltfu      = all_with_stipw$inches_ltfu,
                          inches_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="naive_lme")
},
error =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  #naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  ;
  naive_lme <- data.table(newid            = all_with_stipw$newid,
                          age              = all_with_stipw$age,
                          ses              = all_with_stipw$ses,
                          inches           = all_with_stipw$inches,
                          inches_wtd       = all_with_stipw$inches_wtd,
                          inches_ltfu      = all_with_stipw$inches_ltfu,
                          inches_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="naive_lme")
})
## quick checks:
# summary(naive_lme)
# ggplot(data=naive_lme[newid %in% c(2)], 
#       aes(x=age, y=inches_predicted, color=factor(newid)))+
#   geom_point()+
#   geom_point(aes(y=inches), color='black')



## I have two wtd_lme chunks -- only have one uncommented at a time!
## the one immediately preceding this comment is lme(of wtd inches);
## whereas the other one below it are lme4::lmer of inches.
# wtd_lme<-tryCatch(
# {
#   wtd_lme_model<-lme(inches_wtd ~ bs(age, df=15), random=~1|newid, data=wtd_trajectories,
#                      na.action=na.omit);
#   ##predict(wtd_lme_model, newdata=data.frame(age=age_vec), level=0)
#   wtd_lme <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), level=0), approach="wtd_lme")
# },
# warning =function(cond){
#   wtd_lme <- data.frame(age=age_vec, V1=NA, approach="wtd_lme")  ;
#   wtd_lme
# },
# error =function(cond){
#   wtd_lme <- data.frame(age=age_vec, V1=NA, approach="wtd_lme")  ;
#   wtd_lme
# })
# summary(wtd_lme)
#

wtd_lme<-tryCatch(
{
  #   wtd_lme_model<-lme(inches_wtd ~ bs(age, df=15), random=~1|newid, data=wtd_trajectories,
  #                      na.action=na.omit);
  wtd_lme_model<-lmer(inches_wtd ~ bs(age, df=15) + (1|newid), data=wtd_trajectories,
                      na.action=na.omit);
  ##predict(wtd_lme_model, newdata=data.frame(age=age_vec), level=0)
  ##wtd_lme <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), level=0), approach="wtd_lme")
  ## note: below, use re.form=~0 in lme4:predict is equivalent to level=0 in nlme:lme
  ## re.form=~0
  ##wtd_lme <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="wtd_lme")
  wtd_lme_pop_mean <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="wtd_lme")
  
  ## re.form=~1
  wtd_lme1  <- data.table(newid            = all_with_stipw$newid,
                          age              = all_with_stipw$age,
                          ses              = all_with_stipw$ses,
                          inches           = all_with_stipw$inches,
                          inches_wtd       = all_with_stipw$inches_wtd,
                          inches_ltfu      = all_with_stipw$inches_ltfu,
                          inches_predicted = predict(wtd_lme_model, 
                                                     newdata=all_with_stipw),
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_lme")
  
  wtd_lme.test<- wts[wtd_lme1]
  ## see how many 0's
  wtd_lme.test[,table(round(stipw01.n,2))]
  ## change 0's to NA
  wtd_lme.test[stipw01.n==0, stipw01.n:=NA]
  wtd_lme.test[, inches_predicted_weighted:= inches_predicted]
  wtd_lme.test[!is.na(stipw01.n), inches_predicted_deweighted:= inches_predicted/stipw01.n]
  ## get the wtd_population_mean in there:
  wtd_lme.test[,wtd_pop_mean:=wtd_lme_pop_mean$V1,by=newid]
  ##for now, if don't have observed data there, we just imputed the weighted mean
  ## kinda lame, think on it.
  wtd_lme.test[is.na(stipw01.n), inches_predicted_deweighted:= wtd_pop_mean, by=newid]
  ## end extra steps:
  setkey(wtd_lme.test, newid, age)#, inches_predicted_deweighted)
  setnames(wtd_lme.test, "inches_predicted", "inches_predicted_old")
  setnames(wtd_lme.test, "inches_predicted_deweighted", "inches_predicted")
  wtd_lme <- subset(wtd_lme.test, select=selected)  
  
  },
warning =function(cond){
  write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  #wtd_lme <- data.frame(age=age_vec, V1=cond, approach="wtd_lme")  ;
  wtd_lme   <- data.table(newid            = all_with_stipw$newid,
                          age              = all_with_stipw$age,
                          ses              = all_with_stipw$ses,
                          inches           = all_with_stipw$inches,
                          inches_wtd       = all_with_stipw$inches_wtd,
                          inches_ltfu      = all_with_stipw$inches_ltfu,
                          inches_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_lme")
},
error =function(cond){
  write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  #wtd_lme <- data.frame(age=age_vec, V1=cond, approach="wtd_lme")  ;
  wtd_lme   <- data.table(newid            = all_with_stipw$newid,
                          age              = all_with_stipw$age,
                          ses              = all_with_stipw$ses,
                          inches           = all_with_stipw$inches,
                          inches_wtd       = all_with_stipw$inches_wtd,
                          inches_ltfu      = all_with_stipw$inches_ltfu,
                          inches_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_lme")
})
## quick checks:
# summary(wtd_lme)
# ggplot(data=wtd_lme[newid %in% c(1,2,5)], 
#       aes(x=age, y=inches_predicted, color=factor(newid)))+
#   geom_point()+
#   geom_point(aes(y=inches_ltfu), color='black')
# ## see the original predictions, before 'deweighting'
# ggplot(data=wtd_lme.test[newid %in% c(1,2,5)], 
#        aes(x=age, y=inches_predicted_old, color=factor(newid)))+
#   geom_point()+
#   geom_point(aes(y=inches_ltfu), color='black')





#rbind it!!!!!!!!!!!!!



## plot each approach's mean
## means <- rbind(true_avg, naive_non_parm_avg, wtd_non_parm_avg, naive_fpc, naive_fpc_pabw, weighted_fpc, naive_lme, wtd_lme)
# note: we're simplifying to (w)fpca and (w)lme, so we rbind() fewer than line above (if you decided to add more approaches later,
#         make sure they have same format):
means <- rbind(naive_fpc, weighted_fpc, naive_lme, wtd_lme)

# from older files:
# means$approach<-factor(means$approach, levels=unique(means$approach))
# colnames(means)[colnames(means)=="V1"]<- "inches"


## next three chunks deal with missingness, then we rbind it us with means....

## this is a useful chunk to check the range of 
## probality of being censored (prob.cens)
#   ddply(censored, .(age), function(w) sum(w$instudy==1) )
#   melt.prob.cens=ddply(censored, .(newid,age), function(w) w$prob.cens )
#   dcast.prob.cens=dcast(melt.prob.cens, newid~age, value.var="V1")
#   apply(dcast.prob.cens, 2, function(w) round(range(w),2))
#   head(censored,18)

## do the following to get precent missing at age 18:
## overall:
dcast.wtd.trajectories<-dcast(calculate_wtd_trajectories(calculate_stipw(censored,"keep")), newid~age, value.var="stipw")
percent.missing.at.age.18=sum(is.na(dcast.wtd.trajectories["18"]))/length(unlist(dcast.wtd.trajectories["18"]))
percent.missing = colSums(is.na(dcast.wtd.trajectories[,-1]))/length(unlist(dcast.wtd.trajectories["18"]))

## below/above median SES:
dcast.wtd.trajectories<-dcast(calculate_wtd_trajectories(calculate_stipw(censored,"keep")), newid+ses~age, value.var="stipw")
medianSES<-median(dcast.wtd.trajectories$ses, na.rm=TRUE)
subbie<-subset(dcast.wtd.trajectories, ses <= medianSES,"18")
percent.missing.at.age.18.below.median=sum(is.na(subbie))/nrow(subbie)
subbie<-subset(dcast.wtd.trajectories, ses <= medianSES, select=c(-1,-2))
percent.missing.below.median=colSums(is.na(subbie))/nrow(subbie)
subbie<-subset(dcast.wtd.trajectories, ses  > medianSES,"18")
percent.missing.at.age.18.above.median=sum(is.na(subbie))/nrow(subbie)
subbie<-subset(dcast.wtd.trajectories, ses  > medianSES, select=c(-1,-2))
percent.missing.above.median=colSums(is.na(subbie))/nrow(subbie)

## note: this cbind() works because every age is present for each dataset in `means`
##       and the only non-scalars are vectors that are same length as number of age-levels
## return:
results <- cbind(
        sim_seed = sim_seed,
        sample_size = sample_size,
        sim_slope = sim_slope,
        sim_intercept = sim_intercept,
        sim_ses_coef = sim_ses_coef,
        sim_age_coef = sim_age_coef,
        perc_ltfu_18      =percent.missing.at.age.18,
        percent_missing = percent.missing,
        percent_missing_below_median = percent.missing.below.median,
        percent_missing_above_median = percent.missing.above.median,
        means)




## started tinkering with Out Sample predictions...please see `outsample_sim.R`
# 
# ##a) naive-nonparametric: that is, just na.rm=TRUE and turn the mean() crank
# ## the best we can do to predict for those missing is to predict the mean curve for everyone
# 
# 
# 
# 
# 
# ## DO SOME PREDICTIONS
# ## put the following as function inputs... then remove!
# oos_sample_size=100
# oos_age=5
# oos_ind=as.numeric(names(d[,-1])) <= oos_age
# oos_timepoints=as.numeric(names(d[,-1]))[oos_ind]
# oos<-sample_data_fpc(d, oos_sample_size, seed=sim_seed, timepoints=as.numeric(names(d[,-1])))
# oos.early <- oos 
# oos.early[,!c(TRUE, oos_ind)] <- NA
# ##oos.early <- rbind(c(999999,fpca_wtd_fncs$mu), oos.early)
# 
#   oos.early <- oos[,c(TRUE, oos_ind)]
# 
# 
# oos.late <- oos[,c(TRUE, !oos_ind)]
# 
# ##e) weighted-FPC.
# ##wtd_fncs <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches_wtd")
# fpca_wtd_fncs_pred_oos <- fpca.face(Y=as.matrix(wtd_fncs[,-1]), Y.pred=as.matrix(oos.early[,-1]), argvals=age_vec, knots=10)
# 
# weighted_fpc_pred_oos <- data.frame(age=age_vec, V1=fpca_wtd_fncs_pred_oos$mu, approach="weighted_fpc_pred_oos")
# 
# fpca_wtd_fncs_pred_oos$Yhat
# fpca_wtd_fncs_pred_oos$scores
# 
# oos<-sample_data_fpc(d[,c(TRUE,oos_ind)], oos_sample_size, seed=sim_seed, timepoints=oos_timepoints, knots.face=4)
# 
# 
# 








##  currently not using; going to put some of these in RDS itself
# label=paste("sim_seed", sim_seed,
#       "sample_size", sample_size,
#       "sim_slope", sim_slope,
#       "sim_intercept", sim_intercept,
#       "sim_ses_coef", sim_ses_coef,
#       "sim_age_coef", sim_age_coef,
#       sep="_")

##saveRDS(results,paste0("results_",label,".RDS" ))
saveRDS(results,paste0("results_",abs(round(10000000*rnorm(1))), 
                       sample(letters)[1],abs(round(10000000*rnorm(1))), ".RDS" ))
}
