#' !!!!Achtung: Weighted methods are not done and are omitted.  BJS 5/20/2015  -- Run a simulation and then perform prediction.  Based on insample_sim() / insample_sim.R
#'
#' This function comes after many months of running readme_sim and talks with Bryan Lau on 
#' what we need to demonstrate for the method.
#'
#' @param sim_seed passed to set.seed()
#' @param sample_size (defaults to 1000) and must be a multiple of 100
#' @param sample_size2 (defaults to 1) the out of sample sample size
#' @param sim_slope see slope in calculate_ses()
#' @param sim_intercept see intercept in calculate_ses()
#' @param sim_ses_coef see ses_coef in apply_censoring()
#' @param sim_age_coef see age_coef in apply_censoring()
#' @export
#' @return results a data frame with rmse for each approach for that simulated dataset
#' @examples
#' ---
outsample_sim <- function(sim_seed=101, sample_size=1000, sample_size2=1, sim_slope=100,
                          sim_intercept=12, sim_ses_coef=.01, sim_age_coef=.01){
  ## quick start with the default values as variables for testing.  Loads packages.
  #test_prep_script()
  ## get the data prepped, that is simulated and censored and calculate missing. 
  #simulate_censor_summarize()

  ## just for testing; comment out before github commits
  sim_seed=101; sample_size=1000; sample_size2=1; sim_slope=100; sim_intercept=12; sim_ses_coef=.01; sim_age_coef=.01;

  
  
## This is really, really inefficient, redundant coding, but so be it.  I need a "data0" and a "data1".
## So, I take the code that generates data from insample_sim() and copy and paste and then put a "0" suffix
## on the data that will generate the models and "1" suffix on those for the out of sample.  
## Ah, let me take on moment to see if I can functionalize this bastard.


##http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value/15140507#15140507
##https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

## here are data that will be fitted, and whose fittings will be applied to a new, out of sample dataset
list[
  wtd_trajectories,
  unwtd_fncs,
  wtd_fncs,
  observed_with_stipw,
  all_with_stipw,
  censored,
  age_vec
  ] <- all_things_sim_data(sim_seed, sample_size, sim_slope,
                           sim_intercept, sim_ses_coef, sim_age_coef)


## the new, out of sample dataset objects (all with suffix "2"):
list[
  wtd_trajectories2,
  unwtd_fncs2,
  wtd_fncs2,
  observed_with_stipw2,
  all_with_stipw2,
  censored2,
  age_vec2
  ] <- all_things_sim_data(sim_seed, sample_size, sim_slope,
                           sim_intercept, sim_ses_coef, sim_age_coef)


## now snake out sample_size2 people... and change newid
newid2 <- sample(1:sample_size, sample_size2)
list[
  wtd_trajectories2,
  unwtd_fncs2,
  wtd_fncs2,
  observed_with_stipw2,
  all_with_stipw2,
  censored2,
  age_vec2
  ] <- 
  list(
    wtd_trajectories2[newid %in% newid2],
    unwtd_fncs2[newid %in% newid2],
    wtd_fncs2[newid %in% newid2],
    observed_with_stipw2[newid %in% newid2],
    all_with_stipw2[newid %in% newid2],
    censored2[newid %in% newid2],
    age_vec2
  ) 
dim(wtd_fncs2)
## change newid
## max.newid <- 22200000 + wtd_trajectories[,max(newid)]
newid2.addon <- 22200000
wtd_trajectories2[, newid := newid+newid2.addon]
unwtd_fncs2[, newid := newid+newid2.addon]
wtd_fncs2[, newid := newid+newid2.addon]
observed_with_stipw2[, newid := newid+newid2.addon]
all_with_stipw2[, newid := newid+newid2.addon]
censored2[, newid := newid+newid2.addon]
## don't know why I lost key but put it back
     setkeyv(all_with_stipw2, key(all_with_stipw))
setkeyv(observed_with_stipw2, key(observed_with_stipw))
   setkeyv(wtd_trajectories2, key(wtd_trajectories))
           setkeyv(wtd_fncs2, key(wtd_fncs))
         setkeyv(unwtd_fncs2, key(unwtd_fncs))
           setkeyv(censored2, key(censored))





## standardize the names selected across all approaches and their resultant  data sets
selected<-c("newid","age", "ses", 
            
            "inches", "inches_wtd", "inches_ltfu", "inches_predicted",
            
            "remaining",
            "stipw",
            "stipw01",
            "stipw01.n",
            
            "approach")


##c)  naive-FPC,
## Here we get the mu and the rho (where rho are eigenfunctions):
fpca_unwtd_fncs <- fpca.face(Y=as.matrix(unwtd_fncs)[,-1], argvals=age_vec, knots=26)
naive_fpc       <- data.frame(age=age_vec, V1=fpca_unwtd_fncs$mu, approach="naive_fpc")

mu  <- fpca_unwtd_fncs$mu
rho <- fpca_unwtd_fncs$eigenvectors

## set up new scores matrix (corresponds to dataset objects suffix "2")
scores2 <- matrix(NA, nrow=nrow(unwtd_fncs2), ncol=ncol(rho))

## for loop for now; calculate new scores!
ii <-1
for(i in unwtd_fncs2$newid){
  ## let's take one out of sample trajectory:
  y_i <- unwtd_fncs2[newid == i,-1,with=FALSE]
  ## now across the rows of rho  ;-p 
  for(j in 1:ncol(scores2)){
    integrand_i <- as.numeric((y_i - mu)*rho[,j])
    scores2[ii,j] <- simpson_nint(x=age_vec, fx=integrand_i, n.pts=256)$value
  }
ii <- ii + 1 ## i know this is awful coding but just get it going
}
## put 'em together:
fpca_unwtd_fncs$Yhat2 <- t( rho%*%t(scores2) + mu )

## quick checks on `fpca_unwtd_fncs$Yhat2`:
theguy <- 1
# predicted for the guy in blackdots
qplot(x=age_vec, y=fpca_unwtd_fncs$Yhat2[theguy,]) + 
 ## cyan line: the first guy (different guy, but similar based on SES etc) in the original data set
 geom_path(y=fpca_unwtd_fncs$Yhat[ 734,], color="cyan", size=1.5) +
 ## truth for the guy in pink 
 geom_path(y=as.numeric(unwtd_fncs2[ theguy,-1,with=F]), color="pink", size=1.5) 


## the above approach doesn't return a bunch of NAs, like this Y.pred attempt:
##
## Below doesn't work
##
## using unwtd_fncs eigenvalues, approximate Y.pred
# fpca_unwtd_fncs <- fpca.face(Y      = as.matrix(unwtd_fncs )[,-1], 
#                              Y.pred = as.matrix(unwtd_fncs2)[,-1], 
#                              argvals=age_vec, knots=26)
# 
# naive_fpc       <- data.frame(age=age_vec, V1=fpca_unwtd_fncs$mu, approach="naive_fpc")
##
## Above doesn't work  -- limited to only those with full data.
##

## combine with long for the prediction:
## a little different than previous examples; now we do have individual level curves
## need to extract them (Yhat) and rename them and data.table them
naive_fpc_indiv <- as.data.frame(cbind(unwtd_fncs2$newid,
                                       fpca_unwtd_fncs$Yhat2))
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
##setkey(all_with_stipw2, newid, age)
naive_fpc  <- naive_fpc[all_with_stipw2]
setkey(naive_fpc, newid, age)
naive_fpc <- subset(naive_fpc, select=selected)

## visual checks
#ggplot(data=naive_fpc, aes(x=age, y=inches_predicted, color=factor(newid)))+geom_point()+
#  geom_point(aes(y=inches), color='black')


# 
# ## PLEASE SEE TBD COMMENTS BELOW FOR WEIGHTED-FPC APPROACH.
# ## CURRENTLY ON HOLD FOR OUTSAMPLE_SIM.R   --- BJS 2015/05/20
# ##e) weighted-FPC.
# ## Here we get the mu and the rho (where rho are eigenfunctions):
# fpca_wtd_fncs <- fpca.face(Y=as.matrix(wtd_fncs)[,-1], argvals=age_vec, knots=26)
# weighted_fpc <- data.frame(age=age_vec, V1=fpca_wtd_fncs$mu, approach="weighted_fpc")
# 
# 
# mu  <- fpca_wtd_fncs$mu
# rho <- fpca_wtd_fncs$eigenvectors
# 
# 
# ## TBD: grab for loop for scores from naive (above) and rename stuff:
# ## TBD: grab stuff just below for loop for scores but above/cinluding naive_fpc[all_with_stipw2]
# ##      and rename appropriately.
# ## TBD:  then iron out the weights with one of the two options (chunks) below 
# ##         The hold up was how to generate weights for out of sample people --
#   
# 
# 
# 
# 
# 
# ## should be okay above this line...
# ## should be okay above this line...
# ## should be okay above this line...
# ## should be okay above this line...
# ## should be okay above this line...
# ## should be okay above this line...
# 
# ## JAM!  --> which weights do I use?  Might be useful to think of 
# ##           this simulation as original data:  100 people, new data 1 person.
# ##       --> could rewrite inputs to reflect this...
# 
# 
# 
# ## begin extra steps: (weighted inputs average out for mean, but need to be de-weighted for individual)
# setDT(wtd_trajectories)
# wts <- subset(wtd_trajectories, select=c("newid","age","stipw01.n"))
# setkey(wts, newid, age)
# ## can't do curve completion with these three knuckleheads
# #weighted_fpc.test<-weighted_fpc[wts]
# #weighted_fpc.test[, inches_predicted_weighted:= inches_predicted]
# #weighted_fpc.test[, inches_predicted_deweighted:= inches_predicted/stipw01.n]
# weighted_fpc.test<-wts[weighted_fpc]
# ## see how many 0's
# weighted_fpc.test[,table(round(stipw01.n,2))]
# ## change 0's to NA
# weighted_fpc.test[stipw01.n==0, stipw01.n:=NA]
# weighted_fpc.test[, inches_predicted_weighted:= inches_predicted]
# weighted_fpc.test[!is.na(stipw01.n), inches_predicted_deweighted:= inches_predicted/stipw01.n]
# ## get the wtd_population_mean in there (currently from original data):
# weighted_fpc.test[,wtd_pop_mean:=fpca_wtd_fncs$mu,by=newid]
# ##for now, if don't have observed data there, we just imputed the weighted mean
# ## kinda lame, think on it.
# weighted_fpc.test[is.na(stipw01.n), inches_predicted_deweighted:= wtd_pop_mean, by=newid]
# ## end extra steps:
# setkey(weighted_fpc.test, newid, age)#, inches_predicted_deweighted)
# setnames(weighted_fpc.test, "inches_predicted", "inches_predicted_old")
# setnames(weighted_fpc.test, "inches_predicted_deweighted", "inches_predicted")
# weighted_fpc <- subset(weighted_fpc.test, select=selected)
# 
# 
# ## old stuff (from insample_sim.R):
# ## compare to chunk above, where  i started editing...
# ## begin extra steps: (weighted inputs average out for mean, but need to be de-weighted for individual)
# setDT(wtd_trajectories)
# wts <- subset(wtd_trajectories, select=c("newid","age","stipw01.n"))
# setkey(wts, newid, age)
# ## can't do curve completion with these three knuckleheads
# #weighted_fpc.test<-weighted_fpc[wts]
# #weighted_fpc.test[, inches_predicted_weighted:= inches_predicted]
# #weighted_fpc.test[, inches_predicted_deweighted:= inches_predicted/stipw01.n]
# weighted_fpc.test<-wts[weighted_fpc]
# ## see how many 0's
# weighted_fpc.test[,table(round(stipw01.n,2))]
# ## change 0's to NA
# weighted_fpc.test[stipw01.n==0, stipw01.n:=NA]
# weighted_fpc.test[, inches_predicted_weighted:= inches_predicted]
# weighted_fpc.test[!is.na(stipw01.n), inches_predicted_deweighted:= inches_predicted/stipw01.n]
# ## get the wtd_population_mean in there:
# weighted_fpc.test[,wtd_pop_mean:=fpca_wtd_fncs$mu,by=newid]
# ##for now, if don't have observed data there, we just imputed the weighted mean
# ## kinda lame, think on it.
# weighted_fpc.test[is.na(stipw01.n), inches_predicted_deweighted:= wtd_pop_mean, by=newid]
# ## end extra steps:
# setkey(weighted_fpc.test, newid, age)#, inches_predicted_deweighted)
# setnames(weighted_fpc.test, "inches_predicted", "inches_predicted_old")
# setnames(weighted_fpc.test, "inches_predicted_deweighted", "inches_predicted")
# weighted_fpc <- subset(weighted_fpc.test, select=selected)
# 
# 
# # ggplot(data=weighted_fpc[newid %in% c(1,2,5)], aes(x=age, y=inches_predicted, color=factor(newid)))+geom_point()+
# #   geom_point(aes(y=inches_ltfu), color='black')
# 


##f) lme
naive_lme<-tryCatch(
{
  #naive_lme_model<-lme(inches ~ bs(age, df=15), random=~1|newid, data=observed_with_stipw);
  #naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), level=0), approach="naive_lme")
  
  naive_lme_model<-lmer(inches ~ bs(age, df=15) + (1|newid), data=observed_with_stipw);
  ## re.form=~0
  #naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="naive_lme")
  ## re.form=~1
  naive_lme <- data.table(newid            = all_with_stipw2$newid,
                          age              = all_with_stipw2$age,
                          ses              = all_with_stipw2$ses,
                          inches           = all_with_stipw2$inches,
                          inches_wtd       = all_with_stipw2$inches_wtd,
                          inches_ltfu      = all_with_stipw2$inches_ltfu,
                          inches_predicted = predict(naive_lme_model, 
                                                     newdata=all_with_stipw2, re.form=~0),
                          remaining        = all_with_stipw2$remaining,
                          stipw            = all_with_stipw2$stipw,
                          stipw01          = all_with_stipw2$stipw01,
                          stipw01.n        = all_with_stipw2$stipw01.n,
                          approach="naive_lme")
},
warning =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  ##naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  ;
  naive_lme <- data.table(newid            = all_with_stipw2$newid,
                          age              = all_with_stipw2$age,
                          ses              = all_with_stipw2$ses,
                          inches           = all_with_stipw2$inches,
                          inches_wtd       = all_with_stipw2$inches_wtd,
                          inches_ltfu      = all_with_stipw2$inches_ltfu,
                          inches_predicted = NA,
                          remaining        = all_with_stipw2$remaining,
                          stipw            = all_with_stipw2$stipw,
                          stipw01          = all_with_stipw2$stipw01,
                          stipw01.n        = all_with_stipw2$stipw01.n,
                          approach="naive_lme")
},
error =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  #naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  ;
  naive_lme <- data.table(newid            = all_with_stipw2$newid,
                          age              = all_with_stipw2$age,
                          ses              = all_with_stipw2$ses,
                          inches           = all_with_stipw2$inches,
                          inches_wtd       = all_with_stipw2$inches_wtd,
                          inches_ltfu      = all_with_stipw2$inches_ltfu,
                          inches_predicted = NA,
                          remaining        = all_with_stipw2$remaining,
                          stipw            = all_with_stipw2$stipw,
                          stipw01          = all_with_stipw2$stipw01,
                          stipw01.n        = all_with_stipw2$stipw01.n,
                          approach="naive_lme")
})
setkeyv(naive_lme, key(naive_fpc))
## quick checks:
# summary(naive_lme)
# ggplot(data=naive_lme[newid %in% c(2)], 
#       aes(x=age, y=inches_predicted, color=factor(newid)))+
#   geom_point()+
#   geom_point(aes(y=inches), color='black')



## PLEASE NOTE THAT THE WTD METHODS ARE ON HOLD. PLEASE SEE WTD_FPC COMMENTS ABOVE.
## HOW DO I GENERATE STIPW FOR OUT OF SCAMPLE.  ALSO, IF I DECIDE TO NOT PURSUE
## THIS KIND OF WT-ING, THEN NO NEED TO DUMP IN HOURS AND HOURS OF WORK.
## ON  HOLD FOR NOW.  -BJS 2015/05/20
# wtd_lme<-tryCatch(
# {
#   #   wtd_lme_model<-lme(inches_wtd ~ bs(age, df=15), random=~1|newid, data=wtd_trajectories,
#   #                      na.action=na.omit);
#   wtd_lme_model<-lmer(inches_wtd ~ bs(age, df=15) + (1|newid), data=wtd_trajectories,
#                       na.action=na.omit);
#   ##predict(wtd_lme_model, newdata=data.frame(age=age_vec), level=0)
#   ##wtd_lme <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), level=0), approach="wtd_lme")
#   ## note: below, use re.form=~0 in lme4:predict is equivalent to level=0 in nlme:lme
#   ## re.form=~0
#   ##wtd_lme <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="wtd_lme")
#   wtd_lme_pop_mean <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="wtd_lme")
#   
#   ## re.form=~1
#   wtd_lme1  <- data.table(newid            = all_with_stipw$newid,
#                           age              = all_with_stipw$age,
#                           ses              = all_with_stipw$ses,
#                           inches           = all_with_stipw$inches,
#                           inches_wtd       = all_with_stipw$inches_wtd,
#                           inches_ltfu      = all_with_stipw$inches_ltfu,
#                           inches_predicted = predict(wtd_lme_model, 
#                                                      newdata=all_with_stipw),
#                           remaining        = all_with_stipw$remaining,
#                           stipw            = all_with_stipw$stipw,
#                           stipw01          = all_with_stipw$stipw01,
#                           stipw01.n        = all_with_stipw$stipw01.n,
#                           approach="wtd_lme")
#   
#   wtd_lme.test<- wts[wtd_lme1]
#   ## see how many 0's
#   wtd_lme.test[,table(round(stipw01.n,2))]
#   ## change 0's to NA
#   wtd_lme.test[stipw01.n==0, stipw01.n:=NA]
#   wtd_lme.test[, inches_predicted_weighted:= inches_predicted]
#   wtd_lme.test[!is.na(stipw01.n), inches_predicted_deweighted:= inches_predicted/stipw01.n]
#   ## get the wtd_population_mean in there:
#   wtd_lme.test[,wtd_pop_mean:=wtd_lme_pop_mean$V1,by=newid]
#   ##for now, if don't have observed data there, we just imputed the weighted mean
#   ## kinda lame, think on it.
#   wtd_lme.test[is.na(stipw01.n), inches_predicted_deweighted:= wtd_pop_mean, by=newid]
#   ## end extra steps:
#   setkey(wtd_lme.test, newid, age)#, inches_predicted_deweighted)
#   setnames(wtd_lme.test, "inches_predicted", "inches_predicted_old")
#   setnames(wtd_lme.test, "inches_predicted_deweighted", "inches_predicted")
#   wtd_lme <- subset(wtd_lme.test, select=selected)  
#   
#   },
# warning =function(cond){
#   write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
#   #wtd_lme <- data.frame(age=age_vec, V1=cond, approach="wtd_lme")  ;
#   wtd_lme   <- data.table(newid            = all_with_stipw$newid,
#                           age              = all_with_stipw$age,
#                           ses              = all_with_stipw$ses,
#                           inches           = all_with_stipw$inches,
#                           inches_wtd       = all_with_stipw$inches_wtd,
#                           inches_ltfu      = all_with_stipw$inches_ltfu,
#                           inches_predicted = NA,
#                           remaining        = all_with_stipw$remaining,
#                           stipw            = all_with_stipw$stipw,
#                           stipw01          = all_with_stipw$stipw01,
#                           stipw01.n        = all_with_stipw$stipw01.n,
#                           approach="wtd_lme")
# },
# error =function(cond){
#   write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
#   #wtd_lme <- data.frame(age=age_vec, V1=cond, approach="wtd_lme")  ;
#   wtd_lme   <- data.table(newid            = all_with_stipw$newid,
#                           age              = all_with_stipw$age,
#                           ses              = all_with_stipw$ses,
#                           inches           = all_with_stipw$inches,
#                           inches_wtd       = all_with_stipw$inches_wtd,
#                           inches_ltfu      = all_with_stipw$inches_ltfu,
#                           inches_predicted = NA,
#                           remaining        = all_with_stipw$remaining,
#                           stipw            = all_with_stipw$stipw,
#                           stipw01          = all_with_stipw$stipw01,
#                           stipw01.n        = all_with_stipw$stipw01.n,
#                           approach="wtd_lme")
# })
# ## quick checks:
# # summary(wtd_lme)
# # ggplot(data=wtd_lme[newid %in% c(1,2,5)], 
# #       aes(x=age, y=inches_predicted, color=factor(newid)))+
# #   geom_point()+
# #   geom_point(aes(y=inches_ltfu), color='black')
# # ## see the original predictions, before 'deweighting'
# # ggplot(data=wtd_lme.test[newid %in% c(1,2,5)], 
# #        aes(x=age, y=inches_predicted_old, color=factor(newid)))+
# #   geom_point()+
# #   geom_point(aes(y=inches_ltfu), color='black')
# # 
# 
# 


#rbind it!!!!!!!!!!!!!



## plot each approach's mean
## means <- rbind(true_avg, naive_non_parm_avg, wtd_non_parm_avg, naive_fpc, naive_fpc_pabw, weighted_fpc, naive_lme, wtd_lme)
# note: we're simplifying to (w)fpca and (w)lme, so we rbind() fewer than line above (if you decided to add more approaches later,
#         make sure they have same format):
##means <- rbind(naive_fpc, weighted_fpc, naive_lme, wtd_lme)
## NOTE:  WTD_METHODS ON HOLD (SEE COMMENTS ABOVE FOR WTD_LME AND WTD_FPC)
means <- rbind(naive_fpc, naive_lme)

# from older files:
# means$approach<-factor(means$approach, levels=unique(means$approach))
# colnames(means)[colnames(means)=="V1"]<- "inches"



## this is a useful chunk to check the range of 
## probality of being censored (prob.cens)
#   ddply(censored, .(age), function(w) sum(w$instudy==1) )
#   melt.prob.cens=ddply(censored, .(newid,age), function(w) w$prob.cens )
#   dcast.prob.cens=dcast(melt.prob.cens, newid~age, value.var="V1")
#   apply(dcast.prob.cens, 2, function(w) round(range(w),2))
#   head(censored,18)


censored  <- as.data.frame(censored)
censored2 <- as.data.frame(censored2)
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
