#' Run a simulation as displayed in the readMe
#'
#' Runs all the commands in the readMe to generate rmse for each method.
#' Run several instances of this to do a "simulation study"
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
readme_sim <- function(sim_seed=101, sample_size=1000, sim_slope=100, sim_intercept=12, 
                       sim_ses_coef=.01, sim_age_coef=.01){
#library(fda)
#library(ggplot2)
#library(reshape2)
#library(refund)
#library(nlme)
#library(devtools)
#install_github("swihart/wfpca")
#library(wfpca)
d<-prep_data()
#head(d)
over_samp_mat<-sample_data_fpc(d, sample_size, seed=sim_seed, timepoints=as.numeric(names(d[,-1])))
with_ses <- calculate_ses(over_samp_mat, slope=sim_slope, intercept=sim_intercept)
long <-make_long(with_ses)
#ggplot(subset(long, newid %in% 1:1000), aes(x=age,y=inches, colour=factor(newid)))+geom_line()
age_vec <- c(sort(unique(long$age)))
#head(long)
censored <- apply_censoring(long, ses_coef=sim_ses_coef, age_coef=sim_age_coef)
ddply(censored, .(age), function(w) sum(w$instudy==1) )
melt.prob.cens=ddply(censored, .(newid,age), function(w) w$prob.cens )
dcast.prob.cens=dcast(melt.prob.cens, newid~age, value.var="V1")
apply(dcast.prob.cens, 2, function(w) round(range(w),2))
head(censored,18)
## decided to add/sub 1/8 inch; "measurement error (?)"
censored$inches <- censored$inches + runif(length(censored$inches), -1/8,1/8)
observed_with_stipw <- calculate_stipw(censored,"omit")
wtd_trajectories <- calculate_wtd_trajectories(observed_with_stipw)
#head(wtd_trajectories)
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






## truth, all data:
true_avg <- ddply(long, .(age), function(w)  mean(w$inches))
true_avg$approach<- "true_avg"


## compare the bias of mean trajectories of the following approaches applied to the 
## censored/observed data, not the truth:
##a) naive-nonparametric: that is, just na.rm=TRUE and turn the mean() crank
naive_non_parm_avg <- ddply(observed_with_stipw, .(age), function(w)  mean(w$inches))
naive_non_parm_avg$approach <- "naive_non_parm_avg"
##b) weighted-nonparametric: a little more sophisticated, na.rm=TRUE and turn the weighted mean() crank 
wtd_non_parm_avg <- ddply(observed_with_stipw, .(age), function(w) sum(w$stipw*w$inches)/sum(w$stipw))
wtd_non_parm_avg$approach <- "wtd_non_parm_avg"
##c)  naive-FPC, 
age_vec <- c(sort(unique(long$age)))
unwtd_fncs <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches")
fpca_unwtd_fncs <- fpca.face(Y=as.matrix(unwtd_fncs[,-1]), argvals=age_vec, knots=26)
naive_fpc <- data.frame(age=age_vec, V1=fpca_unwtd_fncs$mu, approach="naive_fpc")
##d) naive-FPC-post-adjusted-by-weights, and 
pabw <- dcast(data=wtd_trajectories, formula= newid~age, value.var="stipw01.n")
pabw[pabw==0]<-NA
pabw_fpca_unwtd_fncs <- fpca_unwtd_fncs$Yhat*pabw[,-1] ## element-wise
naive_fpc_pabw_avg <- colMeans(pabw_fpca_unwtd_fncs,na.rm=TRUE)
naive_fpc_pabw <- data.frame(age=age_vec, V1=naive_fpc_pabw_avg, approach="naive_fpc_pabw")
##e) weighted-FPC. 
wtd_fncs <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches_wtd")
fpca_wtd_fncs <- fpca.face(Y=as.matrix(wtd_fncs[,-1]), argvals=age_vec, knots=26)
weighted_fpc <- data.frame(age=age_vec, V1=fpca_wtd_fncs$mu, approach="weighted_fpc")
##f) lme
#library(lme4)
#naive_lme<-lmer(inches ~ ns(age, df=5) + (1+age|newid), data=observed_with_stipw)
#predict(naive_lme, newdata=data.frame(age=age_vec), level=0)
library(nlme)
library(lme4)
naive_lme<-tryCatch(
{
  #naive_lme_model<-lme(inches ~ bs(age, df=15), random=~1|newid, data=observed_with_stipw);
  #naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), level=0), approach="naive_lme")
  
  naive_lme_model<-lmer(inches ~ bs(age, df=15) + (1|newid), data=observed_with_stipw);
  naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="naive_lme")

  },
  warning =function(cond){
    write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
    naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  ;
    naive_lme
  },
error =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  ;
  naive_lme
  })
summary(naive_lme)

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
  wtd_lme <- data.frame(age=age_vec, V1=predict(wtd_lme_model, newdata=data.frame(age=age_vec), re.form=~0), approach="wtd_lme")
  
},
warning =function(cond){
  write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  wtd_lme <- data.frame(age=age_vec, V1=cond, approach="wtd_lme")  ;
  wtd_lme
},
error =function(cond){
  write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  wtd_lme <- data.frame(age=age_vec, V1=cond, approach="wtd_lme")  ;
  wtd_lme
})
summary(wtd_lme)




## plot each approach's mean
means <- rbind(true_avg, naive_non_parm_avg, wtd_non_parm_avg, naive_fpc, naive_fpc_pabw, weighted_fpc, naive_lme, wtd_lme)
means$approach<-factor(means$approach, levels=unique(means$approach))
colnames(means)[colnames(means)=="V1"]<- "inches"
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




# library(ggplot2);
# ggplot(means, aes(x=age,y=inches, colour=approach))+geom_point()+geom_path()
# #zoom!
# ggplot(means, aes(x=age,y=inches, colour=approach))+geom_point()+geom_path()+coord_cartesian(xlim=c(14.9,18.1),ylim=c(68,75))
# #zoom! +facetting reveals overplotting
# naive_non_parm_avg == naive_fpc
#   wtd_non_parm_avg == weighted_fpc
#   naive_fpc_pabw is distinct but in ball park
##ggplot(means, aes(x=age,y=inches, colour=approach))+geom_point()+geom_path()+coord_cartesian(xlim=c(14.9,18.1),ylim=c(69,75)) + facet_grid(.~approach) 
## transform data for straightforward rmse calculations:
#one_row_per_age=dcast(means, age~approach, value.var="inches")
## calculate rmse against the true_avg
# a=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_non_parm_avg"])^2))
# b=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"wtd_non_parm_avg"])^2))
# c=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_fpc"])^2))
# d=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_fpc_pabw"])^2))
# e=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"weighted_fpc"])^2))
# f=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_lme"])^2))
# g=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"wtd_lme"])^2))
# results=data.frame(naive_non_parm_avg=a,
#                    wtd_non_parm_avg  =b,
#                    naive_fpc         =c,
#                    naive_fpc_pabw    =d,
#                    weighted_fpc      =e,
#                    naive_lme         =f,
#                    wtd_lme           =g,
#                    perc_ltfu_18      =percent.missing.at.age.18,
#                    sim_seed = sim_seed, 
#                    sample_size = sample_size,
#                    sim_slope = sim_slope,
#                    sim_intercept = sim_intercept, 
#                    sim_ses_coef = sim_ses_coef,
#                    sim_age_coef = sim_age_coef)
# 
# ## throw in missing at each age:
# results <- as.data.frame(t(unlist(c(results, as.data.frame(t(percent.missing))))))




## as opposed to calculating rmse against the true_avg as above;
## we just calculate squared error so we have it for each age below:
# a=(one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_non_parm_avg"])^2
# b=(one_row_per_age[,"true_avg"]-one_row_per_age[,"wtd_non_parm_avg"])^2
# c=(one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_fpc"])^2
# d=(one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_fpc_pabw"])^2
# e=(one_row_per_age[,"true_avg"]-one_row_per_age[,"weighted_fpc"])^2
# f=(one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_lme"])^2
# g=(one_row_per_age[,"true_avg"]-one_row_per_age[,"wtd_lme"])^2
# results=data.frame(age               = one_row_per_age[,"age"],
#                    naive_non_parm_avg=a,
#                    wtd_non_parm_avg  =b,
#                    naive_fpc         =c,
#                    naive_fpc_pabw    =d,
#                    weighted_fpc      =e,
#                    naive_lme         =f,
#                    wtd_lme           =g,
#                    perc_ltfu_18      =percent.missing.at.age.18,
#                    sim_seed = sim_seed, 
#                    sample_size = sample_size,
#                    sim_slope = sim_slope,
#                    sim_intercept = sim_intercept, 
#                    sim_ses_coef = sim_ses_coef,
#                    sim_age_coef = sim_age_coef,
#                    percent_missing = percent.missing,
#                    percent_missing_below_median = percent.missing.below.median,
#                    percent_missing_above_median = percent.missing.above.median)
# 

## throw in missing at each age:
##results <- as.data.frame(t(unlist(c(results, as.data.frame(t(percent.missing))))))


## instead of calculating errors at this stage, just return the predictions.
## do errors in the post-simulation fashion
# a=(one_row_per_age[,"true_avg"]
# b=one_row_per_age[,"wtd_non_parm_avg"]
# c=one_row_per_age[,"naive_fpc"]
# d=one_row_per_age[,"naive_fpc_pabw"]
# e=one_row_per_age[,"weighted_fpc"]
# f=one_row_per_age[,"naive_lme"]
# g=one_row_per_age[,"wtd_lme"]
# h=one_row_per_age[,"naive_non_parm_avg"]
# results=data.frame(age               = one_row_per_age[,"age"],
#                    true_avg          =a,
#                    naive_non_parm_avg=h,
#                    wtd_non_parm_avg  =b,
#                    naive_fpc         =c,
#                    naive_fpc_pabw    =d,
#                    weighted_fpc      =e,
#                    naive_lme         =f,
#                    wtd_lme           =g,
#                    perc_ltfu_18      =percent.missing.at.age.18,
#                    sim_seed = sim_seed, 
#                    sample_size = sample_size,
#                    sim_slope = sim_slope,
#                    sim_intercept = sim_intercept, 
#                    sim_ses_coef = sim_ses_coef,
#                    sim_age_coef = sim_age_coef,
#                    percent_missing = percent.missing,
#                    percent_missing_below_median = percent.missing.below.median,
#                    percent_missing_above_median = percent.missing.above.median)




# qplot(x=age_vec, y=wtd_lme$V1) + 
#   geom_point(aes(x=age_vec,y=true_avg$V1),color="red") + 
#   geom_point(aes(x=age_vec,y=true_avg$V1),color="blue")+
#   geom_point(aes(x=age_vec,y=naive_lme$V1),color="orange")
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
saveRDS(results,paste0("results_",abs(round(10000000*rnorm(1))), sample(letters)[1],abs(round(10000000*rnorm(1))), ".RDS" ))
}
