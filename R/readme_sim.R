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
library(fda)
library(ggplot2)
library(reshape2)
library(refund)
library(nlme)
library(devtools)
#install_github("swihart/wfpca")
library(wfpca)
d<-prep_data()
#head(d)
over_samp_mat<-sample_data(d, sample_size, seed=sim_seed)
with_ses <- calculate_ses(over_samp_mat, slope=sim_slope, intercept=sim_intercept)
long <-make_long(with_ses)
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
dcast.wtd.trajectories<-dcast(wtd_trajectories, newid~age, value.var="inches")
percent.missing.at.age.18=sum(is.na(dcast.wtd.trajectories["18"]))/length(unlist(dcast.wtd.trajectories["18"]))




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
tryCatch(
{naive_lme_model<-lme(inches ~ ns(age, df=5), random=~age|newid, data=observed_with_stipw)
##predict(naive_lme_model, newdata=data.frame(age=age_vec), level=0)
naive_lme <- data.frame(age=age_vec, V1=predict(naive_lme_model, newdata=data.frame(age=age_vec), level=0), approach="naive_lme")
},
warning ={
  naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  
},
error ={
  naive_lme <- data.frame(age=age_vec, V1=NA, approach="naive_lme")  
})




## plot each approach's mean
means <- rbind(true_avg, naive_non_parm_avg, wtd_non_parm_avg, naive_fpc, naive_fpc_pabw, weighted_fpc, naive_lme)
means$approach<-factor(means$approach, levels=unique(means$approach))
colnames(means)[colnames(means)=="V1"]<- "inches"
##library(ggplot2);
##ggplot(means, aes(x=age,y=inches, colour=approach))+geom_point()+geom_path()
## zoom!
##ggplot(means, aes(x=age,y=inches, colour=approach))+geom_point()+geom_path()+coord_cartesian(xlim=c(14.9,18.1),ylim=c(68,75))
## zoom! +facetting reveals overplotting
## naive_non_parm_avg == naive_fpc
##   wtd_non_parm_avg == weighted_fpc
##   naive_fpc_pabw is distinct but in ball park
##ggplot(means, aes(x=age,y=inches, colour=approach))+geom_point()+geom_path()+coord_cartesian(xlim=c(14.9,18.1),ylim=c(69,75)) + facet_grid(.~approach) 
## transform data for straightforward rmse calculations:
one_row_per_age=dcast(means, age~approach, value.var="inches")
## calculate rmse against the true_avg
a=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_non_parm_avg"])^2))
b=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"wtd_non_parm_avg"])^2))
c=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_fpc"])^2))
d=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_fpc_pabw"])^2))
e=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"weighted_fpc"])^2))
f=sqrt(mean((one_row_per_age[,"true_avg"]-one_row_per_age[,"naive_lme"])^2))
results=data.frame(naive_non_parm_avg=a,
                   wtd_non_parm_avg  =b,
                   naive_fpc         =c,
                   naive_fpc_pabw    =d,
                   weighted_fpc      =e,
                   naive_lme         =f,
                   perc_ltfu_18      =percent.missing.at.age.18,
                   sim_seed = sim_seed, 
                   sample_size = sample_size,
                   sim_slope = sim_slope,
                   sim_intercept = sim_intercept, 
                   sim_ses_coef = sim_ses_coef,
                   sim_age_coef = sim_age_coef)


## currently not using; going to put some of these in RDS itself
label=paste("sim_seed", sim_seed, 
      "sample_size", sample_size,
      "sim_slope", sim_slope,
      "sim_intercept", sim_intercept, 
      "sim_ses_coef", sim_ses_coef,
      "sim_age_coef", sim_age_coef,
      sep="_")

##saveRDS(results,paste0("results_",label,".RDS" ))
saveRDS(results,paste0("results_",abs(round(10000000*rnorm(1))), sample(letters)[1],abs(round(10000000*rnorm(1))), ".RDS" ))
}
