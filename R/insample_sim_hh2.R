#' Hopkins Hybrid 1: Run a simulation and then perform prediction.  Based on insample_sim().  Censoring is
#' fabricated a la the SES of the Growth data.  We will make a _hh2 to do a more realistic censoring.
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
#' @param censoring_coef used to multiply the prob.cens from the simplog in hh2 runs. Default value is 1.  will try 0.95 and 1.05 for range.
#' @param hh_rds the path in single quotes to the local copy of hopkins_hybrid.RDS
#' @export
#' @return results a data frame with rmse for each approach for that simulated dataset
#' @examples
#' ---
insample_sim_hh2 <- function(sim_seed=101, sample_size=1000, sim_slope=100,
                                  sim_intercept=12, sim_ses_coef=.01, sim_age_coef=.01, censoring_coef=1,
                            hh_rds='./data_local/hopkins_hybrid.RDS', 
                            hh_rds_simplog='./data_local/hopkins_hybrid_prob_ltfu_coeffs.RDS', 
                            hh_rds_long='./data_local/hopkins_hybrid_long.RDS'){
  ## quick start with the default values as variables for testing.  Loads packages.
  #test_prep_script()
  ## get the data prepped, that is simulated and censored and calculate missing.
  #simulate_censor_summarize()

  ## just for testing; comment out before github commits
  #library(devtools); library(roxygen2); install(); document(); 
  #sample_size=1000; censoring_coef=1.50; hh_rds_simplog='./data_local/hopkins_hybrid_prob_ltfu_coeffs.RDS';hh_rds_long='./data_local/hopkins_hybrid_long.RDS'; hh_rds='./data_local/hopkins_hybrid.RDS'; sim_seed=101;  sim_slope=100; sim_intercept=12; sim_ses_coef=.01; sim_age_coef=.01;
  
  ## d will be a 187x13 matrix based on the extraction of ./data_local/hopkins_hybrid_prep.R
  ## we oversample d based on the fpc as well extract out the times of measurement
  d<-readRDS(hh_rds)
  time_vec <- c(as.numeric(colnames(d[,-1,with=FALSE])))
  over_samp_mat<-sample_data_fpc(as.matrix(d), sample_size, seed=sim_seed, timepoints=time_vec)

  ## we calculate ses on the oversampled dataset and turn the matrix into a long dataset.
  #with_ses <- calculate_ses(over_samp_mat, slope=sim_slope, intercept=sim_intercept)
  #long <-make_long(with_ses)
  ##
  ## previous two steps were for BG and hh1.  Here, we need to be more
  ## sophisticated about covariates and their relation to censoring.
  ## 1.  read in long, which has covariates in long format
  ## 2.  take unique() 
  ## 3.  Sample, then rbind covars with over_samp_mat (akin to `with_ses`)
  ## 4.  make it long (akin to `long` from `make_long(with_ses)`)
  ## implement:
  ## 1.
  long2<- readRDS(hh_rds_long)
  setDT(long2)
  setkey(long2, newid)
  long2[,id:=.GRP, by=newid]
  ## 2.
  covars <- subset(unique(long2),
                   select=c("id", "age",
                            "sex", "race",
                            "hetero","msm","ivdu"))
  ## 3.  
  covars_over_sample_ids <- sample(covars$id, 1000, replace=TRUE)
  covars_over_sample_mat <- covars[covars_over_sample_ids]
  with_covars <- cbind(over_samp_mat[,-1], covars_over_sample_mat[,-id, with=FALSE])
  ## 4.
  long <- melt(with_covars, id.vars=c("id", "age",
                                           "sex", "race",
                                           "hetero","msm","ivdu"),
                    variable="time", value="cd4", value.name="cd4")
  long$time <- as.numeric(long$time)
  long$newid <- long$id
  
  ## In this chunk we apply censoring.
  ## censored <- apply_censoring(long, ses_coef=sim_ses_coef, time_coef=sim_time_coef, protected=1:4)

  ## Measurement error: within 1/8 inch.  Can comment out or change.
  ##censored$cd4 <- censored$cd4 + runif(length(censored$cd4), -1/8,1/8)
  ## step above was developed for BG and hh1.
  ## steps below apply the simple logit glm in hopkins_hybrid_prep.R
  simplog<- readRDS(hh_rds_simplog)
  
  whamcast <- long
  setDT(whamcast)
  
  ## calculate the prob.cens multiplied by the censoring_coef
  whamcast[, prob.cens:=predict.glm(simplog, newdata=whamcast, type="response")*censoring_coef]
  whamcast[, summary(prob.cens)]
  whamcast[, subject_time_index:={ c(1:.N)}, by="newid"]
  protected <- 1:4
  whamcast[subject_time_index %in% protected, prob.cens:=0]
  whamcast[, instudy.sim:= 1-rbinom(.N,1,prob.cens)]
  whamcast[, instudy:=cumprod(instudy.sim), by="newid"]
  
  
  censored <- as.data.frame(whamcast)
  
  
  ## observed_with_stipw has the standardized inverse probability weights (stipw)
  ## wtd_trajectories has same info as observed_with_stipw but has cd4_wtd_hadamard as well
  ## we calculate all_with_stipw with NAs
       all_with_stipw <- calculate_stipw_hh2(censored,"keep")
  observed_with_stipw <- calculate_stipw_hh2(censored,"omit")
  wtd_trajectories    <- calculate_wtd_trajectories_hh2(observed_with_stipw)




  ## use data.table where possible speeds up future rbinds post simulation
#setDT(long)
#setkey(long, id, time)

setDT(all_with_stipw)
setkey(all_with_stipw, newid, time, cd4)

setDT(wtd_trajectories)
setkey(wtd_trajectories, newid, time, cd4)

interim <- wtd_trajectories[all_with_stipw]


## BEGIN:  new step...DEAN:
key.interim <- unique(interim[,c("time","remaining","denom"), with=FALSE])[!is.na(remaining) & !is.na(remaining)]
setkey(key.interim, time)
setkey(interim, time)
holder<-key.interim[interim]
holder[is.na(stipw01.n), stipw01.n := remaining*(stipw2/denom) ]
## END:  new step...DEAN:


## Dean step: change below to holder (formerly interim)
all_with_stipw<-subset(holder,
                       select=c(names(all_with_stipw),
                                "cd4_wtd_hadamard",
                                "remaining",
                                "stipw",
                                "stipw01",
                                "stipw01.n"))


## DEAN step: reset keys;
setkey(all_with_stipw, newid, time)
setkey(wtd_trajectories, newid, time)



all_with_stipw[            , cd4_ltfu:=cd4 ]
all_with_stipw[is.na(stipw), cd4_ltfu:=NA ]

## calculate "remean" data prep:
all_with_stipw[,
               cd4_wtd_remean:= cd4_ltfu -
                                   mean(cd4_ltfu, na.rm=T) +
                                   mean(cd4_wtd_hadamard , na.rm=T),
               by=time]


## standardize the names selected across all approaches and their resultant  data sets
selected<-c("newid","time", 

            "cd4", "cd4_wtd_hadamard","cd4_wtd_remean", "cd4_ltfu", "cd4_predicted",

            "remaining",
            "stipw",
            "stipw01",
            "stipw01.n",

            "approach")



##c)  naive-FPC,
unwtd_fncs      <- dcast(data=wtd_trajectories, formula= newid~time, value.var="cd4")
naive_fpc_proc_minutes <- system.time(
  fpca_unwtd_fncs <- fpca.face(Y=as.matrix(unwtd_fncs)[,-1], argvals=time_vec, knots=7)
)[3]/60
naive_fpc       <- data.frame(time=time_vec, V1=fpca_unwtd_fncs$mu, approach="naive_fpc")
## combine with long for the prediction:
## a little different than previous examples; now we do have individual level curves
## need to extract them (Yhat) and rename them and data.table them
naive_fpc_indiv <- as.data.frame(cbind(1:nrow(fpca_unwtd_fncs$Yhat),fpca_unwtd_fncs$Yhat))
colnames(naive_fpc_indiv) <- colnames(unwtd_fncs)
setDT(naive_fpc_indiv)
naive_fpc <- melt(naive_fpc_indiv,
                  id.vars=c("newid"),
                  variable.name = "time",
                  variable.factor=FALSE,
                  value.name="cd4_predicted")
naive_fpc[,approach:="naive_fpc",]
naive_fpc[,time:=as.numeric(time)]
setDT(naive_fpc)
setkey(naive_fpc, newid, time)
naive_fpc  <- naive_fpc[all_with_stipw]
setkey(naive_fpc, newid, time)
naive_fpc <- subset(naive_fpc, select=selected)
## add these on post dean:  minutes and number of principle components (npc)
naive_fpc[, minutes:=naive_fpc_proc_minutes]
naive_fpc[, number_pc:=fpca_unwtd_fncs$npc]



# visual checks
# ggplot(data=naive_fpc[newid %in% c(1,2,5)], 
#        aes(x=time, y=cd4_predicted, color=factor(newid)))+
#   geom_path()+
#   geom_point()+
#   geom_line(aes(y=cd4_ltfu, id=factor(newid)), color='black' )+
#   geom_point(aes(y=cd4_ltfu, id=factor(newid), shape=factor(newid)), color='black')




##e) remean - weighted-FPC.
wtd_remean_fncs <- dcast(data=subset(all_with_stipw,
                                     select=c("newid","time","cd4_wtd_remean")),
                         formula= newid~time, value.var="cd4_wtd_remean")
weighted_remean_fpc_proc_minutes <- system.time(
  fpca_wtd_remean_fncs <- fpca.face(Y=as.matrix(wtd_remean_fncs)[,-1], 
                                    argvals=time_vec, knots=7)
)[3]/60
#weighted_fpc <- data.frame(time=time_vec, V1=fpca_wtd_fncs$mu, approach="weighted_fpc")
## combine with long for the prediction:
## a little different than previous examples; now we do have individual level curves
## need to extract them (Yhat) and rename them and data.table them
weighted_remean_fpc_indiv <- as.data.frame(cbind(1:nrow(fpca_wtd_remean_fncs$Yhat),
                                                 fpca_wtd_remean_fncs$Yhat))
colnames(weighted_remean_fpc_indiv) <- colnames(wtd_remean_fncs)
setDT(weighted_remean_fpc_indiv)
weighted_remean_fpc <- melt(weighted_remean_fpc_indiv,
                     id.vars=c("newid"),
                     variable.name = "time",
                     variable.factor=FALSE,
                     value.name="cd4_predicted")
weighted_remean_fpc[,approach:="wtd_remean_fpc",]
weighted_remean_fpc[,time:=as.numeric(time)]
setDT(weighted_remean_fpc)
setkey(weighted_remean_fpc, newid, time)
weighted_remean_fpc <- weighted_remean_fpc[all_with_stipw]
setkey(weighted_remean_fpc, newid, time)
weighted_remean_fpc <- subset(weighted_remean_fpc, select=selected)
## add these on post dean:  minutes and number of principle components (npc)
weighted_remean_fpc[, minutes:=weighted_remean_fpc_proc_minutes]
weighted_remean_fpc[, number_pc:=fpca_wtd_remean_fncs$npc]
## visual checks
# ggplot(data=weighted_remean_fpc[newid %in% c(1,2,5)], 
#        aes(x=time, y=cd4_predicted, color=factor(newid)))+
#   geom_path()+
#   geom_point()+
#   geom_line(aes(y=cd4_ltfu, id=factor(newid)), color='black' )+
#   geom_point(aes(y=cd4_ltfu, id=factor(newid), shape=factor(newid)), color='black')


##e) weighted-FPC.
wtd_fncs <- dcast(data=wtd_trajectories, formula= newid~time, value.var="cd4_wtd_hadamard")
weighted_fpc_proc_minutes <- system.time(
  fpca_wtd_fncs <- fpca.face(Y=as.matrix(wtd_fncs)[,-1], argvals=time_vec, knots=7)
)[3]/60
#weighted_fpc <- data.frame(time=time_vec, V1=fpca_wtd_fncs$mu, approach="weighted_fpc")
## combine with long for the prediction:
## a little different than previous examples; now we do have individual level curves
## need to extract them (Yhat) and rename them and data.table them
weighted_fpc_indiv <- as.data.frame(cbind(1:nrow(fpca_wtd_fncs$Yhat),fpca_wtd_fncs$Yhat))
colnames(weighted_fpc_indiv) <- colnames(wtd_fncs)
setDT(weighted_fpc_indiv)
weighted_fpc <- melt(weighted_fpc_indiv,
                     id.vars=c("newid"),
                     variable.name = "time",
                     variable.factor=FALSE,
                     value.name="cd4_predicted")
weighted_fpc[,approach:="wtd_hadamard_fpc",]
weighted_fpc[,time:=as.numeric(time)]
setDT(weighted_fpc)
setkey(weighted_fpc, newid, time)
weighted_fpc <- weighted_fpc[all_with_stipw]
## begin extra steps: (weighted inputs avertime out for mean, but need to be de-weighted for individual)
setDT(wtd_trajectories)
## pre-DEAN: wts <- subset(wtd_trajectories, select=c("newid","time","stipw01.n"))
## post-DEAN: subset all_with_stipw
wts <- subset(all_with_stipw, select=c("newid","time","stipw01.n"))
setkey(wts, newid, time)

## can't do curve completion with these three knuckleheads
#weighted_fpc.test<-weighted_fpc[wts]
#weighted_fpc.test[, cd4_predicted_weighted:= cd4_predicted]
#weighted_fpc.test[, cd4_predicted_deweighted:= cd4_predicted/stipw01.n]

weighted_fpc.test<-wts[weighted_fpc]
## see how many 0's
weighted_fpc.test[,table(round(stipw01.n,2))]
## change 0's to NA
weighted_fpc.test[stipw01.n==0, stipw01.n:=NA]
weighted_fpc.test[, cd4_predicted_weighted:= cd4_predicted]
weighted_fpc.test[!is.na(stipw01.n), cd4_predicted_deweighted:= cd4_predicted/stipw01.n]
## get the wtd_population_mean in there:
## skip this time:  ##weighted_fpc.test[,wtd_pop_mean:=fpca_wtd_fncs$mu,by=newid]
##for now, if don't have observed data there, we just imputed the weighted mean
## kinda lame, think on it.
## skip this time:  ##weighted_fpc.test[is.na(stipw01.n), cd4_predicted_deweighted:= wtd_pop_mean, by=newid]
## end extra steps:
setkey(weighted_fpc.test, newid, time)#, cd4_predicted_deweighted)
setnames(weighted_fpc.test, "cd4_predicted", "cd4_predicted_old")
setnames(weighted_fpc.test, "cd4_predicted_deweighted", "cd4_predicted")
weighted_fpc <- subset(weighted_fpc.test, select=selected)
## add these on post dean:  minutes and number of principle components (npc)
weighted_fpc[, minutes:=weighted_fpc_proc_minutes]
weighted_fpc[, number_pc:=fpca_wtd_fncs$npc]
## visual checks
# ggplot(data=weighted_fpc[newid %in% c(1,2,5, 500, 999,1000)], 
#        aes(x=time, y=cd4_predicted, color=factor(newid)))+
#   geom_path()+
#   geom_point()+
#   geom_line(aes(y=cd4_ltfu, id=factor(newid)), color='black' )+
#   geom_point(aes(y=cd4_ltfu, id=factor(newid), shape=factor(newid)), color='black')


##f) lme
naive_lme<-tryCatch(
{
  #naive_lme_model<-lme(cd4 ~ bs(time, df=15), random=~1|newid, data=observed_with_stipw);
  #naive_lme <- data.frame(time=time_vec, V1=predict(naive_lme_model, newdata=data.frame(time=time_vec), level=0), approach="naive_lme")

  ## for 7*12 timepts and 1000 subjects, df=7 gives warning.  Stay at df=5.
  ## for 7*12 timepts and 5000 subjects, df=5 gives FATAL ERROR for ses_coef=0.05
  ## for 7*12 timepts and 5000 subjects, df=5    is OK          for ses_coef=0.01
    naive_lme_proc_minutes <- system.time(
  naive_lme_model<-lmer(cd4 ~ bs(time, df=5) + (bs(time, df=5)|newid), data=observed_with_stipw)
          )[3]/60
  ## re.form=~0
  #naive_lme <- data.frame(time=time_vec, V1=predict(naive_lme_model, newdata=data.frame(time=time_vec), re.form=~0), approach="naive_lme")
  ## re.form=~1
  naive_lme <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean         = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = predict(naive_lme_model,
                                                     newdata=all_with_stipw),
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="naive_lme",
                          minutes          = naive_lme_proc_minutes,
                          number_pc        = NA)
},
warning =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  ##naive_lme <- data.frame(time=time_vec, V1=NA, approach="naive_lme")  ;
  naive_lme <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean         = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="naive_lme",
                          minutes          = naive_lme_proc_minutes,
                          number_pc        = NA)
},
error =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  #naive_lme <- data.frame(time=time_vec, V1=NA, approach="naive_lme")  ;
  naive_lme <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean         = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="naive_lme",
                          minutes          = naive_lme_proc_minutes,
                          number_pc        = NA)
})
setkey(naive_lme, newid, time)
# quick checks:
# summary(naive_lme)
## visual checks
# ggplot(data=naive_lme[newid %in% c(1,2,5, 500, 999,1000)], 
#        aes(x=time, y=cd4_predicted, color=factor(newid)))+
#   geom_path()+
#   geom_point()+
#   geom_line(aes(y=cd4_ltfu, id=factor(newid)), color='black' )+
#   geom_point(aes(y=cd4_ltfu, id=factor(newid), shape=factor(newid)), color='black')


##f) weighted_remean_lme
wtd_remean_lme<-tryCatch(
{
  #naive_lme_model<-lme(cd4 ~ bs(time, df=15), random=~1|newid, data=observed_with_stipw);
  #naive_lme <- data.frame(time=time_vec, V1=predict(naive_lme_model, newdata=data.frame(time=time_vec), level=0), approach="naive_lme")



#   wtd_remean_lme_model<-lmer(cd4_wtd_remean ~ bs(time, df=15) + (1|newid), data=all_with_stipw,
#                              na.action=na.omit);
  ## for 7*12 timepts and 1000 subjects, df=7 gives warning.  Stay at df=5.
      wtd_remean_proc_minutes <- system.time(
  wtd_remean_lme_model<-lmer(cd4_wtd_remean ~ bs(time, df=5) + (bs(time, df=5)|newid), data=all_with_stipw,
                             na.action=na.omit)
            )[3]/60
  ## re.form=~0
  #naive_lme <- data.frame(time=time_vec, V1=predict(naive_lme_model, newdata=data.frame(time=time_vec), re.form=~0), approach="naive_lme")
  ## re.form=~1
  wtd_remean_lme <-
               data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean         = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = predict(wtd_remean_lme_model,
                                                     newdata=all_with_stipw),
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_remean_lme",
                          minutes = wtd_remean_proc_minutes,
                          number_pc = NA)
},
warning =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  ##naive_lme <- data.frame(time=time_vec, V1=NA, approach="naive_lme")  ;
  wtd_remean_lme <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean         = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_remean_lme",
                          minutes = wtd_remean_proc_minutes,
                          number_pc = NA)
},
error =function(cond){
  write.csv(observed_with_stipw, paste0("data_that_failed_nlme_lme_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  ## re.form=~0
  #naive_lme <- data.frame(time=time_vec, V1=NA, approach="naive_lme")  ;
  wtd_remean_lme <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean         = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_remean_lme",
                          minutes = wtd_remean_proc_minutes,
                          number_pc = NA)
})
setkey(wtd_remean_lme, newid, time)
## quick checks:
# summary(wtd_remean_lme)
## visual checks
# ggplot(data=wtd_remean_lme[newid %in% c(1,2,5, 500, 999,1000)], 
#        aes(x=time, y=cd4_predicted, color=factor(newid)))+
#   geom_path()+
#   geom_point()+
#   geom_line(aes(y=cd4_ltfu, id=factor(newid)), color='black' )+
#   geom_point(aes(y=cd4_ltfu, id=factor(newid), shape=factor(newid)), color='black')



## I have two wtd_lme chunks -- only have one uncommented at a time!
## the one immediately preceding this comment is lme(of wtd cd4);
## whereas the other one below it are lme4::lmer of cd4.
# wtd_lme<-tryCatch(
# {
#   wtd_lme_model<-lme(cd4_wtd_hadamard ~ bs(time, df=15), random=~1|newid, data=wtd_trajectories,
#                      na.action=na.omit);
#   ##predict(wtd_lme_model, newdata=data.frame(time=time_vec), level=0)
#   wtd_lme <- data.frame(time=time_vec, V1=predict(wtd_lme_model, newdata=data.frame(time=time_vec), level=0), approach="wtd_lme")
# },
# warning =function(cond){
#   wtd_lme <- data.frame(time=time_vec, V1=NA, approach="wtd_lme")  ;
#   wtd_lme
# },
# error =function(cond){
#   wtd_lme <- data.frame(time=time_vec, V1=NA, approach="wtd_lme")  ;
#   wtd_lme
# })
# summary(wtd_lme)
#

wtd_lme<-tryCatch(
{
  #   wtd_lme_model<-lme(cd4_wtd_hadamard ~ bs(time, df=15), random=~1|newid, data=wtd_trajectories,
  #                      na.action=na.omit);

#   wtd_lme_model<-lmer(cd4_wtd_hadamard ~ bs(time, df=15) + (1|newid), data=wtd_trajectories,
#                       na.action=na.omit)
  wtd_lme_proc_minutes <- system.time(
  wtd_lme_model<-lmer(cd4_wtd_hadamard ~ bs(time, df=5) + (bs(time, df=5)|newid), data=wtd_trajectories,
                      na.action=na.omit)
  )[3]/60
  ##predict(wtd_lme_model, newdata=data.frame(time=time_vec), level=0)
  ##wtd_lme <- data.frame(time=time_vec, V1=predict(wtd_lme_model, newdata=data.frame(time=time_vec), level=0), approach="wtd_lme")
  ## note: below, use re.form=~0 in lme4:predict is equivalent to level=0 in nlme:lme
  ## re.form=~0
  ##wtd_lme <- data.frame(time=time_vec, V1=predict(wtd_lme_model, newdata=data.frame(time=time_vec), re.form=~0), approach="wtd_lme")
  #wtd_lme_pop_mean <- data.frame(time=time_vec, V1=predict(wtd_lme_model, newdata=data.frame(time=time_vec), re.form=~0), approach="wtd_lme")

  ## re.form=~1
  wtd_lme1  <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean       = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = predict(wtd_lme_model,
                                                     newdata=all_with_stipw),
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_hadamard_lme")

  wtd_lme.test<- wts[wtd_lme1]
  ## see how many 0's
  wtd_lme.test[,table(round(stipw01.n,2))]
  ## change 0's to NA
  wtd_lme.test[stipw01.n==0, stipw01.n:=NA]
  wtd_lme.test[, cd4_predicted_weighted:= cd4_predicted]
  wtd_lme.test[!is.na(stipw01.n), cd4_predicted_deweighted:= cd4_predicted/stipw01.n]
  ## get the wtd_population_mean in there:
  ## skip this time ## wtd_lme.test[,wtd_pop_mean:=wtd_lme_pop_mean$V1,by=newid]
  ##for now, if don't have observed data there, we just imputed the weighted mean
  ## kinda lame, think on it.
  ## skip this time ## wtd_lme.test[is.na(stipw01.n), cd4_predicted_deweighted:= wtd_pop_mean, by=newid]
  ## end extra steps:
  setkey(wtd_lme.test, newid, time)#, cd4_predicted_deweighted)
  setnames(wtd_lme.test, "cd4_predicted", "cd4_predicted_old")
  setnames(wtd_lme.test, "cd4_predicted_deweighted", "cd4_predicted")
  wtd_lme <- subset(wtd_lme.test, select=selected)
  wtd_lme[,minutes:=wtd_lme_proc_minutes]
  wtd_lme[,number_pc:=NA]
  },
warning =function(cond){
  write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  #wtd_lme <- data.frame(time=time_vec, V1=cond, approach="wtd_lme")  ;
  wtd_lme   <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean       = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_hadamard_lme",
                          minutes=wtd_lme_proc_minutes,
                          number_pc=NA)
},
error =function(cond){
  write.csv(wtd_trajectories, paste0("data_that_failed_lme4_lmer_fit_",abs(rnorm(1,100,100)),".csv"), row.names=FALSE)
  #wtd_lme <- data.frame(time=time_vec, V1=cond, approach="wtd_lme")  ;
  wtd_lme   <- data.table(newid            = all_with_stipw$newid,
                          time              = all_with_stipw$time,
                          cd4           = all_with_stipw$cd4,
                          cd4_wtd_hadamard       = all_with_stipw$cd4_wtd_hadamard,
                          cd4_wtd_remean       = all_with_stipw$cd4_wtd_remean,
                          cd4_ltfu      = all_with_stipw$cd4_ltfu,
                          cd4_predicted = NA,
                          remaining        = all_with_stipw$remaining,
                          stipw            = all_with_stipw$stipw,
                          stipw01          = all_with_stipw$stipw01,
                          stipw01.n        = all_with_stipw$stipw01.n,
                          approach="wtd_hadamard_lme",
                          minutes=wtd_lme_proc_minutes,
                          number_pc=NA)
})
setkey(wtd_lme, newid, time)
## quick checks:
# summary(wtd_lme)
## visual checks
# ggplot(data=wtd_lme[newid %in% c(1,2,5, 500, 999,1000)], 
#        aes(x=time, y=cd4_predicted, color=factor(newid)))+
#   geom_path()+
#   geom_point()+
#   geom_line(aes(y=cd4_ltfu, id=factor(newid)), color='black' )+
#   geom_point(aes(y=cd4_ltfu, id=factor(newid), shape=factor(newid)), color='black')




#rbind it!!!!!!!!!!!!!



## plot each approach's mean
## means <- rbind(true_avg, naive_non_parm_avg, wtd_non_parm_avg, naive_fpc, naive_fpc_pabw, weighted_fpc, naive_lme, wtd_lme)
# note: we're simplifying to (w)fpca and (w)lme, so we rbind() fewer than line above (if you decided to add more approaches later,
#         make sure they have same format):
means <- rbind(naive_fpc, weighted_fpc, weighted_remean_fpc, naive_lme, wtd_lme, wtd_remean_lme)
means[,     newid := as.integer(newid)]
means[, remaining := as.integer(remaining)]
setkey(means, "newid","time")



didya<-dcast(means, newid+time ~ approach, value.var="cd4_predicted")
## after rbind multiple instances, use this to melt it
##melt.didya <- melt(didya, id.vars=c("newid","time"), variable="approach", value="cd4_predicted")


# ## make additions columnwise, not rbind-wise.  Save those GBs.  Can melt once in memory.
# key(naive_fpc)
# key(weighted_fpc)
# key(weighted_remean_fpc)
# key(naive_lme)
# key(wtd_lme)
# key(wtd_remean_lme)
#
base.select <-
            c("newid","time", 
              "cd4", "cd4_wtd_hadamard","cd4_wtd_remean", "cd4_ltfu",
#               "cd4_predicted",
               "remaining",
               "stipw",
               "stipw01",
               "stipw01.n"#,
#               "approach"
              )

base <- subset(naive_fpc, select=base.select)

base_means <- base[didya]




# from older files:
# means$approach<-factor(means$approach, levels=unique(means$approach))
# colnames(means)[colnames(means)=="V1"]<- "cd4"


## next three chunks deal with missingness, then we rbind it us with means....

## this is a useful chunk to check the range of
## probality of being censored (prob.cens)
#   ddply(censored, .(time), function(w) sum(w$instudy==1) )
#   melt.prob.cens=ddply(censored, .(newid,time), function(w) w$prob.cens )
#   dcast.prob.cens=dcast(melt.prob.cens, newid~time, value.var="V1")
#   apply(dcast.prob.cens, 2, function(w) round(range(w),2))
#   head(censored,18)

## do the following to get precent missing at time 18:
## overall:
dcast.wtd.trajectories<-dcast(calculate_wtd_trajectories_hh2(calculate_stipw_hh2(censored,"keep")), newid~time, value.var="stipw")
last_time<-as.character(time_vec[length(time_vec)])
percent.missing.at.time.18=sum(is.na(dcast.wtd.trajectories[last_time]))/length(unlist(dcast.wtd.trajectories[last_time]))
percent.missing = colSums(is.na(dcast.wtd.trajectories[,-1]))/length(unlist(dcast.wtd.trajectories[last_time]))
percent.missing
## below/above median linear predictor...(legacy code `medianSES`):
dcast.wtd.trajectories<-dcast(calculate_wtd_trajectories_hh2(calculate_stipw_hh2(censored,"keep")), 
                              newid+age+sex+race+hetero+msm+ivdu~time, value.var="stipw")

## calculate linear predictors
##whamcast[, linear.pred:=predict.glm(simplog, newdata=whamcast)]
## set medianSES to linear predictor median -- could be used in future.  For NOW,
## let's just make the split between ivdu and non-ivdu
#medianSES<-median(whamcast$linear.pred, na.rm=TRUE)
medianmsm<-0
## for last_time, historically "18" for BG:
subbie<-subset(dcast.wtd.trajectories, as.numeric(msm!="no") <= medianmsm, last_time)
percent.missing.at.time.18.below.median=sum(is.na(subbie))/nrow(subbie)
## across time, for those below median / msm==0
subbie<-subset(dcast.wtd.trajectories, as.numeric(msm!="no") <= medianmsm, select=c(-1,-2))
percent.missing.below.median=colSums(is.na(subbie))/nrow(subbie)
## for last_time, historically "18" for BG:
subbie<-subset(dcast.wtd.trajectories, as.numeric(msm=="yes")  > medianmsm, last_time)
percent.missing.at.time.18.above.median=sum(is.na(subbie))/nrow(subbie)
## across time, for those above median / msm==1
subbie<-subset(dcast.wtd.trajectories, as.numeric(msm=="yes")  > medianmsm, select=c(-1,-2))
percent.missing.above.median=colSums(is.na(subbie))/nrow(subbie)

## note: this cbind() works because every time is present for each dataset in `means`
##       and the only non-scalars are vectors that are same length as number of time-levels
## return:
results <- cbind(
        sim_seed                     = as.integer(sim_seed),
        sample_size                  = as.integer(sample_size),
        censoring_coef               = censoring_coef,
        perc_ltfu_18                 = percent.missing.at.time.18,
        percent_missing              = percent.missing,
        percent_missing_below_median = percent.missing.below.median[as.character(time_vec)],
        percent_missing_above_median = percent.missing.above.median[as.character(time_vec)],
        #means,
        base_means)
## note:  choose means OR base_means.  means is from rbind() above and multiplies rows by 6.



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
# oos_time=5
# oos_ind=as.numeric(names(d[,-1])) <= oos_time
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
# ##wtd_fncs <- dcast(data=wtd_trajectories, formula= newid~time, value.var="cd4_wtd_hadamard")
# fpca_wtd_fncs_pred_oos <- fpca.face(Y=as.matrix(wtd_fncs[,-1]), Y.pred=as.matrix(oos.early[,-1]), argvals=time_vec, knots=10)
#
# weighted_fpc_pred_oos <- data.frame(time=time_vec, V1=fpca_wtd_fncs_pred_oos$mu, approach="weighted_fpc_pred_oos")
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
#       "sim_time_coef", sim_time_coef,
#       sep="_")


idnum1<-abs(round(10000000*rnorm(1)))
midletter<-sample(letters)[1]
idnum2<-abs(round(10000000*rnorm(1)))
##saveRDS(results,paste0("results_",label,".RDS" ))
saveRDS(results,paste0("results_",idnum1, midletter,idnum2, ".RDS"))



setkey(means, approach)
proc_minutes_number_pcs<-cbind(subset(unique(means),
                                      select=c("approach", "minutes","number_pc")),
                               sample_size=unique(results$sample_size),
                               censoring_coef=unique(results$censoring_coef)
)


write.csv(proc_minutes_number_pcs,
        paste0("proc_minutes_number_pcs_",idnum1, midletter,idnum2, ".csv"),
        row.names=FALSE)


}
