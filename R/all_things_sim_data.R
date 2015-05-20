#' all_things_sim_data -- samples, censors, stipw's -- gets outsample_sim() up and running quickly
#'
#' not implemented in insample_sim() (yet).  Same inputs as in/outsample_sim()
#' 
#' @param sim_seed passed to set.seed()
#' @param sample_size (defaults to 1000) and must be a multiple of 100
#' @param sim_slope see slope in calculate_ses()
#' @param sim_intercept see intercept in calculate_ses()
#' @param sim_ses_coef see ses_coef in apply_censoring()
#' @param sim_age_coef see age_coef in apply_censoring()
#' @export
#' @return return(list(a =   wtd_trajectories, b = unwtd_fncs, c =   wtd_fncs, d = observed_with_stipw, e =      all_with_stipw, f = censored, g = age_vec))
#' @examples
#' ---
all_things_sim_data <- function(sim_seed, sample_size, sim_slope,
                                sim_intercept, sim_ses_coef, sim_age_coef){
  
  
  
  
  
  
  
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
  
  ## for fpca methods
  unwtd_fncs      <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches")
  wtd_fncs        <- dcast(data=wtd_trajectories, formula= newid~age, value.var="inches_wtd")
  
  
  setDT(observed_with_stipw)
  setDT(censored)
  
  return(list(
    a =   wtd_trajectories,
    b = unwtd_fncs,
    c =   wtd_fncs,
    d = observed_with_stipw,
    e =      all_with_stipw,
    f = censored,
    g = age_vec
  ))
}