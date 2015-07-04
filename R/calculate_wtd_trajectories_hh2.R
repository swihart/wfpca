#' Calculate the ST-STIPW based on for graphical viewing as well as input to fpca methods
#'
#' Takes a calculate_stipw(,na.action="omit") object and calculates the weighted trajectories
#' where the weights are standardized (within observation time) STIPWs so that
#' the (regular, unweighted) mean of the weighted trajectories is the same as the 
#' weighted mean of the original trajectories.
#' Names are hard coded, so stick
#' to script in the examples.
#' @param data_in an object returned from calculate_stipw()
#' @export
#' @return a dataset similar to data_in with stipw in columns; and potentially fewer rows or NA-filled rows
#' for induced censoring.  
#' @examples
#' d<-prep_data()
#' head(d)
#' over_samp_mat<-sample_data(d,1000)
#' with_ses <- calculate_ses(over_samp_mat)
#' long <-make_long(with_ses)
#' head(long)
#' censored <- apply_censoring(long)
#' head(censored,18)
#' observed_with_stipw <- calculate_stipw(censored,"omit")
#' wtd_trajectories <- calculate_wtd_trajectories(observed_with_stipw)
#' head(wtd_trajectories)
calculate_wtd_trajectories_hh2 <- function(data_in=NULL){
## BJS developed this in externalized.R while
## BL  was creating externalized_BL2.R.
## my hope is that it carries over well to visualizing what the wts are doing
## to the curves now that we have better weights thanks to BL edits.
## @knitr std.wts.0.1
## Attempt: to see what weighted curves look like before mean.  Issue:
## To do a weighted average the denominator is the sum of all weights,
## so intuition is out.  What if: I stdize all weights to add to one.
## Then multiply by total number of curves.
head(data_in)
denom <- ddply(data_in, .(time), function(w) sum(w$stipw))
names(denom)[2] <- "denom"
remaining <- ddply(data_in, .(time), function(w) sum(w$instudy))
names(remaining)[2] <- "remaining"
cbind(remaining[,2], denom[,2])
long01.int <- join(data_in,denom, by="time")
long01 <- join(long01.int, remaining, by="time")

# long01[,stipw01 := stipw/denom]
# long01[,stipw01.n := remaining*stipw01]
# long01[,inches_wtd_hadamard := cd4*stipw01.n]
# long01[inches_wtd_hadamard==0,inches_wtd_hadamard:= NA]  ## uncomment; 

long01$stipw01 <- long01$stipw/long01$denom
long01$stipw01.n <- long01$remaining*long01$stipw01
long01$cd4_wtd_hadamard <- long01$cd4*long01$stipw01.n
long01$cd4_wtd_hadamard[long01$cd4_wtd_hadamard==0] <- NA  ## uncomment; 

## calculate as you did above for wtd_avg; avg, avg_bias
##avg_biased01 <- ddply(subset(long01,instudy==1), .(time), function(w)  mean(w$inches_wtd))
##names(avg_biased01)[2] <- "ab01"
long01
}