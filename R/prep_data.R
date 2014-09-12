#' Prep the Berkeley Growth Study Data
#'
#' This returns "boys".
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' d<-prep_data()
#' head(d)
prep_data <- function(){
  library(fda)
  library(reshape2)
  ## @knitr boysAndGirls
  ## in Functional Data Analysis, matplot() is used often for convenience.
  ## we'll opt to "melt" datasets and make nice(r) graphics in ggplot
  ## 'growth' is the data object from fda and contains Berkeley Study Data
  ##with(growth, matplot(age, hgtf, hgtm, type=c("l","l"), col=c("pink","blue")))
  
  ## @knitr ggDataPrep
  ## get into ggplot2 data format.
  ## there are probably nicer/more succinct ways to do this, but the following works.
  ## take 'growth' and make 'gg.growth' which takes matrices, takes out the columns, and stacks said columns
  ## we also add an $id variable, as well as an $inches variable
  (total.sub <- with(growth,ncol(hgtf)+ncol(hgtm)))
  gg.growth <- with(growth, data.frame(age=rep(age, ncol(hgtf)+ncol(hgtm)),
                                       hgt=c(hgtf,hgtm),
                                       sex=c(rep("female", length(age)*ncol(hgtf)),
                                             rep("male"  , length(age)*ncol(hgtm)))))
  gg.growth$id <- with(growth,rep(1:(ncol(hgtf)+ncol(hgtm)), each=length(age)))
  gg.growth$inches <- gg.growth$hgt / 2.54
  
  ## @knitr aFew
  ## here we plot a few from gg.growth to see the ggplot style graphics
  ##ggplot(subset(gg.growth, id %in% c(1,40,88,93 )), aes(x=age, y=hgt, colour=sex, group=id), alpha=.5) + geom_line()
  
  ## @knitr all
  ## we can assign ggplot() to objects and then add options to that object for display.
  ## here we assign to data.view (no display), and then we display data.view with an overlay of group averages by loess
  ##data.view <- ggplot(gg.growth, aes(x=age, y=inches, colour=sex, group=id)) + geom_line(alpha=.25) + geom_rug(sides="b",color="black")
  #data.view+geom_smooth(aes(group=sex), method = "loess", size = 1.5, alpha=1)
  ##data.view+geom_smooth(aes(group=sex), method = "loess", size = 1.5, alpha=1,span=.2)
  
  
  ## @knitr extractBoys
  ## in our protocol we isolate the boys.  This can be done with dcast,
  ## which takes the "melted" gg.growth for only males and makes it 1 row to 1 id
  ## and the age over the columns (so we cast it into the matrix, 'boys',  as in the traditional
  ## Functional Data Analysis way of doing things).  dcast() is in library(reshape2)
  boys <- dcast(subset(gg.growth, sex=="male"), id~age)
  #boys[1:6,1:6]
  boys
}