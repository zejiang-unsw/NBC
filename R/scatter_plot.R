#' Scatter Plot for statistics check
#'
#' @param obs.sta observed statistics
#' @param mod.sta modelled statistics
#'
#' @return
#' @export
#'
#' @examples
scatterPlot <- function(obs.sta, mod.sta){

  type <- c("Means","Standard Deviations","LAG1 correlations")
  par(mfcol=c(3,1),
      mar=c(3,4,2,2),# margin of the plot
      oma = c(2, 1, 1, 2), # move plot to the right and up
      mgp=c(2, 0.5, 0), # move axis labels closer to axis
      bg = "transparent",
      pty="m", # maximal plotting region
      cex.lab=2)

  for(i in 1:3){
    #i <- 1
    #monthly statistics
    x1 <- obs.sta[[i+6]]
    y1 <- mod.sta[[i+6]]

    #annual statistics
    x2 <- obs.sta[[i+3]]
    y2 <- mod.sta[[i+3]]


    all = c(x1,y1,x2,y2)
    range = c(min(all), max(all))
    plot(x1,y1,lty=1,xlim=range, ylim=range, lwd=1,pch=0,col=c("red"), xlab= "Observed", ylab="Simulated")
    points(x2,y2, lty=1, lwd=1,pch=1, col=c("blue"));abline (0, 1)

    legend("topleft",title=type[i],cex=1.5,
           legend=c("Monthly","Annual"),pch=c(0,1),col=c("red","blue"),bg = "transparent",box.lty=0)


  }

  return(recordPlot())

}
