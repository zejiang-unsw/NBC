context("Nested Bias Correction")
library(scales)
#--------------------------------------------
test_that("NBC", {
  ###load raw data
  for(j in 1:2){
  variable= switch(j, "uwnd","vwnd","air","shum","hgt","slp")
  i=1 #pressure level
  GCM <- switch(1, "CanESM2.", "CSIRO-Mk3-6-0.")

  load(system.file("extdata", "land.mask.rda", package="NBC"))
  load(system.file("extdata", paste0("NCEP.",variable,".rda"), package="NBC"))
  load(system.file("extdata", paste0(GCM,variable,".rda"), package="NBC"))
  load(system.file("extdata", paste0(GCM,"RCP.",variable,".rda"), package="NBC"))

  data.NCEP <- eval(parse(text=paste0("data.NCEP.",variable)))
  NCEP.var <- data.NCEP$data[[i]]
  NCEP.grid<- data.NCEP$xy

  data.GCM <- eval(parse(text=paste0("data.GCM.",variable)))
  GCM.var <- data.GCM$data[[i]]

  data.GCM.RCP <- eval(parse(text=paste0("data.GCM.RCP.",variable)))
  GCM.RCP.var <- data.GCM.RCP$data[[i]]

  Date <- format(data.NCEP$Date, "%Y %m")
  Date.f <- format(data.GCM.RCP$Date, "%Y %m")
  start.yr <- 1951; end.yr <- 2005
  start.yr.f <- 2021; end.yr.f <- 2075

  ###overview of bias correction at a sampled grid NBC
  samples=sample(1:nrow(NCEP.grid), 4)
  samples=c(2149,1255,2135,1937)
  samples=c(2250,2392,2410,2456)
  grid.xy=NCEP.grid[samples,]

  # #overview - study grid
  # par(mfrow=c(1,1),pty="m")
  # plot(land.mask);points(grid.xy, col="red", pch=16)
  # text(grid.xy, labels=rownames(grid.xy), pos=3, col="red", pch=16)

  #------------------------------------------------------------------------------
  ###compute statistics
  #NCEP current
  obs.cali <- ts(NCEP.var[,samples], start=c(start.yr,1), frequency = 12)
  obs.sta <- obs.stats(n=ncol(obs.cali),obs.cali, fun="mean")

  #GCM current
  mod.cali <- ts(GCM.var[,samples], start=c(start.yr,1), frequency = 12)
  mod.sta <- obs.stats(n=ncol(mod.cali),mod.cali,fun="mean")

  #GCM future
  mod.vali <- ts(GCM.RCP.var[,samples], start=c(start.yr.f,1), frequency = 12)
  mod.sta.f <- obs.stats(n=ncol(mod.vali),mod.vali,fun="mean")

  #------------------------------------------------------------------------------
  ###NBC
  #GCM current
  mod.bcc <- nest.mod.cal(n=ncol(mod.cali),mod.cali, obs.sta,
                          start.yr, end.yr, fun="mean", tol=0.1)

  #GCM future
  mod.bcf <- nest.mod.val(n=ncol(mod.vali),mod.vali, obs.sta, mod.sta,
                          start.yr.f, end.yr.f, fun="mean", tol=0.1)

  #----------------------------------
  par(mfrow=c(4,2), mar=c(0,3,2,1),# margin of the plot
      oma = c(2, 1, 1, 2), # move plot to the right and up
      mgp=c(1, 0.5, 0), # move axis labels closer to axis
      bg = "transparent", pty="m", # maximal plotting region
      ps=8)
  for(grid in seq_along(samples)){
    #if(TRUE){ ###NBC


    ts.plot(obs.cali[,grid], mod.cali[,grid], mod.bcc$gcm.cor[,grid],
            gpars=list(xaxs="i", xlab=NA, ylab=paste0(variable,": G",grid),
                       xlim=c(start.yr-1,end.yr+1),
                       col=c(alpha("black",0.6),alpha("blue",1),alpha("red",1)),
                       lwd=c(2,1,1)))
    legend("bottom", inset=c(0,-0.02), box.lty=0, bg=NA, ncol=3,
           legend=c("Observed","GCM sim","GCM NBC sim"),
           lty=c(1,1,1),col=c("black","blue","red"))

    ts.plot(mod.vali[,grid], mod.bcf$gcm.cor[,grid],
            gpars=list(xaxs="i", xlab=NA, ylab=paste0(variable,": G",grid),
                       xlim=c(start.yr.f-1,end.yr.f+1),
                       col=c(alpha("blue",1),alpha("red",1)),
                       lwd=c(1,1)))

    legend("bottom", inset=c(0,-0.02), box.lty=0, bg=NA, ncol=3,
           legend=c(NA,"GCM sim","GCM NBC sim"),
           lty=c(NA,1,1),col=c(NA,"blue","red"))

  }
  }


})
