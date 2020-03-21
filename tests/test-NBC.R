######################################################################
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

### Set working directory
setwd("~/OneDrive - UNSW/R_Package/NBC/tests")

library(NBC)
#------------------------------------------------------------------------------
flag.plt = 1

iteration = 3 # RNBC iterations
start.yr <- 1951; end.yr <- 2005
start.yr.f <- 2021; end.yr.f <- 2075

#----------------------------------
#import basemap
map <- "../data-raw/Global_110m_land.shp"
land.mask <- readOGR(map)

###load raw data
variable="slp"; i=1
data(NCEP.slp)
data(CanESM2.slp)
data(CanESM2.RCP.slp)
Date <- format(data.NCEP.slp$Date, "%Y %m")
Date.f <- format(data.GCM.RCP.slp$Date, "%Y %m")

data.NCEP <- eval(parse(text=paste0("data.NCEP.",variable)))
NCEP.var <- data.NCEP$data[[i]]
NCEP.grid<- data.NCEP$xy

data.GCM <- eval(parse(text=paste0("data.GCM.",variable)))
GCM.var <- data.GCM$data[[i]]

data.GCM.RCP <- eval(parse(text=paste0("data.GCM.RCP.",variable)))
GCM.RCP.var <- data.GCM.RCP$data[[i]]
#------------------------------------------------------------------------------
###compute statistics
#NCEP current
obs.cali <- ts(NCEP.var, start=c(start.yr,1), frequency = 12)
obs.sta <- obs.stats(n=ncol(obs.cali),obs.cali, fun="mean")

#GCM current
mod.cali <- ts(GCM.var, start=c(start.yr,1), frequency = 12)
mod.sta <- obs.stats(n=ncol(mod.cali),mod.cali,fun="mean")

#GCM future
mod.vali <- ts(GCM.RCP.var, start=c(start.yr.f,1), frequency = 12)
mod.sta.f <- obs.stats(n=ncol(mod.vali),mod.vali,fun="mean")

#------------------------------------------------------------------------------
###NBC
#GCM current
mod.bcc <- nest.mod.cal(n=ncol(mod.cali),mod.cali, obs.sta,
                        start.yr, end.yr, fun="mean", tol=0.1)

#GCM future
mod.bcf <- nest.mod.val(n=ncol(mod.vali),mod.vali, obs.sta, mod.sta,
                        start.yr.f, end.yr.f, fun="mean", tol=0.1)

#------------------------------------------------------------------------------
###RNBC
mod.bcc.r <- rnbc.mod.cal(n=ncol(mod.cali),mod.cali, obs.sta,
                          start.yr, end.yr, fun="mean", tol=0.1, iteration)
sum(abs(mod.bcc$gcm.cor-mod.bcc.r))

# mod.bcf.r <- rnbc.mod.val(n=ncol(mod.vali),mod.vali, obs.sta, mod.sta,
#                           start.yr.f, end.yr.f, fun="mean", tol=0.1, iteration)
# sum(abs(mod.bcf$gcm.cor-mod.bcf.r))
#------------------------------------------------------------------------------
###overview of bias correction at a sampled grid NBC
grid=sample(1:nrow(NCEP.grid), 1)
grid=2149   # 1937; 1255; 2135; 2149
grid.xy=NCEP.grid[grid,]

#overview - study grid
par(mfrow=c(1,1),pty="m")
plot(land.mask);points(grid.xy, col="red", pch=16)

#----------------------------------
if(TRUE){ ###NBC
par(mfcol=c(2,1), mar=c(0,3,2,1),# margin of the plot
    oma = c(2, 1, 1, 2), # move plot to the right and up
    mgp=c(1, 0.5, 0), # move axis labels closer to axis
    bg = "transparent", pty="m", # maximal plotting region
    ps=8)

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

p1 <- plot_grid(recordPlot())

}
# #----------------------------------
# if(TRUE){ ###RNBC
#   par(mfcol=c(2,1), mar=c(0,3,2,1),# margin of the plot
#       oma = c(2, 1, 1, 2), # move plot to the right and up
#       mgp=c(1, 0.5, 0), # move axis labels closer to axis
#       bg = "transparent", pty="m", # maximal plotting region
#       ps=8)
#
#   ts.plot(obs.cali[,grid], mod.cali[,grid], mod.bcc.r[,grid],
#           gpars=list(xaxs="i", xlab=NA, ylab=paste0(variable,": G",grid),
#                      xlim=c(start.yr-1,end.yr+1),
#                      col=c(alpha("black",0.6),alpha("blue",1),alpha("red",1)),
#                      lwd=c(2,1,1)))
#
#   ts.plot(mod.vali[,grid], mod.bcf.r[,grid],
#           gpars=list(xaxs="i", xlab=NA, ylab=paste0(variable,": G",grid),
#                      xlim=c(start.yr.f-1,end.yr.f+1),
#                      col=c(alpha("blue",1),alpha("red",1)),
#                      lwd=c(1,1)))
#
#   p2 <- plot_grid(recordPlot())
#
# }

#------------------------------------------------------------------------------
if(flag.plt){

  ggsave(paste0("Grid",grid,"_",variable,"_NBC.jpg"),
         width = 115, height = 115, units = "mm", dpi = 600, p1)

  # ggsave(paste0("Grid",grid,"_",variable,"_RNBC.jpg"),
  #        width = 115, height = 115, units = "mm", dpi = 600, p2)

}



