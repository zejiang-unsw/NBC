## code to prepare `DATASET` dataset goes here
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

library(zoo)
library(devtools)

GCM <- switch(1, "CanESM2", "CSIRO-Mk3-6-0")
types = c("GCM","GCM.RCP")
vars = c("uwnd","vwnd","air","shum","hgt","slp")

setwd("C:/Users/z5154277/OneDrive - UNSW/R_Package/NBC/data-raw/")
flag.save = 1
flag.ncep = 0

#---------------------------------------------------
load("NCEP.vwnd.Rdata")
load("NCEP.uwnd.Rdata")
load("NCEP.air.Rdata")
load("NCEP.hgt.Rdata")
load("NCEP.shum.Rdata")
load("NCEP.slp.Rdata")

load("Ind_NCEP.Rdata")
is.na(data.NCEP.slp$data)

setwd(paste0("./CMIP5_",GCM,"/"))
load(paste0(GCM,".vwnd.Rdata"))
load(paste0(GCM,".uwnd.Rdata"))
load(paste0(GCM,".air.Rdata"))
load(paste0(GCM,".hgt.Rdata"))
load(paste0(GCM,".shum.Rdata"))
load(paste0(GCM,".slp.Rdata"))

load(paste0(GCM,".RCP.vwnd.Rdata"))
load(paste0(GCM,".RCP.uwnd.Rdata"))
load(paste0(GCM,".RCP.air.Rdata"))
load(paste0(GCM,".RCP.hgt.Rdata"))
load(paste0(GCM,".RCP.shum.Rdata"))
load(paste0(GCM,".RCP.slp.Rdata"))

summary(data.NCEP.slp)
summary(data.GCM.slp)
summary(data.GCM.RCP.slp)

#---------------------------------------------------
#subset:
# study period: current - 1951-2005; future - 2021-2075
# study area:
# start.yr=c(1951,1); end.yr=c(2005,12)
# start.yr.f=c(2021,1); end.yr.f=c(2075,12)
start.yr="1951-01-01"; end.yr="2005-12-01"
start.yr.f="2021-01-01"; end.yr.f="2075-12-01"

###NCEP: 1951-2005
for(var in vars){

  data.list <- eval(parse(text=paste0("data.NCEP.",var)))

  data <- data.list$data; date <- data.list$Date
  date.n <- seq(as.Date(start.yr), as.Date(end.yr), by="month")
  index <- which(date %in% date.n)

  data.n <- lapply(data, function(df) df[index,])

  data.list$data <- data.n
  data.list$Date <- date.n

  assign(paste0("data.NCEP.",var),data.list)

}

###GCM current: 1951-2005
for(var in vars){

  data.list <- eval(parse(text=paste0("data.GCM.",var)))

  data <- data.list$data; date <- data.list$Date
  date.n <- seq(as.Date(start.yr), as.Date(end.yr), by="month")
  index <- which(date %in% date.n)

  data.n <- lapply(data, function(df) df[index,])

  data.list$data <- data.n
  data.list$Date <- date.n

  assign(paste0("data.GCM.",var),data.list)

}

###GCM future: 2021-2075
for(var in vars){

  data.list <- eval(parse(text=paste0("data.GCM.RCP.",var)))

  data <- data.list$data; date <- data.list$Date
  date.n <- seq(as.Date(start.yr.f), as.Date(end.yr.f), by="month")
  index <- which(date %in% date.n)

  data.n <- lapply(data, function(df) df[index,])

  data.list$data <- data.n
  data.list$Date <- date.n

  assign(paste0("data.GCM.RCP.",var),data.list)

}

#---------------------------------------------------
#overview at sample grids
#1659; 1685; 1131 - CSIRO   2149 - CanESM2
sample = sample(Ind_NCEP,1)
sample = switch(GCM,"CanESM2"=2149,"CSIRO-Mk3-6-0"=1131)
for(var in vars){ #variables

  df1 <- eval(parse(text=paste0("data.NCEP.",var)))$data
  df2 <- eval(parse(text=paste0("data.GCM.",var)))$data
  df3 <- eval(parse(text=paste0("data.GCM.RCP.",var)))$data

  if(TRUE){
      jpeg(filename = paste0(GCM,"_",var,"_Sample_G",sample,".jpg"),
           width = 115, height = 95, units = "mm", res = 600)

      par(mfrow=c(length(df1),1), mar=c(2,3,1,2), mgp=c(1, 0.5, 0),
          bg = "transparent", ps=8)
      ts.plot(ts(cbind(df1[[1]][,sample],df2[[1]][,sample],df3[[1]][,sample])),
              col=c("black","blue","red"), ylab=var, main=paste0("Sample Grid: ",sample),
              xlab=NA)
      legend("top", inset=c(0,-0.02), box.lty=0, bg=NA, ncol=3,
             legend=c("Observed","GCM current","GCM future"),
             lty=c(1,1,1),col=c("black","blue","red"))

      if(length(df1)>1){
          par(mar=c(2,3,1,2))
          for(i in 2:length(df1)){ #levels
              ts.plot(ts(cbind(df1[[i]][,sample],df2[[i]][,sample],df3[[i]][,sample])),
                      col=c("black","blue","red"), ylab=var,
                      xlab=NA)
          }
      }

      dev.off()
  }

}


#---------------------------------------------------
#save to data folder
if(flag.save){
  #GCM current and future data
  for(type in types){
    for(var in vars){
      name = paste0("data.",type,".",var)
      do.call("use_data",list(as.name(name),overwrite=T))
    }
  }

  if(flag.ncep){ #NCEP current data
      for(var in vars){
          name = paste0("data.NCEP.",var)
          do.call("use_data",list(as.name(name),overwrite=T))
      }
  }


}










