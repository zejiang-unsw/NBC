################################################################
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

setwd("~/OneDrive - UNSW/PhD Research Work/5th Paper/WRR_2020/data-raw")

library(ncdf4)
library(raster)
library(rgdal)
source("data_NC.R")

#---------------------------------------------------------------------------------------------
flag.save = 1

#import filenames 
path <- "../../Real Data Case - Downscaling/"
GCM <- switch(2, "CanESM2", "CSIRO-Mk3-6-0")

ncfile1 <- paste0(path,"va_Amon_",GCM,"_rcp85_r1i1p1_200601-210012.nc")
ncfile2 <- paste0(path,"ua_Amon_",GCM,"_rcp85_r1i1p1_200601-210012.nc")
ncfile3 <- paste0(path,"ta_Amon_",GCM,"_rcp85_r1i1p1_200601-210012.nc")
ncfile4 <- paste0(path,"zg_Amon_",GCM,"_rcp85_r1i1p1_200601-210012.nc")
ncfile5 <- paste0(path,"hus_Amon_",GCM,"_rcp85_r1i1p1_200601-210012.nc")
ncfile6 <- paste0(path,"psl_Amon_",GCM,"_rcp85_r1i1p1_200601-210012.nc")

#---------------------------------------------------------------------------------------------
#import basemap
map <- "Global_110m_land.shp"
land.mask <- readOGR(map)

ncep.land <- "Global_Land_NCEP_Grid.shp"
ncep.land.grid <- readOGR(ncep.land)

plot(land.mask); plot(ncep.land.grid, add=T)

#---------------------------------------------------------------------------------------------
#read ncfile 
###overview of the ncfile
nc_tmp <- nc_open(ncfile5)
attributes(nc_tmp)
var <- attributes(nc_tmp$var)$names;var
attributes(nc_tmp$dim)

nc_lon <- ncvar_get(nc_tmp, "lon")
nc_lat <- ncvar_get(nc_tmp, "lat")
nc_plev <- ncvar_get(nc_tmp, "plev");nc_plev #unit Pa
level <- which(as.vector(nc_plev/100) %in% c(500,700,850)) #index of pressure at 500, 700, and 850 hPa

nc_time <- ncvar_get(nc_tmp, "time")
tunits <- ncatt_get(nc_tmp,"time","units");tunits
date.char <- seq(as.Date("2006-01-01"), by="month", length.out = length(nc_time))
head(date.char); tail(date.char)
nc_close(nc_tmp)

# ncfile = ncfile6; var="psl"; level=1
# nc2ras <- sapply(level, function(i) rotate(brick(ncfile,varname = var, level=i)))
# nc2ras
# 
# df.ras.GCM <- raster(SpatialPixelsDataFrame(ncep.land.grid, ncep.land.grid@data))
# 
# df.ras.GCM <- lapply(nc2ras, function(ras) mask(resample(ras, df.ras.GCM, method="bilinear"),land.mask))
# df.ras.GCM
# image(df.ras.GCM[[1]])
# 
# df.list <- lapply(df.ras.GCM, function(ls) raster::as.data.frame(ls, xy = TRUE, na.rm=T))
# data.GCM <- lapply(df.list, function(df) matrix(t(subset(df, select=-c(x,y))), ncol=nrow(df)))
# 
# lat_lon <- data.frame(longitude=df.list[[1]]$x, latitude=df.list[[1]]$y)
# summary(lat_lon)

#------------------------------------------------
start="2006-01-01"; by ="month"

data.GCM.RCP.vwnd <- data_cmip(ncfile1, var="va", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.RCP.uwnd <- data_cmip(ncfile2, var="ua", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.RCP.air  <- data_cmip(ncfile3, var="ta", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.RCP.hgt  <- data_cmip(ncfile4, var="zg", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.RCP.shum <- data_cmip(ncfile5, var="hus",level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.RCP.slp  <- data_cmip(ncfile6, var="psl",level=1, start, by, grid=ncep.land.grid, mask=land.mask)

lat_lon.cmip5 <- data.GCM.RCP.slp$xy
summary(lat_lon.cmip5)

#------------------------------------------------
data.list <- c("data.GCM.RCP.vwnd","data.GCM.RCP.uwnd","data.GCM.RCP.air",
               "data.GCM.RCP.hgt","data.GCM.RCP.shum","data.GCM.RCP.slp")
for(data in data.list){
  # identify NaN grid
  data.GCM <- eval(parse(text=data))
  Ind_GCM <- which(apply(data.GCM$data[[1]], 2, function(m) sum(is.na(m)))==0)
  Ind_GCM_NaN <- which(apply(data.GCM$data[[1]], 2, function(m) sum(is.na(m)))==nrow(data.GCM$data[[1]]))
  length(Ind_GCM)+length(Ind_GCM_NaN)
  
  print(sum(!is.na(data.GCM$data[[1]][,Ind_GCM[1]])))
  print(sum(!is.na(data.GCM$data[[1]][,Ind_GCM[length(Ind_GCM)]])))
  
}

#------------------------------------------------
#load and compare
load("lat-lon_NCEP.Rdata")
load("data.NCEP.slp.Rdata")
load("Ind_NCEP.Rdata")
# summary(lat_lon)
unique(lat_lon$longitude);unique(lat_lon$latitude)

sum(abs(Ind_GCM-Ind_NCEP))
unique(lat_lon.cmip5$longitude)-unique(lat_lon$longitude)
unique(lat_lon.cmip5$latitude)-unique(lat_lon$latitude)

#------------------------------------------------
if(flag.save){
  save(data.GCM.RCP.vwnd, file=paste0(GCM,".RCP.vwnd.Rdata"))
  save(data.GCM.RCP.uwnd, file=paste0(GCM,".RCP.uwnd.Rdata"))
  save(data.GCM.RCP.air,  file=paste0(GCM,".RCP.air.Rdata"))
  save(data.GCM.RCP.hgt,  file=paste0(GCM,".RCP.hgt.Rdata"))
  save(data.GCM.RCP.shum, file=paste0(GCM,".RCP.shum.Rdata"))
  save(data.GCM.RCP.slp,  file=paste0(GCM,".RCP.slp.Rdata"))
  
}


