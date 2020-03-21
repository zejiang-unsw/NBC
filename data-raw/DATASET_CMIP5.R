################################################################
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

setwd("~/OneDrive - UNSW/R_Package/NBC/data-raw")

library(ncdf4)
library(raster)
library(rgdal)

#---------------------------------------------------------------------------------------------
flag.save = 0

#import filenames
path <- "~/OneDrive - UNSW/PhD Research Work/5th Paper/Real Data Case - Downscaling/"
GCM <- switch(2, "CanESM2", "CSIRO-Mk3-6-0")

ncfile1 <- paste0(path,"va_Amon_",GCM,"_historical_r1i1p1_185001-200512.nc")
ncfile2 <- paste0(path,"ua_Amon_",GCM,"_historical_r1i1p1_185001-200512.nc")
ncfile3 <- paste0(path,"ta_Amon_",GCM,"_historical_r1i1p1_185001-200512.nc")
ncfile4 <- paste0(path,"zg_Amon_",GCM,"_historical_r1i1p1_185001-200512.nc")
ncfile5 <- paste0(path,"hus_Amon_",GCM,"_historical_r1i1p1_185001-200512.nc")
ncfile6 <- paste0(path,"psl_Amon_",GCM,"_historical_r1i1p1_185001-200512.nc")

#---------------------------------------------------------------------------------------------
#land mask
map <- "Global_110m_land.shp"
land.mask <- readOGR(map)

#resample grids
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
date.char <- seq(as.Date("1850-01-01"), by="month", length.out = length(nc_time))

nc_close(nc_tmp)

# ncfile = ncfile5; var="hus"; level=6
# nc2ras <- sapply(level, function(i) rotate(raster(ncfile,varname = var, level=i)))
# image(nc2ras[[1]])
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
start="1850-01-01"; by ="month"

data.GCM.vwnd <- data_cmip(ncfile1, var="va", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.uwnd <- data_cmip(ncfile2, var="ua", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.air  <- data_cmip(ncfile3, var="ta", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.hgt  <- data_cmip(ncfile4, var="zg", level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.shum <- data_cmip(ncfile5, var="hus",level, start, by, grid=ncep.land.grid, mask=land.mask)
data.GCM.slp  <- data_cmip(ncfile6, var="psl",level=1, start, by, grid=ncep.land.grid, mask=land.mask)

lat_lon.cmip5 <- data.GCM.slp$xy
summary(lat_lon.cmip5)

#------------------------------------------------
data.list <- c("data.GCM.vwnd","data.GCM.uwnd","data.GCM.air","data.GCM.hgt","data.GCM.shum","data.GCM.slp")
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
  save(data.GCM.vwnd, file=paste0(GCM,".vwnd.Rdata"))
  save(data.GCM.uwnd, file=paste0(GCM,".uwnd.Rdata"))
  save(data.GCM.air,  file=paste0(GCM,".air.Rdata"))
  save(data.GCM.hgt,  file=paste0(GCM,".hgt.Rdata"))
  save(data.GCM.shum, file=paste0(GCM,".shum.Rdata"))
  save(data.GCM.slp,  file=paste0(GCM,".slp.Rdata"))

}


