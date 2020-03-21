## code to prepare `DATASET` dataset goes here

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

ncfile1 <- paste0(path,"vwnd.mon.mean.nc")
ncfile2 <- paste0(path,"uwnd.mon.mean.nc")
ncfile3 <- paste0(path,"air.mon.mean.nc")
ncfile4 <- paste0(path,"hgt.mon.mean.nc")
ncfile5 <- paste0(path,"shum.mon.mean.nc")
ncfile6 <- paste0(path,"slp.mon.mean.nc")

#---------------------------------------------------------------------------------------------
#import basemap
map <- "Global_110m_land.shp"
land.mask <- readOGR(map)

plot(land.mask)
#---------------------------------------------------------------------------------------------
#read ncfile
###overview of the ncfile
nc_tmp <- nc_open(ncfile5)
attributes(nc_tmp)
var <- attributes(nc_tmp$var)$names;var
attributes(nc_tmp$dim)

nc_lon <- ncvar_get(nc_tmp, "lon")
nc_lat <- ncvar_get(nc_tmp, "lat")
nc_plev <- ncvar_get(nc_tmp, "level");nc_plev #unit hPa
level <- which(as.vector(nc_plev) %in% c(500,700,850)) #index of pressure at 500, 700, and 850 hPa

nc_time <- ncvar_get(nc_tmp, "time")
tunits <- ncatt_get(nc_tmp,"time","units");tunits
date.char <- seq(as.Date("1949-01-01"), by="month", length.out = length(nc_time))

nc_close(nc_tmp)

# ncfile=ncfile6; var="slp"; level=1
# nc2ras <- sapply(level, function(i) mask(rotate(brick(ncfile, varname = var, level=i)),land.mask))
# nc2ras
# 
# df.list <- lapply(nc2ras, function(ls) raster::as.data.frame(ls, xy = TRUE, na.rm=T))
# data.NCEP <- lapply(df.list, function(df) matrix(t(subset(df, select=-c(x,y))), ncol=nrow(df)))
# 
# lat_lon <- data.frame(longitude=df.list[[1]]$x, latitude=df.list[[1]]$y)
# summary(lat_lon)

#------------------------------------------------
start="1949-01-01"; by ="month"

data.NCEP.vwnd <- data_ncep(ncfile1, var="vwnd", level, start, by, mask=land.mask)
data.NCEP.uwnd <- data_ncep(ncfile2, var="uwnd", level, start, by, mask=land.mask)
data.NCEP.air  <- data_ncep(ncfile3, var="air", level, start, by, mask=land.mask)
data.NCEP.hgt  <- data_ncep(ncfile4, var="hgt", level, start, by, mask=land.mask)
data.NCEP.shum <- data_ncep(ncfile5, var="shum", level, start, by, mask=land.mask)
data.NCEP.slp  <- data_ncep(ncfile6, var="slp", level=1, start, by, mask=land.mask)

lat_lon <- data.NCEP.slp$xy
summary(lat_lon)

#------------------------------------------------
data.list <- c("data.NCEP.vwnd","data.NCEP.uwnd","data.NCEP.air","data.NCEP.hgt","data.NCEP.shum","data.NCEP.slp")
for(data in data.list){
# identify NaN grid
data.NCEP <- eval(parse(text=data))
Ind_NCEP <- which(apply(data.NCEP$data[[1]], 2, function(m) sum(is.na(m)))==0)
Ind_NCEP_NaN <- which(apply(data.NCEP$data[[1]], 2, function(m) sum(is.na(m)))==nrow(data.NCEP$data[[1]]))
length(Ind_NCEP)+length(Ind_NCEP_NaN)

print(sum(!is.na(data.NCEP$data[[1]][,Ind_NCEP[1]])))
print(sum(!is.na(data.NCEP$data[[1]][,Ind_NCEP[length(Ind_NCEP)]])))

}

#------------------------------------------------
if(flag.save){
  save(data.NCEP.vwnd, file="data.NCEP.vwnd.Rdata")
  save(data.NCEP.uwnd, file="data.NCEP.uwnd.Rdata")
  save(data.NCEP.air,  file="data.NCEP.air.Rdata")
  save(data.NCEP.hgt,  file="data.NCEP.hgt.Rdata")
  save(data.NCEP.shum, file="data.NCEP.shum.Rdata")
  save(data.NCEP.slp,  file="data.NCEP.slp.Rdata")
  
  save(Ind_NCEP, file="Ind_NCEP.Rdata")
  save(lat_lon, file="lat-lon_NCEP.Rdata")
  write.csv(lat_lon, file="Global_NCEP_Grid.csv")
}




