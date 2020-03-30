#-----------------------------------------------
#' Function: Load NCEP Reanalysis 1 data
#'
#' @param ncfile  A nc file
#' @param var     Variable to extract
#' @param level   Pressure level
#' @param start   Start Date
#' @param by      Time step
#' @param mask    A shapefile mask, land mask, ocean mask or country mask.
#'
#' @return
#' @export
#'
#' @examples
#'
data_ncep <- function(ncfile, var="slp", level, start="1949-01-01", by="month", mask)
{
  nc2ras <- sapply(level, function(i) mask(rotate(brick(ncfile, varname = var, level=i)), mask))
  #nc2ras

  df.list <- lapply(nc2ras, function(ls) raster::as.data.frame(ls, xy = TRUE, na.rm=T))
  data.NCEP <- lapply(df.list, function(df) matrix(t(subset(df, select=-c(x,y))), ncol=nrow(df)))

  lat_lon <- data.frame(longitude=df.list[[1]]$x, latitude=df.list[[1]]$y)
  #summary(lat_lon)

  time <- seq(as.Date(start), by=by, length.out = nrow(data.NCEP[[1]]))

  return(list(data=data.NCEP,
              xy=lat_lon,
              Date=time,
              levels=level))
}
#-----------------------------------------------
#' Function: Load UDEL Air Temperature & Precipitation
#'
#' @param ncfile  A nc file
#' @param var     Variable to extract
#' @param level   Pressure level
#' @param start   Start Date
#' @param by      Time step
#' @param mask    A shapefile mask, land mask, ocean mask or country mask.
#'
#' @return
#' @export
#'
#' @examples
#'
data_udel <- function(ncfile, var="precip", level, start="1900-01-01", by="month", mask)
{
  nc2ras <- mask(rotate(brick(ncfile,varname = var)),mask)
  #nc2ras

  df <- raster::as.data.frame(nc2ras, xy = TRUE, na.rm=T)
  data.UDEL <- matrix(t(subset(df, select=-c(x,y))), ncol=nrow(df))

  lat_lon <- data.frame(longitude=df$x, latitude=df$y)
  #summary(lat_lon)

  time <- seq(as.Date(start), by=by, length.out = nrow(data.UDEL))

  return(list(data=data.UDEL,
              xy=lat_lon,
              Date=time,
              levels=level))

}
#-----------------------------------------------
#' Function: Load NCEP Reanalysis 1 data
#'
#' @param ncfile  A nc file
#' @param var     Variable to extract
#' @param level   Pressure level
#' @param start   Start Date
#' @param by      Time step
#' @param grid    Resample grids, e.g., NCEP reanalysis grids.
#' @param mask    A shapefile mask, e.g., land mask, ocean mask or country mask.
#'
#' @return
#' @export
#'
#' @examples
#'
data_cmip <- function(ncfile, var="psl", level, start="1850-01-01", by="month", grid, mask)
{
  nc2ras <- sapply(level, function(i) rotate(brick(ncfile, varname = var, level=i)))
  #nc2ras

  df.ras.GCM <- raster(SpatialPixelsDataFrame(grid, grid@data))
  df.ras.GCM <- lapply(nc2ras, function(ras) mask(resample(ras, df.ras.GCM, method="bilinear"),mask))
  #df.ras.GCM

  df.list <- lapply(df.ras.GCM, function(ls) raster::as.data.frame(ls, xy = TRUE, na.rm=T))
  data.GCM <- lapply(df.list, function(df) matrix(t(subset(df, select=-c(x,y))), ncol=nrow(df)))

  lat_lon <- data.frame(longitude=df.list[[1]]$x, latitude=df.list[[1]]$y)
  #summary(lat_lon)

  time <- seq(as.Date(start), by=by, length.out = nrow(data.GCM[[1]]))

  return(list(data=data.GCM,
              xy=lat_lon,
              Date=time,
              levels=level))
}
