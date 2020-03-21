# Nesting model for GCM bias correction - for calibration and validation periods
# R 3.6.0
# Ze Jiang - ze.jiang@hotmail.com
# 10/01/2020
#-------------------------------------------------------------------------------------------------
find.val=function(x,val)
{ # function to find index(es) of data in x matching 'val'
  temp=seq(along=x)[x==val]
  temp=temp[is.na(temp)!=T]
  temp
}
#-------------------------------------------------------------------------------------------------
#' Function: Recursive Nested bias correction for calibration - atmospheric variables
#'
#' @param n           The number of columns (e.g., stations)
#' @param y           A matrix of time series atmospheric variables
#' @param obs.sta     Observed (Target) data statistics
#' @param start.yr    Start year
#' @param end.yr      End year
#' @param fun         A function to compute the summary statistics, "mean" or "sum".
#' @param tol         Tolerance for small values, default 0.1
#' @param iteration   The number of iteration to run the bias correction.
#'
#' @return
#' @export
#'
#' @examples
#'
rnbc.mod.cal <- function(n,y,obs.sta,start.yr,end.yr,fun="mean",tol=0.1, iteration=3){

  rnbc.cal <- function(y,x){
    nest.mod.cal(n=ncol(y),y,obs.sta, start.yr, end.yr, fun, tol)$gcm.cor
  }

  Reduce('rnbc.cal', 1:iteration, init = y, accumulate = FALSE)

}
#-------------------------------------------------------------------------------------------------
#' Function: Recursive Nested bias correction for validation - atmospheric variables
#'
#' @param n           The number of columns (e.g., stations)
#' @param y           A matrix of time series atmospheric variables
#' @param obs.sta     Observed (Target) data statistics
#' @param mod.sta     Modeled data statistics at calibration
#' @param start.yr    Start year
#' @param end.yr      End year
#' @param fun         A function to compute the summary statistics, "mean" or "sum".
#' @param tol         Tolerance for small values, default 0.1
#' @param iteration   The number of iteration to run the bias correction.
#'
#' @return
#' @export
#'
#' @examples
#'
rnbc.mod.val <- function(n,y,obs.sta,mod.sta,start.yr,end.yr,fun="mean",tol=0.1, iteration=3){

  rnbc.val <- function(y,x){
    nest.mod.val(n=ncol(y),y,obs.sta, mod.sta, start.yr, end.yr, fun, tol)$gcm.cor
  }

  Reduce('rnbc.val', 1:iteration, init = y, accumulate = FALSE)

}
#-------------------------------------------------------------------------------------------------
#' Function: Calculate atmospheric variables statistics
#'
#' @param n     The number of columns (e.g., stations)
#' @param y     A matrix of time series atmospheric variables
#' @param fun   A function to compute the summary statistics, "mean" or "sum".
#'
#' @return
#' @export
#'
#' @examples
#'
obs.stats <- function(n,y,fun)
{
    y.rho <- array(0,dim=c(n,12))                               # Array to save results for monthly lag one model autocorrelations
    for (i in 1:12)
    {
    if (i==1) ind.i=find.val(cycle(y),i)[-1] else ind.i=find.val(cycle(y),i)     # Find which month of the year we're in
    ind.i1=ind.i-1
    for (j in 1:n)
    {
        if (n==1) y.rho[i] <- cor(as.vector(y[ind.i]),as.vector(y[ind.i1]),use="pair") else    # if only one grid cell, calculate correlation between month i and i-1
        {
            y.rho[j,i] <- cor(as.vector(y[ind.i,j]),as.vector(y[ind.i1,j]),use="pair")         # calculate correlation at multiple grid cells
        }
    }
    }

	# Calculate monthly mean and standard deviation
    y.mean <- array(NA,dim=c(n,12))                                # Calculate model monthly mean (need to save for validation period)
    y.sd <- array(NA,dim=c(n,12))																	 # Calculate model monthly std dev (need to save for validation period)
    for (i in 1:12)
    {
        ind <- find.val(cycle(y),i)                                # find month i
        if (n==1) temp <- y[ind] else temp <- y[ind,]                                            # Find all rainfalls in month i
        if (n ==1) y.mean[i] <- mean(temp) else y.mean[,i] <- apply(temp,2,mean,na.rm=T)                   # calculate monthly mean
        if (n ==1) y.sd[i] <- sd(temp) else y.sd[,i] <- sqrt(apply(temp,2,var,use="pair"))             # Calculate monthly std deviation
    }

	# Aggregate to annual
    #z.hat <- aggregate(y,freq=1)                                                    # Create annual rainfall time series
    z.hat <- aggregate(y,nfrequency=1,FUN=fun)

	# Calculate model lag 1 autocorrelations
    z.rho <- rep(NA,n)
    if (n ==1) z.rho <- acf(z.hat,lag=1,plot=FALSE,na.action=na.pass)$acf[2,,]
    if (n > 1)
    {
        for (i in 1:n) z.rho[i] <- acf(z.hat[,i],lag=1,plot=F,na.action=na.pass)$acf[2,,]  # Calculate lag one autocorrelation for each grid cell
    }
	# Standardise annual totals
    if (n==1) z.mean <- mean(z.hat) else z.mean <- apply(z.hat,2,mean,na.rm=T)          # Find annual mean
    if (n==1) z.sd <- sd(z.hat) else z.sd <- sqrt(apply(z.hat,2,var,use="pair"))        # Find annual std dev

	# List of outputs
  # out <- list(mon.mean = y.mean,mon.sd = y.sd,mon.rho = y.rho,yr.mean = z.mean,yr.sd = z.sd,yr.rho = z.rho)
    out <- list(mon.mean = y.mean,mon.sd = y.sd,mon.rho = y.rho,yr.mean = z.mean,yr.sd = z.sd,yr.rho = z.rho,
                m.mu=apply(as.matrix(y),2,mean),
                m.sd=apply(as.matrix(y),2,sd),
                m.rho= apply(as.matrix(y),2, function(x) acf(x,lag=1,plot=FALSE,na.action=na.pass)$acf[2,,]))
    out
}
#-------------------------------------------------------------------------------------------------
#' Function: Nested bias correction for calibration - atmospheric variables
#'
#' @param n           The number of columns (e.g., stations)
#' @param y           A matrix of time series atmospheric variables
#' @param obs.sta     Observed (Target) data statistics
#' @param start.yr    Start year
#' @param end.yr      End year
#' @param fun         A function to compute the summary statistics, "mean" or "sum".
#' @param tol         Tolerance for small values, default 0.1
#'
#' @return
#' @export
#'
#' @examples
#'
nest.mod.cal <- function(n,y,obs.sta,start.yr,end.yr,fun="mean",tol=0.1)
{
  obs.mon.mean=obs.sta$mon.mean;obs.mon.sd=obs.sta$mon.sd;obs.mon.cor=obs.sta$mon.rho
  obs.yr.mean=obs.sta$yr.mean;obs.yr.sd=obs.sta$yr.sd;obs.yr.cor=obs.sta$yr.rho

	# Define zero tolerance
  # tol = 0.1                                                 # Used to define months and years with ~ 0 rainfall
	nyrs <- end.yr +1 - start.yr									              # Number of years

	# Check for modelled months all zero for any month -
	# add some slight noise to remove any numerical problems
	for (i in 1:12)
	{
    ind.i <- find.val(cycle(y),i)     # Find which month of the year we're in
		if (n > 1)
		{
        for (j in 1:n)
    		{ # If all months are zero and not missing then replace with noise -
          # uniform distn from 0 to tolerance value
    			if (sum(y[ind.i,j]==0&!is.na(y[ind.i,j])) == nyrs) y[ind.i,j] <- runif(nyrs,0,tol)
    		}
		}
	}

	#-----------------------------------------------------------------------------------------------------------monthly
	# Calculate monthly model autocorrelations
    m.rho.mod <- array(0,dim=c(n,12))                               # Array to save results for monthly lag one model autocorrelations
    for (i in 1:12)
    {
        if (i==1) ind.i=find.val(cycle(y),i)[-1] else ind.i=find.val(cycle(y),i)     # Find which month of the year we're in
        ind.i1=ind.i-1
        for (j in 1:n)
        {
            if (n==1) m.rho.mod[i] <- cor(as.vector(y[ind.i]),as.vector(y[ind.i1]),use="pair") else    # if only one grid cell, calculate correlation between month i and i-1
            {
				      cor.temp <- cor(as.vector(y[ind.i,j]),as.vector(y[ind.i1,j]),use="pair")			   # calculate correlation at each grid cells
				      if (is.na(cor.temp)==TRUE) m.rho.mod[j,i] <- 0 else m.rho.mod[j,i] <- cor.temp
      #				m.rho.mod[j,i] <- cor(as.vector(y[ind.i,j]),as.vector(y[ind.i1,j]),use="pair")
            }
        }
    }

	# Standardise data by modelled means and standard deviation
    y.mean <- array(NA,dim=c(n,12))                                # Calculate model monthly mean (need to save for validation period)
    y.sd <- array(NA,dim=c(n,12))																	 # Calculate model monthly std dev (need to save for validation period)
    y.std <- y                                                     # Create matrix to save standardised results in
    for (i in 1:12)
    {
        ind <- find.val(cycle(y),i)                                # find month i
        if (n==1) temp <- y[ind] else temp <- y[ind,]                                            # Find all rainfalls in month i
        if (n==1) y.mean[i] <- mean(temp) else y.mean[,i] <- apply(temp,2,mean,na.rm=T)          # calculate monthly mean
        if (n==1) y.sd[i] <- sd(temp) else y.sd[,i] <- sqrt(apply(temp,2,var,use="pair"))        # Calculate monthly std deviation

        if (n==1) temp2 <- (temp-y.mean[i])/y.sd[i]                               # Minus mean and divide by std dev (scaling function)
        if (n > 1) temp2 <- apply(temp,2,scale)                               # Minus mean and divide by std dev (scaling function)
        if (n==1) y.std[ind] <- temp2 else y.std[ind,] <- temp2                                  # save results into new matrix
    }
    gc() # clean up memory

  #-----------------------------------------------------------------------------------------------------------monthly	correct
	# Remove the autocorrelations already in the model
    y.std.uncor <- y.std                                           # Create matrix to save standardised results in
    if (n ==1) leny <- length(y)
    if (n > 1) leny <- dim(y)[1]
    for (i in 2:leny)
    {
        j <- cycle(y.std)[i]                                         # find month i
        if (n==1) y.std.uncor[i] <- (y.std[i]-m.rho.mod[j]*y.std[i-1])/sqrt(1-m.rho.mod[j]^2)    # take out model autocorrelation
        if (n >1) y.std.uncor[i,] <- (y.std[i,]-m.rho.mod[,j]*y.std[i-1,])/sqrt(1-m.rho.mod[,j]^2)    # take out model autocorrelation
    }

	# Apply Lag one autocorrelation of monthly aggregated data
    y.hat <- y.std.uncor                                            # Create matrix to save standardised results in
    for (i in 2:leny)
    {
        j <- cycle(y.std.uncor)[i]                                  # find month i
        if (n==1) y.hat[i] <- obs.mon.cor[j]*y.hat[i-1]+sqrt(1-obs.mon.cor[j]^2)*y.std.uncor[i] else           # Add in observed monthly autocorrelations (for one grid cell)
        {
            y.hat[i,] <- obs.mon.cor[,j]*y.hat[i-1,]+sqrt(1-obs.mon.cor[,j]^2)*y.std.uncor[i,]                   # Add in observed monthly autocorrelations (if more than one grid cell)
        }
    }

	# Unpack by scaling with observed monthly means and standard deviations
    j <- cycle(y.hat)
    for (i in 1:12)
    {
        ind <- find.val(j,i)                                  # find month i
        if (n==1) y.hat[ind] <- t(t(y.hat[ind])*obs.mon.sd[i]+obs.mon.mean[i]) else    # multiply by obs standard deviation and add obs mean (one grid cell)
        {
            y.hat[ind,] <- t(t(y.hat[ind,])*obs.mon.sd[,i]+obs.mon.mean[,i])             # multiply by obs standard deviation and add obs mean (more than one grid cell)
        }
    }

    # If any values are less than zero, set to zero like rainfall
#    y.hat[(is.na(y.hat)!=T)&(y.hat<0)] <- 0
#    y.hat[(is.na(y.hat)!=T)&(abs(y.hat)<tol)] <- tol
    y.hat <- ts(y.hat,start=c(start.yr,1),end=c(end.yr,12),freq=12)											 # make matrix into time series so we can add to annual data
    gc() # clean up memory

	#-----------------------------------------------------------------------------------------------------------annual
	# Sum to form annual totals from corrected monthly data
    #z.hat <- aggregate(y.hat,nfrequency=1)                                                     # Create annual rainfall time series
    z.hat <- aggregate(y.hat,nfrequency=1,FUN=fun)

	# Calculate model lag 1 autocorrelations
    y.rho.mod <- rep(NA,n)
    if (n==1) y.rho.mod <- acf(z.hat,lag=1,plot=F,na.action=na.pass)$acf[2,,]  # Calculate lag one autocorrelation for each grid cell
    if (n >1)
    {
        for (i in 1:n) y.rho.mod[i] <- acf(z.hat[,i],lag=1,plot=F,na.action=na.pass)$acf[2,,]  # Calculate lag one autocorrelation for each grid cell
    }
	# Standardise annual totals
    z.std <- scale(z.hat)																																# Standarise by removing mean and standard deviation
    if (n==1) z.mean <- mean(z.hat) else z.mean <- apply(z.hat,2,mean,na.rm=T)          # Find annual mean
    if (n==1) z.sd <- sd(z.hat) else z.sd <- sqrt(apply(z.hat,2,var,use="pair"))        # Find annual std dev

    #-----------------------------------------------------------------------------------------------------------annual correct
    # Remove any lag 1 autocorrelations from annual data
    z.std.uncor <- z.std
    for (i in 2:dim(z.std)[1])
    {
        if (n==1)	z.std.uncor[i] <- (z.std[i]-y.rho.mod*z.std[i-1])/sqrt(1-y.rho.mod^2)          # Remove model lag one auto-correlation
        if (n >1)
        {
            z.std.uncor[i,] <- (z.std[i,]-y.rho.mod*z.std[i-1,])/sqrt(1-y.rho.mod^2)          # Remove model lag one auto-correlation
        }
    }

	# Correct standardised yearly values
    z.tick <- z.std.uncor
    for (i in 2:dim(z.tick)[1])
    {
        if (n==1) z.tick[i] <- obs.yr.cor*z.tick[i-1]+sqrt(1-obs.yr.cor^2)*z.std.uncor[i] else       # Add in observed annual lag one auto-cor (one grid cell)
        {
           z.tick[i,]=obs.yr.cor*z.tick[i-1,]+sqrt(1-obs.yr.cor^2)*z.std.uncor[i,]                   # Add in observed annual lag one auto-cor (one grid cell)
        }
    }

	# Unpack by scaling with observed annual mean and standard deviation
    z.tick <- t(t(z.tick)*obs.yr.sd+obs.yr.mean)                                    # multiply by obs standard deviation and add obs mean
    z.tick <- ts(z.tick,start=start.yr,end=end.yr)                                  # Make into time series

    gc() # clean up memory

   	#-----------------------------------------------------------------------------------------------------------correction factor
	# Create y tick
    y[between(as.vector(y),0,tol)]=tol; y[between(as.vector(y),-tol,0)]= -tol
    month.const <- as.matrix(y.hat)/as.matrix(y)                                    # Monthly correction factor - new monthly total/old monthly total
#   month.const[y < tol] <- 0                                                       # if monthly total (old) is too small - make factor zero (don't want to divide by zero)
    month.const[is.na(month.const)] <- 0                                            # If monthly total is missing, make 0

    z.hat[between(as.vector(z.hat),0,tol)]=tol; z.hat[between(as.vector(z.hat),-tol,0)]= -tol
    year.const <- z.tick/z.hat                                                      # Yearly correction factor - new Yearly total/old Yearly total from corrected monthly data?!
#   year.const[z.hat < tol] <- 0                                                    # if yearly total (old) is too small - make factor zero (don't want to divide by zero)
    year.const[is.na(year.const)] <- 0                                              # If yearly total is missing, make 0

    const <- month.const                                                            # We now need to multiple monthly constant by yearly constant
    for (j in 1:12)
    {
        ind <- find.val(cycle(y.hat),j)                                             # Find all month i's
        if (n==1) const[ind]<- month.const[ind]*year.const else const[ind,] <- month.const[ind,]*year.const  # Multiply each month i by corresponding year constant
    }

    if (n==1) y.tick <- as.vector(const)*y                                                               # Correct monthly values
    if (n >1) y.tick <- const*y
    y.tick <- ts(y.tick,start=start(y),end=end(y),freq=12)                          # Make timeseries
    gc()
    #     out.list <- list(gcm.cor=y.tick,gcm.raw=y,
    #                      mod.mon.mean=y.mean,mod.mon.sd=y.sd,mod.mon.cor=m.rho.mod,
    # 											  mod.yr.mean=z.mean,mod.yr.sd=z.sd,mod.yr.cor=y.rho.mod) #statictis from corrected monthly data

    out.list <- list(gcm.cor=y.tick,gcm.raw=y,
                     mod.sta=list(mon.mean=y.mean,mon.sd=y.sd,mon.rho=m.rho.mod, #model statistics for validation
                                  yr.mean=z.mean,yr.sd=z.sd,yr.rho=y.rho.mod)
                    )
    out.list                                                                        # Write out results list
}
#-------------------------------------------------------------------------------------------------
#' Function: Nested bias correction for validation - atmospheric variables
#'
#' @param n           The number of columns (e.g., stations)
#' @param y           A matrix of time series atmospheric variables
#' @param obs.sta     Observed (Target) data statistics
#' @param mod.sta     Modeled data statistics at calibration
#' @param start.yr    Start year
#' @param end.yr      End year
#' @param fun         A function to compute the summary statistics, "mean" or "sum".
#' @param tol         Tolerance for small values, default 0.1
#'
#' @return
#' @export
#'
#' @examples
#'
nest.mod.val <- function(n,y,obs.sta,mod.sta,start.yr,end.yr,fun="mean",tol=0.1)
{
  obs.mon.mean=obs.sta$mon.mean;obs.mon.sd=obs.sta$mon.sd;obs.mon.cor=obs.sta$mon.rho
  obs.yr.mean=obs.sta$yr.mean;obs.yr.sd=obs.sta$yr.sd;obs.yr.cor=obs.sta$yr.rho

  mod.mon.mean=mod.sta$mon.mean;mod.mon.sd=mod.sta$mon.sd;mod.mon.cor=mod.sta$mon.rho
  mod.yr.mean=mod.sta$yr.mean;mod.yr.sd=mod.sta$yr.sd;mod.yr.cor=mod.sta$yr.rho

	# Define zero tolerance
  # tol <- 0.1                                                   # Used to define months and years with ~ 0 rainfall
  nyrs <- end.yr +1 - start.yr									                 # Number of years

	# Check for modelled months all zero for any month - add some slight noise to remove any numerical problems
	for (i in 1:12)
	{
    ind.i <- find.val(cycle(y),i)     # Find which month of the year we're in
		if (n > 1)
		{
        for (j in 1:n)
    		{
    			if (sum(y[ind.i,j]==0&!is.na(y[ind.i,j])) == nyrs) y[ind.i,j] <- runif(nyrs,0,tol)		# If all months are zero and not missing then replace with noise - uniform distn from 0 to tolerance value
    		}
		}
	}

	# Standardise data by modelled means and standard deviation for calibration period
    y.std <- y
    for (i in 1:12)
    {
        # if (n==1) if(mod.mon.sd[i] < 1) mod.mon.sd[i] <- 1 # for small value like rainfall
        # if (n > 1) mod.mon.sd[,i][mod.mon.sd[,i] < 1] <- 1 # for small value like rainfall
        ind <- find.val(cycle(y),i)
        if (n==1) temp <- y[ind] else temp <- y[ind,]
        if (n==1) temp2 <- scale(temp,center=as.vector(mod.mon.mean[i]),scale=mod.mon.sd[i])
        if (n >1) temp2 <- scale(temp,center=as.vector(mod.mon.mean[,i]),scale=mod.mon.sd[,i])
        if (n ==1) y.std[ind] <- temp2 else y.std[ind,] <- temp2
    }
    gc()

	# Remove the autocorrelations from calibration period already in the model
    y.std.uncor <- y.std
    if (n ==1) leny <- length(y)
    if (n > 1) leny <- dim(y)[1]
    for (i in 2:leny)
    {
        j <- cycle(y.std)[i]
        if (n==1) y.std.uncor[i] <- (y.std[i]-mod.mon.cor[j]*y.std[i-1])/sqrt(1-mod.mon.cor[j]^2)
        if (n >1) y.std.uncor[i,] <- (y.std[i,]-mod.mon.cor[,j]*y.std[i-1,])/sqrt(1-mod.mon.cor[,j]^2)
    }

    # Create y hat
    # Calculate correlation for future observed
    # Apply Lag one autocorrelation of monthly aggregated data
    y.hat <- y.std.uncor
    for (i in 2:leny)
    {
        j <- cycle(y.std.uncor)[i]
        if (n==1) y.hat[i] <- obs.mon.cor[j]*y.hat[i-1]+sqrt(1-obs.mon.cor[j]^2)*y.std.uncor[i] else
        {
            y.hat[i,] <- obs.mon.cor[,j]*y.hat[i-1,]+sqrt(1-obs.mon.cor[,j]^2)*y.std.uncor[i,]
        }
    }

    # Unpack by scaling with observed monthly means and standard deviations
    j <- cycle(y.hat)
    for (i in 1:12)
    {
        ind <- find.val(j,i)
        if (n==1) y.hat[ind] <- t(t(y.hat[ind])*obs.mon.sd[i]+obs.mon.mean[i]) else
        {
            y.hat[ind,] <- t(t(y.hat[ind,])*obs.mon.sd[,i]+obs.mon.mean[,i])
        }
    }
#    y.hat[(is.na(y.hat)!=T)&(y.hat<0)]=0
#    y.hat[(is.na(y.hat)!=T)&(abs(y.hat)<tol)]=tol
    y.hat <- ts(y.hat,start=c(start.yr,1),end=c(end.yr,12),freq=12)
    gc()

    # Create z hat
    # Sum to form annual totals
    #z.hat <- aggregate(y.hat,nfrequency=1)
    z.hat <- aggregate(y.hat,nfrequency=1,FUN=fun)

    # Standardise annual totals
    z.std <- scale(z.hat,center=mod.yr.mean,scale=mod.yr.sd)

    # Remove any lag 1 autocorrelations from annual data
    z.std.uncor <- z.std
    for (i in 2:dim(z.std)[1])
    {
        if (n==1) z.std.uncor[i] <- (z.std[i]-mod.yr.cor*z.std[i-1])/sqrt(1-mod.yr.cor^2)
        if (n >1) z.std.uncor[i,] <- (z.std[i,]-mod.yr.cor*z.std[i-1,])/sqrt(1-mod.yr.cor^2)
    }

    # Correct standardised yearly values
    z.tick <- z.std.uncor
    for (i in 2:dim(z.tick)[1])
    {
        if (n==1) z.tick[i] <- obs.yr.cor*z.tick[i-1]+sqrt(1-obs.yr.cor^2)*z.std.uncor[i] else
        {
            z.tick[i,] <- obs.yr.cor*z.tick[i-1,]+sqrt(1-obs.yr.cor^2)*z.std.uncor[i,]
        }
    }

    # Unpack by scaling with observed annual mean and standard deviation
    z.tick <- t(t(z.tick)*obs.yr.sd+obs.yr.mean)                                    # multiply by obs standard deviation and add obs mean
    z.tick <- ts(z.tick,start=start.yr,end=end.yr)                                  # Make into time series

    gc() # clean up memory

    # Create y tick
    y[between(as.vector(y),0,tol)]=tol; y[between(as.vector(y),-tol,0)]= -tol
    month.const <- as.matrix(y.hat)/as.matrix(y)
#   month.const[y < tol] <- 0
    month.const[is.na(month.const)] <- 0

    z.hat[between(as.vector(z.hat),0,tol)]=tol; z.hat[between(as.vector(z.hat),-tol,0)]= -tol
    year.const <- z.tick/z.hat
#   year.const[z.hat < tol] <- 0
    year.const[is.na(year.const)] <- 0

    const <- month.const
    for (j in 1:12)
    {
        ind <- find.val(cycle(y.hat),j)
        if (n==1) const[ind]<- month.const[ind]*year.const else const[ind,] <- month.const[ind,]*year.const
    }

    if (n==1) y.tick <- as.vector(const)*y                                         # Correct monthly values
    if (n >1) y.tick <- const*y
    y.tick <- ts(y.tick,start=start(y),end=end(y),freq=12)                         # Create timeseries
    gc()
    out.list <- list(gcm.cor=y.tick,gcm.raw=y)
    out.list                                                                       # Write out results
}

