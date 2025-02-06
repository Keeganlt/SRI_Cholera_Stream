#Note....the number sign meens the rest of the line is a comment and
#will not be run. Comment your code a lot!!! Your life will be much easier.

#First thing we need to do is to load the library for solving ODEs
library(deSolve)

#Here are some constants representing default values. You may
#want to change these default values instead of changing how you
#call the model.
DEF.BETA <- 1
DEF.GAMMA <- .5
DEF.S <- .999
DEF.I  <- .001
DEF.R <- 0
MAX.TIME <- 365


#This function actually runs your SIR model.
# Inputs -
#    beta - the force of infection
#    gamma - the recovery rate
#    initial.state - a represents the initial state.
#Returns
#    an ode output object. to plot this object use the
#    functions plot.ode and legend.ode
#
#
runSIR <- function(beta=DEF.BETA,
                   gamma=DEF.GAMMA,
                   initial.state= c(S=DEF.S, I=DEF.I, R=DEF.R),
                   max.time = MAX.TIME,
                   freq.dependent=FALSE) {



  #If the model is frequency dependent we modify beta
  #based on the total populations size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["S"]+initial.state["I"]+initial.state["R"],
                         1)


  #create the parameter vector.
  param <- c(beta=beta/beta.divisor, gamma=gamma)

 #Sequence of times at which we want estimates..
 #   ..here we say daily until max.time
  times <- seq(0,max.time,1)

  #Run the ODE-Solver
  sir.output <- lsoda(initial.state, times, dx.dt.SIR, param)

  #now we will return the output.
  return(sir.output)
}

#This function does the model time step. We call it dx.dt.SIR for an SIR model.
#Inputs:
#     t - current time
#     y - current state of the system
#     param - the set of parameters. in this case beta and gamma.
#Returns:
#     a list representing the changes in S, I, and R on this time step.
dx.dt.SIR <- function(t, y, param) {
  #This is the equation to calculate the change in susceptibles
  dS <- -param["beta"] * y["S"] * y["I"]

  #Here we calculate the change in infectious individuals
  dI <- param["beta"] * y["S"] * y["I"] - param["gamma"] * y["I"]

  #finally we calculate the changes in R
  dR <- param["gamma"] * y["I"]

  #Return the results!!! If you add something you need to make sure that
  #the deltas are returned in the exact same order as the corresponding
  #compartments in the state vector passed in.
  return(list(c(dS, dI, dR)))
}

#This function actually runs your SIR model, tracking cumulative incidence.
# Inputs -
#    beta - the force of infection
#    gamma - the recovery rate
#    initial.state - a represents the initial state.
#Returns
#    an ode output object. to plot this object use the
#    functions plot.ode and legend.ode
#
#
runIncidenceTrackingSIR <- function(beta=DEF.BETA,
                                    gamma=DEF.GAMMA,
                                    initial.state=
                                    c(S=DEF.S, CI = DEF.I,  I=DEF.I, R=DEF.R),
                                    max.time = MAX.TIME,
                                    freq.dependent=FALSE) {

  #If the model is frequency dependent we modify beta
  #based on the total populations size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["S"]+initial.state["I"]+initial.state["R"],
                         1)


  #print(beta/beta.divisor)
  #create the parameter vector.
  param <- c(beta=beta/beta.divisor, gamma=gamma)

 #Sequence of times at which we want estimates..
 #   ..here we say daily until max.time
  times <- seq(0,max.time,1)


  #Run the ODE-Solver
  sir.output <- lsoda(initial.state, times, dx.dt.IncidenceTrackingSIR, param)

  #now we will return the output.
  return(sir.output)
}

################################################
# Here are some functions that will help you plot
# SIR models. By modifying these functions or writing
# your own you can make your own custom plots.
################################################

#Creates a plot of the output of an ODE
#Input:
#    lsoda.output = the output from an ODE solver.
#    tdivisor = the divisor for the time axis,
#          useful for changing the ploted time scale
plot.ode <- function(lsoda.output, t.divisor=1)  {
  #check if this looks like lsoda.output,
  #if not, through an error.
  if(!is.matrix(lsoda.output))
    stop("'lsoda.output' must be output from odesolve")

  #make an empty plot with the appropriate time axis, etc.
  plot(0, 0, xlim=range(lsoda.output[,"time"])/t.divisor,
       ylim=range(lsoda.output[,2:dim(lsoda.output)[2]]),
       xlab="Time",
       ylab="State",
       type="n")

  # plot the lsoda.output using the matlines function
  # whcih takes a single x vector and matrix of multiple y
  # values.
  matlines(lsoda.output[,"time"]/t.divisor,
           lsoda.output[,2:dim(lsoda.output)[2]], lwd=2)
}


#Plots the legend for an ODE.
legend.ode <- function(x , y , lsoda.output)
    {
    #how many states are in the model? (the number of columns excluding time)
    num.state <- dim(lsoda.output)[2] - 1

    #draw a legend where names are derived from the columns
    #of lsoda.output
    legend(x, y, legend=colnames(lsoda.output)[2:dim(lsoda.output)[2]],
           col=rep(1:6, l=num.state), lty=rep(1:6, l=num.state), lwd=2)
}



#This function does the model time step. We call it dx.dt.SIR for an SIR model.
#Inputs:
#     t - current time
#     y - current state of the system
#     param - the set of parameters. in this case beta and gamma.
#Returns:
#     a list representing the changes in S, CI, I, and R on this time step.
dx.dt.IncidenceTrackingSIR <- function(t,y,param) {
  
  #this tracks the change in cumulative incidence
  dCI <- param["beta"] * y["S"] * y["I"]
  
  #get the normal dx.dt
  sirRes <-  dx.dt.SIR(t,y,param)
  
  #Return the results!!! If you add something you need to make sure that
  #the deltas are returned in the exact same order as the corresponding
  #compartments in the state vector passed in.
  rc <- list(c(sirRes[[1]][1], dCI, sirRes[[1]][2], sirRes[[1]][3]))
  
  return(rc)
  
}


#Function for using OPTIM to fit beta, gamma, or both
#
#Parameters -
#    act.incidences - a single incidence curve
#                or a list of incidence curves to fit
#    initial.states - a single initial state or a list there of
#    beta - the starting value for beta if beta is being fit,
#             the actual value if it is not.
#    gamma - the starting value of gamma if it is being fit,
#             the actual value if it is not.
#    fit - what to fit, can be  "beta", "gamma" or "both",
#             default is "beta"
#    freq.dependent - use frequency dependent transmission? default is
#             FALSE
#    pretty.graphics - show us some pretty graphics as we go?
#
#Return -
#    optimal values of beta and gamma
fitSIR <- function(act.incidences,
                   initial.states,
                   beta,
                   gamma,
                   fit="beta",
                   freq.dependent = FALSE,
                   pretty.graphics = FALSE) {
  #ensure the listiness of our data
  if (is.data.frame(act.incidences)) {
    act.incidences <- list(act.incidences)
    initial.states <- list(initial.states)
  }
  
  #do the appropriate optim fit.
  if (fit=="beta") {
    
    #Run optim!!!
    opt.rc <- nlminb(start=beta,
                     objective=function(beta) {
                       rc <- 0
                       
                       for (i in 1:length(act.incidences)){
                         
                         rc <- rc+sir.MSE(beta=beta, gamma=gamma,
                                          initial.state=initial.states[[i]],
                                          freq.dependent=freq.dependent,
                                          act.incidence=act.incidences[[i]],
                                          pretty.graphics=pretty.graphics)
                       }
                       return(rc)
                     }, lower=0)
    
    beta <- opt.rc$par
    best.mse <- opt.rc$objective
  } else if (fit=="gamma") {
    #TODO: Implement
  } else {
    #TODO: Implement
  }
  
  return(c(beta=beta, gamma=gamma, mse=best.mse))
}




###################################
# Utility Methods
###################################

#Function that creates an object representing epidemic
#curve data.
#
#Parameters
#  onset.time  - a date object representing the time of onset
#  case.status - an indicator of case status, by default one
#      for everyone listed (non-cases are only needed for final
#                           size calculation)
#
#Returns
#  an epi.curve object, any data.frame with an
#    onset time column and a case status column
#    can be cast to this class
epi.curve <- function(onset.time, case.status = 1) {
  rc <- data.frame(onset.time = onset.time)
  rc$is.case <- case.status
  
  attr(rc,"class") <- c("epi.curve", "data.frame")
  return(rc)
}


#function creates an object represent an incidence curve
#
#Parameters
# time - the times
# incidence - the number of people infected on that day.
#
#Return -
# a incidence curve object
incidence.curve <- function(time, incidence) {
  
  
  rc <- data.frame(time=time, incidence=incidence)
  
  attr(rc,"class") <- c("incidence.curve", "data.frame")
  return(rc)
}

#turns an incidence curve into an epi curve
#
#Parameters -
#  incidence.curve - an incidence curve
#
#Return -
#  an epidemci curve
incidence.curve.to.epi.curve <- function(ic,start.date=as.Date("2007-01-01")) {
  #make there be whole numbers of individuals
  ic <- round(ic)
  
  #Now turn this into a row per person data set.
  onset.times <- c()
  
  for (i in 1:nrow(ic)) {
    onset.times <-
      c(onset.times, rep(ic[i,"time"], ic[i,"incidence"]))
  }
  
  #make it have an actual time
  onset.times <-start.date+onset.times
  
  #turn into an epi.curve object
  ec <- epi.curve(onset.times)  
}

epi.curve.to.incidence.curve <- function(ec, max.day=NULL) {
  #find the start date
  start.date <- min(ec$onset.time)-1
  
  #convert incidence curve to daily numerics
  ec$onset.time <- as.numeric(ec$onset.time-start.date)
  
  #create daily counts
  if (is.null(max.day)) {
    max.day <- max(ec$onset.time)
  }
  
  ic <- matrix(nrow=max.day, ncol=2)
  ic[,1] <- 1:max.day
  for (i in 1:max.day) {
    ic[i,2] <- sum(ec$onset.time==i)
  }
  
  return(incidence.curve(ic[,1], ic[,2]))
}



#Function that plots an epi curve
#
#Parameters
# epi.curve - the curve to plot.
# ... - other paramteters passed to the underlying plot command
#
#Returns
# Nothing
plot.epi.curve <- function(epi.curve, freq=TRUE,
                           xlab="day", ylab="# of cases",
                           main="Epidemic Curve",
                           days.of.interest=NULL,
                           ...) {
  
  #Days of the epidemic
  if (is.null(days.of.interest)){
    days.of.interest <- seq(min(epi.curve$onset.time, na.rm=TRUE),
                            max(epi.curve$onset.time, na.rm=TRUE)+1,
                            "day")
  }
  
  #plot it in a histogram with appropriate axes, etc.
  tmp.hist <- hist(epi.curve$onset.time, breaks = days.of.interest-.5,
                   plot=FALSE)
  
  
  plot(tmp.hist,
       freq=freq,
       xlab=xlab,
       ylab=ylab,
       axes=FALSE,
       main=main,...)
  
  
  axis(1, at = days.of.interest,
       labels= as.character(days.of.interest,"%b%d"), las=2, cex.axis=.75)
  axis(2)
  
}

#Gets the daily cumulative cases
#
#parameters
#  epi.curve - the epidemic curve
#  prop - should this be proportional?
#
#return
#  cum.cases - the cumulative cases
cum.epi.curve <- function(epi.curve, prop = FALSE) {
  #Days of the epidemic
  days.of.interest <- seq(min(epi.curve$onset.time, na.rm=TRUE),
                          max(epi.curve$onset.time, na.rm=TRUE),
                          "day")
  
  #cumulative cases plot
  cumulative.cases <- sapply(days.of.interest, function(x) {
    sum(epi.curve$onset.time<=x, na.rm=TRUE)
  })
  
  if (prop) {
    cumulative.cases <- cumulative.cases/sum(epi.curve$is.case)
  }
  
  return(data.frame(day=days.of.interest, cum.cases = cumulative.cases))
}

#Function plots the cumulative cases of the epidemic curve.
#Useful for deciding when exponential growth might have stopped
#
#parameters
#  epi.curve - the epidemic curve to plot
#  prop - plot the proportion instead of absolute #?
#  ... - additional parameters to plot
cumCasesPlot <- function (epi.curve, prop=FALSE,...) {
  
  cum.curve <- cum.epi.curve(epi.curve, prop)
  
  plot(cum.curve$day, cum.curve$cum.cases,
       axes=FALSE,
       ylab="cumulative # of cases",
       xlab="day",
       ...)
  
  axis(1, at = cum.curve$day,
       labels= as.character(cum.curve$day,"%b%d"), las=2, cex.axis=.75)
  axis(2)
}

#Function plot the exponential growth rate for
#the epidemic over time. This is the slope of
#the log curve, we are usually intersted in the
#smoothed version
#
#parameters -
#   epi.curve - the epidemic curve to analyze
#   ... -  additional parameters to plot
growthRatePlot <- function(epi.curve, spar=.5,...) {
  #make curve with log values and log slope
  cum.curve <-cum.epi.curve(epi.curve)
  cum.curve$logval <- log(cum.curve$cum.cases)
  cum.curve$logslope <- 0
  cum.curve$logslope[1:(nrow(cum.curve)-1)] <-
    cum.curve$logval[2:(nrow(cum.curve))] -
    cum.curve$logval[1:(nrow(cum.curve)-1)]
  
  plot(cum.curve$day, cum.curve$logslope, axes=FALSE,
       ylab="r", ...)
  
  axis(1, at = cum.curve$day,
       labels= as.character(cum.curve$day,"%b%d"), las=2, cex.axis=.75)
  axis(2)
  
  #now lets calculate the smooth spline version
  cum.curve$smooth.logval <-smooth.spline(cum.curve$day,
                                          log(cum.curve$cum.cases),
                                          spar=spar)$y
  cum.curve$smooth.logslope <- 0
  cum.curve$smooth.logslope[1:(nrow(cum.curve)-1)] <-
    cum.curve$smooth.logval[2:(nrow(cum.curve))] -
    cum.curve$smooth.logval[1:(nrow(cum.curve)-1)]
  
  lines(cum.curve$day, cum.curve$smooth.logslope, lty=2, lwd=2)
}



#Compare the MSE of the predicted curve to the
#actual epi.curve passed in.
#Parametes -
#    act.incidence - the actual incidence curve
#    pred.incidence - the predicted incidence curve
#
#Return -
#    sum((act.incidence-pred.inicidence)^2)
incidence.MSE <- function(act.incidence,
                          pred.incidence) {
  #make sure they are equal
  if (nrow(act.incidence)>nrow(pred.incidence)) {
    time <- act.incidence$time
    incidence <- c(pred.incidence$incidence,
                   rep(0,nrow(act.incidence)-nrow(pred.incidence)))
    pred.incidence <- incidence.curve(time,incidence)
  } else if (nrow(act.incidence)<nrow(pred.incidence)) {
    time <- pred.incidence$time
    incidence <- c(act.incidence$incidence,
                   rep(0,nrow(pred.incidence)-nrow(act.incidence)))
    act.incidence <- incidence.curve(time,incidence)
  }
  
  #ASSERT: nrow(act.incidence)==nrow(pred.incidence)
  return(sum((act.incidence$incidence-pred.incidence$incidence)^2))
}


#Takes SIR parameters, generates incidence curves,
#and returns the MSE of the incidence curve compared
#to the data
sir.MSE <-  function(beta, gamma,
                     initial.state,
                     freq.dependent,
                     act.incidence,
                     pretty.graphics = FALSE) {
  
  #Generate a deterministic epidemic based on the current
  #parameters
  itSIR.res <- runIncidenceTrackingSIR(beta=beta,
                                       gamma=gamma,
                                       initial.state=initial.state,
                                       max.time=4*max(act.incidence$time),
                                       freq.dependent=freq.dependent)
  
  #make an incidnece curve from the epidemic
  pred.incidence <- itSIR.to.incidence.curve(itSIR.res)
  
  rc <- incidence.MSE(act.incidence, pred.incidence)
  
  if (pretty.graphics) {
    plot(act.incidence, type="l", col="blue", lwd=2)
    lines(pred.incidence, col="red", lwd=2, lty=2)
  }
  
  return(rc)
}

#transforms the output of an IncidenceTrackingSIR to
#timestep incidence data.
itSIR.to.incidence.curve <- function(itSIR.output) {
  
  rc <- matrix(nrow=nrow(itSIR.output), ncol=2)
  colnames(rc) <- c("time", "incidence")
  rc[,"time"] <- itSIR.output[,"time"]
  
  last.CI <- c(0,itSIR.output[1:(nrow(itSIR.output)-1),"CI"])
  rc[,"incidence"] <- itSIR.output[,"CI"]-last.CI
  
  return(incidence.curve(rc[,"time"], rc[,"incidence"]))
}
