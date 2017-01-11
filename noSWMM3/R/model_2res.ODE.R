######################################################################################
# 2 cascade reservoir model
#
# J. Pastorek, DEC 2016
# based on a model for North England (O.W. - Feb 2016) (Earlier version - Del Giudice)
######################################################################################

# References

# Breinholt, A., F. O. Thordarson, J. K. Mller, M. Grum, P. S. Mikkelsen, and H. Madsen
# 684 (2011), Grey-box modelling of flow in sewer systems with state-dependent diffusion,
# Environmetrics, 22 (8), 946{961, doi:10.1002/env.1135.


model.2res.ODE <- function(par, L, inp.file, out.data)
{
  # model:
  # ------
  #
  #   flow through a 2 cascade reservoir model:
  #
  #   dS1                          1
  #   ----  =  exp(A) * P + a0 -  -- * S1
  #    dt                          K
  #
  #   dS2        1   
  #   ----  = - --  * (S1-S2)
  #    dt        K
  #
  #             1                2*pi*t            4*pi*t            2*pi*t            4*pi*t
  #   flow  =  -- * S2 + si1*sin(------) + si2*sin(------) + co1*cos(------) + co2*cos(------)
  #             K                  24                24                24                24
  #
  # state variables:
  # ----------------
  #
  #   S1      cascade state 1
  #   S2      cascade state 2
  #
  # parameters:
  # -----------
  #
  #   A       log(catchment area)
  #   a0      mean dry weather flow
  #   K       time lag constant
  #   si1,si2,co1,co2     parameters of trigonometric function to describe dry weather variation
  
  start.t <- proc.time()
  
  pack.dir <- substr(system.file("extdata", "gawk.exe", package = Package),1 ,         # path to the package
                     nchar(system.file("extdata", "gawk.exe", package = Package))-9)
  
  # set original parameter values:
  # --------------
  par.orig <- c(A = log ( (1/3) * 1.3 * 1000 * 1000 ),   # log(m^2)
                k = 0.3,  # hr
                S1 = 4, S2 = 4,    # initial cascade states
                nsi1 = 0, nsi2 = 0, nco1 = 0, co2 = 0, a0 = 0
               )
  
  # mulitply the original parameters:
  # --------------
  param <- par[paste("mult.", names(par.orig), sep="")]
  param[ which(is.na(param)) ] <- 1
  param <- as.numeric(param)
  
  par.mod <- par.orig * param
  
  # get rain data:
  # --------------
  rain.path <- paste(substr(inp.file, 1, nchar(inp.file)-23), "data/",
                     substr(inp.file, nchar(inp.file)-22, nchar(inp.file)-4 ), 
                     "_rain.dat", sep = "")
  rain.all  <- read.csv(rain.path, sep=";", header=T)   
  
  Rain     <- rain.all[,(1:2)]
  if (length(rain.all[1,]) > 2) {
    Rain[,2] <- apply( rain.all[,2:length(rain.all[1,])], MARGIN = 1, mean  )   #// using the mean of available RGs  
  }
  
  Rain <- Rain[c(T,F),]   # every second row
  Rain[,2] <- Rain[,2] / 1000 # mm/h -> m/h
  
  # manage time layout:
  # --------------
  Dt = 2/60 # (2 min)
  
  origo   = kimisc::hms.to.seconds(format(as.POSIXct(Rain[1,1], format="%m/%d/%Y %H:%M:%S",tz="UCT"), format="%H:%M:%S"))/3600 #hr
  if ( sysanal.decode(L)$val[1] - origo <= 10/60) {
    origo = sysanal.decode(L)$val[1]
  }
  ult.cal = (nrow(Rain)-1)*Dt + origo #hr -- end of event 1
  t.grid  = format( seq(origo, ult.cal, Dt), nsmall=6 )
  Lm      = paste("Q", t.grid, sep="_")
  time    = sysanal.decode(Lm)$val
  
  Rain <- Rain[,2]
  
  
  co1 = -par.mod["nco1"]; co2 =  par.mod["co2"]; si1 = -par.mod["nsi1"]; si2 = -par.mod["nsi2"]
  WWF = si1*sin(2*pi*time/24) + si2*sin(4*pi*time/24) + co1*cos(2*pi*time/24) + co2*cos(4*pi*time/24)
  

  # STEP I: start simulations inustus_temp times earlier
  # burn-in phase
  #---------------
  inustus_temp = 40 # time steps        
  
  start = max(Dt, time[1]-inustus_temp*(Dt))
  end   = time[length(time)]
  time  = seq(start, end, by=Dt)
  if ((length(time)-inustus_temp) < length(Lm))
  {
    end = end + Dt*(length(Lm)-(length(time)-inustus_temp))
    time =seq(start,end, by=Dt)
  }
  
  Inp.elon = c(rep(0, inustus_temp), Rain) # m/h
  

  # calculate results:
  # ------------------  
  state  <-c(par.mod["S1"], par.mod["S2"])
  res_ode <- deSolve::ode(y = state, times = time, func = cascade2, parms = par.mod, Inp.elon = Inp.elon, tempo = time)
  
  #convert reservoir state to flow
  K = par.mod["k"]
  SWF = 1/K*res_ode[,"S2"]
  SWF = SWF[!is.element(SWF, NA)] # remove a particular group of elements (which are NAs) in a vector
  
  flow <- SWF[(inustus_temp+1):length(SWF)] / 3.6 + WWF   #// m^3/h -> l/s
  res <- as.numeric(flow)
  res[ which(res<0) ]  <- 1e-7
  
  shortLm <- substr(Lm, 1, nchar(Lm)-1 )
  shortL  <- substr(L , 1, nchar(L )-1 )
  matching <- match(shortL, shortLm)
  if ( length(which(is.na(matching)==TRUE)) == length(matching) ) {
    stop("timestamp mismatch")
  }
  
  yM <- as.numeric(as.character(res[matching])) # removes data from timesteps when there is NA in observed data
  names(yM) = L
  
  
  end.t <- proc.time(); time.taken <- end.t - start.t; print(time.taken)
  
  return(yM)
}


####################################################################################
#function defining the rhs of the ODE's
cascade2<-function(t, y, parms, Inp.elon, tempo) {
  with(as.list(c(y, parms)),{
    
    a0 <- parms["a0"]; K <- parms["k"]; A <- parms["A"]; 
    S1 <- y["S1"]; S2 <- y["S2"]
    
    start <- tempo[1]
    end   <- tempo[length(tempo)]
    
    if ( (t>start) && (t<end) ) {
      P <- Inp.elon[ min(which(tempo>=t)) - 1 ]
    } 
    else {
      P <- Inp.elon[1]
    }
    if (t>end) {
      P <- Inp.elon[length(tempo)]
    }
    
    # rate of change
    dS1 <- exp(A)*P + a0 - 1/K*S1
    dS2 <- 1/K*(S1-S2)
    
    # return the rate of change
    list(c(dS1, dS2))
  })
}

