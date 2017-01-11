######################################################################################
# 1 cascade reservoir model
#
# J. Pastorek, DEC 2016
# based on a script by D. DelGiudice
######################################################################################

model.1res <- function(par, L, inp.file, out.data) 
{
  # model:
  # ------
  #
  #   flow through a 1 cascade reservoir model:
  #
  #   dS1                                                 S1
  #   ----  =  Qin - Qout ;     Qin= A*P + bf ;   Qout =  -- 
  #    dt                                                 K
  #
  #
  #     S1  = S1_0*exp(-Dt/K) + Qin[t]*K - Qin[t]*K*exp(-Dt/K)     #//  what is this exactely ???
  # 
  #                            2*pi*t            4*pi*t            2*pi*t            4*pi*t
  #   flow  =   Qout + si1*sin(------) + si2*sin(------) + co1*cos(------) + co2*cos(------)
  #                             24                24                24                24
  #
  # state variables:
  # ----------------
  #
  #   S1      cascade state 1
  #
  # parameters:
  # -----------
  #
  #   A       catchment area [m2]
  #   K       reservoir residence time [hr]
  #   si1,si2,co1,co2     parameters of trigonometric function to describe dry weather variation
  #   Rain    Inp [m/hr]
  #   bf      constant groundwater inflow in the sewer [m3/hr]
  #
  #   S1_0    initial cascade state  
  
  
  # start.t <- proc.time()
  
  pack.dir <- substr(system.file("extdata", "gawk.exe", package = Package),1 ,         # path to the package
                     nchar(system.file("extdata", "gawk.exe", package = Package))-9)
  
  # set original parameter values:
  # --------------
  par.orig <- c(A = (1/3) * 1.3 * 1000 * 1000,  # m2
                k = 0.3,  # hr
                S1_0 = 4,   # initial cascade state
                nsi1 = 0, nsi2 = 0, nco1 = 0, co2 = 0, bf = 0 
#                 nsi1 = 3.044435e-01, nsi2 = 8.289373e-01,
#                 nco1 = 6.874255e-01, co2 = 8.195015e-02, 
#                 bf = 2 #7.389576, #from l/s to m3/hr
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
  
  origo   = kimisc::hms.to.seconds(format(as.POSIXct(Rain[1,1], format="%m/%d/%Y %H:%M:%S",tz="UCT"), format="%H:%M:%S"))/3600#hr
  if ( sysanal.decode(L)$val[1] - origo <= 10/60) {
    origo   = sysanal.decode(L)$val[1]
  }
  ult.cal = (nrow(Rain)-1)*Dt + origo #hr -- end of event 1
  t.grid  = format( seq(origo, ult.cal, Dt), nsmall=6 )
  Lm      = paste("Q", t.grid, sep="_")
  time    = sysanal.decode(Lm)$val
  
  Rain <- Rain[,2]
  
  
  co1 = -par.mod["nco1"]; # if (is.na(co1)) print("par missing")
  co2 =  par.mod["co2"];   # if (is.na(co2)) print("par missing")
  si1 = -par.mod["nsi1"]; # if (is.na(si1)) print("par missing")
  si2 = -par.mod["nsi2"]; # if (is.na(si2)) print("par missing")
  
  WWF = si1*sin(2*pi*time/24)+si2*sin(4*pi*time/24)+co1*cos(2*pi*time/24)+co2*cos(4*pi*time/24)
  
  # STEP I: start simulations inustus_temp times earlier
  # burn-in phase
  #---------------
  
  inustus_temp = 20 # time steps
  
  start = max(Dt, time[1]-inustus_temp*(Dt))
  end   = time[length(time)]
  time  = seq(start, end, by=Dt)
  if ((length(time)-inustus_temp) < length(Lm))
  {
    end = end + Dt*(length(Lm)-(length(time)-inustus_temp))
    time =seq(start,end, by=Dt)
  }
  
  Rain.elon = c(rep(0, inustus_temp), Rain) # m/h
  
  # calculate results:
  # ------------------
  
  Qin  = (par.mod["A"] * Rain.elon)  + par.mod["bf"] # m3/h
  SWF <- Qout_Qin_dSdt(Qin = Qin, K = par.mod["k"], S_0 = par.mod["S1_0"], Dt = Dt, time = time) # m^3/h
  
  flow <- SWF[(inustus_temp+1):length(SWF)] / 3.6 + WWF #l/s
  res  <- as.numeric(flow)
  res[which(res<0)]  <- 1e-7 # eliminate neg val
  
  shortLm <- substr(Lm, 1, nchar(Lm)-1 )
  shortL  <- substr(L , 1, nchar(L )-1 )
  matching <- match(shortL, shortLm)
  if ( length(which(is.na(matching)==TRUE)) == length(matching) ) {
    stop("timestamp mismatch")
  }
  
  yM <- as.numeric(as.character(res[matching])) # removes data from timesteps when there is NA in observed data
  names(yM) = L
  
#   # saves model  parameters and ouptuts (Q)
#   # ------------------
#   out.name <- substr(inp.file, nchar(inp.file)-22, nchar(inp.file)-4)
#   hlp <- as.character(as.numeric(format(param, digits=4)))
#   out.name <- paste( out.name, hlp[1], sep="_")
#   for (i in 2:length(hlp) ) {
#     out.name <- paste(out.name, hlp[i], sep="")
#   }
#   par.out <- data.frame( names = names(par), values = unname( format(par, digits = 17, scientific = TRUE) ) )
#   write.table(par.out, paste(pack.dir,"/data/", out.name, "_par.out", sep=""), col.names=F, row.names=F, sep=";", quote=F)
#   
#   out.time <- kimisc::seconds.to.hms( time[(inustus_temp+1):length(SWF)]*3600 + 1 )
#   yM.out  <- data.frame( time = out.time , values = round(unname(res), digits=3) )
#   write.table(yM.out, paste(pack.dir,"/data/", out.name, "_Q.out", sep=""), col.names=F, row.names=F, sep=";", quote=F)
  
  
  # end.t <- proc.time(); time.taken <- end.t - start.t; print(time.taken)
  
  return(yM)
}


#-----------------------------------------------------------------------------
#
# Qout = Qin - dS/dt
#
Qout_Qin_dSdt <- function(Qin, K, S_0, Dt, time) {
  
  S   <- rep(NA,length(time))
  SWF <- rep(NA,length(time))
  for (t in 1:length(time))
  {
    S[t]  <- S_0*exp(-1/K*Dt) + Qin[t]*K - Qin[t]*K*exp(-1/K*Dt) # analytical sol linear reservoir
    
    SWF[t]<- Qin[t] - (S[t]-S_0)/Dt  
    
    S_0    <- S[t] 
  }
  
  return(SWF)
}


#-----------------------------------------------------------------------------
#
# Qout = S/K
#
Qout_S_K <- function(Qin, K, S_0, Dt, time) {
  
  S   <- rep(NA,length(time))
  SWF <- rep(NA,length(time))
  for (t in 1:length(time))
  {
    SWF[t]  <- S_0/K*exp(-1/K*Dt) + Qin[t] - Qin[t]*exp(-1/K*Dt)  
    
    S_0    <- SWF[t]*K 
  }
  
  return(SWF)
}