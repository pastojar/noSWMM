#########################################
# J. Pastorek, APR 2016
# based on a script by O. Wani
#
# MCMC Bayesian calibration, prediction and uncertainty analysis
#########################################


###################################################################################################################################

Ca <- function(prodata, model, par.hydr, input, transf, seed, runs) {
  
  dataCa <- prodata$Ca
  
  ##############################################################################
  ## Defines error model parameters
  
  par.hydr.init <- as.numeric(par.hydr[,"EV"]); names(par.hydr.init) <- rownames(par.hydr)   # hydrological model parameters defined elsewhere
  par.init  <- c(par.hydr.init,   
                 sd.Eps_Q = 0.5,   
                 sd.B_Q   = 50,    # based on anaylsis of previous data
                 corrlen  = 0.5
                ) 
  
  par.fix   <- c(Del.Max = 0, 
                 ks_Q    = 0,      #// we don't use these
                 Delta   = 0
                 )
  lim.sx.ks_Q = NA
  
  pri_sd <- 0.5 * par.init[ (length(par.hydr.init)+1) : length(par.init)]
  
  
  ##############################################################################
  ## Tranforms error model parameters
  
  ref.out   <- c(y.ref_Q = 1200)
  par.tr    <- transf$par.tr
  
  var="Q"
  
  (if (is.na(par.tr["alpha"])) 
  {  # BC or no transformation
    
    if (!is.na(pri_sd["corrlen"])) {
      par.init[paste("sd.B",var, sep="_")] <- par.init[paste("sd.B",var, sep="_")]*
                                              sysanal.boxcox.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                   par.tr["l1"], par.tr["l2"])
    }  
    
    par.init[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*
                                              sysanal.boxcox.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                   par.tr["l1"], par.tr["l2"])
    
    if (!is.na(par.init[paste("ks",var, sep="_")])) { # we have input dependence via an heteroskedastik EM
      lim.sx.ks_Q = as.numeric(sysanal.boxcox(0,par.tr["l1"],par.tr["l2"]) / max(imber.vel))
      pri_sd["ks_Q"] = (-lim.sx.ks_Q+(sysanal.boxcox(par.init["ks_Q"],par.tr["l1"],par.tr["l2"]) / max(imber.vel)))
      par.init["ks_Q"] = 0.1*pri_sd["ks_Q"]  # start from a small value
    } 
  }  
  
  else  
  { # log-sinh
    
    if (!is.na(pri_sd["corrlen"])) {
      par.init[paste("sd.B",var, sep="_")] <- par.init[paste("sd.B",var, sep="_")]*
                                              sysanal.logsinh.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                    par.tr["alpha"],par.tr["beta"])
    }
    
    par.init[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*
                                              sysanal.logsinh.deriv(ref.out[paste("y.ref",var, sep="_")], 
                                                                    par.tr["alpha"],par.tr["beta"])

    if (!is.na(par.init[paste("ks",var, sep="_")])) {
      lim.sx.ks_Q = as.numeric(sysanal.logsinh(0,par.tr["alpha"],par.tr["beta"]) / max(imber.vel))
      pri_sd["ks_Q"] = (-lim.sx.ks_Q+(sysanal.logsinh(par.init["ks_Q"],par.tr["alpha"],par.tr["beta"]) / max(imber.vel)))
      par.init["ks_Q"] = 0.1*pri_sd["ks_Q"]
    } 
  } 
  )  
  
  
  
  ##############################################################################
  ## Defines prior distributions of the error model parameters 
  
  
  pri_sd[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*50 #// !!!
  pri_sd[paste("sd.B",  var, sep="_")] <- par.init[paste("sd.B",  var, sep="_")]*50 #// !!!
  
  pri_min <- pri_sd; pri_max <- pri_sd
  pri_min[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*0.01    #// !!!
  pri_min[paste("sd.B",  var, sep="_")] <- par.init[paste("sd.B"  ,var, sep="_")]*0.00001 #// !!!
  pri_max[paste("sd.Eps",var, sep="_")] <- par.init[paste("sd.Eps",var, sep="_")]*150     #// !!!
  pri_max[paste("sd.B",  var, sep="_")] <- par.init[paste("sd.B"  ,var, sep="_")]*1e+08   #// !!!
  
  par.hyp    <- data.frame(sd.Eps_Q =  c("NormalTrunc", 
                                                par.init["sd.Eps_Q"], pri_sd["sd.Eps_Q"], 
                                                pri_min[paste("sd.Eps",var, sep="_")], pri_max[paste("sd.Eps",var, sep="_")] ), 
                           sd.B_Q   =  c("NormalTrunc", 
                                                par.init["sd.B_Q"], pri_sd["sd.B_Q"],
                                                pri_min[paste("sd.B",  var, sep="_")],  pri_max[paste("sd.B",  var, sep="_")] ),
                           corrlen  =  c("NormalTrunc", par.init["corrlen"], pri_sd["corrlen"], 0.01, 3 ), # same units as layout
                           ks_Q     =  c("NormalTrunc", lim.sx.ks_Q, pri_sd["ks_Q"], lim.sx.ks_Q, 1e+08 )
                           )
  par.hyp <- t(par.hyp); colnames(par.hyp) <- c("distr", "EV", "sd", "min", "max")
  
  other.distr.types <- list(delta    = unname( c("Exponential", par.init["Delta"]) )                            )
  
  par.NormTrunc <- rbind(par.hydr, par.hyp)  # prior distributions of the hydrological model parameters defined elsewhere
  prior.pbdis <- list()          
  for (i in 1 : length(par.NormTrunc[,1]) ) {
    prior.pbdis[[i]] <- unname(par.NormTrunc[i,])
    names(prior.pbdis)[i] <- rownames(par.NormTrunc)[i]
  }
  
  prior.pbdis <- c(prior.pbdis,  other.distr.types)
  
  # print(prior.pbdis)
  
  
  ##############################################################################
  ## Defines objective function
  
  logposterior.unlim.swmm <- function(par)
  { 
    names(par) <- names(par.init)
    out <- sysanal.logposterior.unlim.swmm( par,
                                            model         = model,
                                            dataCa        = dataCa, 
                                            prior.dist    = "indep",
                                            prior.def     = prior.pbdis,
                                            loglikeli     = sysanal.loglikeli.bias.inp.JA,
                                            par.fix       = par.fix,
                                            par.tr        = par.tr,
                                            Var.Bs        = sysanal.Var.Bs,
                                            sd.Eps        = sysanal.sd.Eps.L
                                           )
    
    if (rnorm(1, mean = 0, sd = 1) > 1.9) { # to monitor the progress during iterations
      print(paste("log post: ", format(out, digits=2)))
      print(par)
    }
    return(out)
  }
  
  
  # Checks ability to compute posterior
  
  logposterior.unlim.swmm(par = par.init)
  
  
  # Prepares calibration
  
  nlogposterior <- function(par) {
    out = -1*logposterior.unlim.swmm(par)
    return(out)
  }
  
  
  ############################################################################## 
  ## RUNS CALIBRATION in 2 steps
  
  set.seed(seed)
  
  # i. Optimization
  par.names <- intersect( names(par.init) , names(par.NormTrunc[,1]) )
  low_ran = as.numeric(par.NormTrunc[par.names,"min"]) * 1.01
  up_ran  = as.numeric(par.NormTrunc[par.names,"max"]) * 0.99

  Opt.precal <- GenSA::GenSA( par = par.init,
                              fn    = nlogposterior, 
                              lower =  low_ran,
                              upper =  up_ran,
                              control= list( max.call = runs[1], verbose=T)
                            )

  par.optim.1  <- Opt.precal$par;  names(par.optim.1) <- names(par.init)


  # ii. Improves jump distribution and keeps optimizing
  
  RAM      <- adaptMCMC::MCMC( p     = logposterior.unlim.swmm,
                               init  = par.optim.1*rnorm(n=length(par.init),mean=1,sd=0.01),
                               scale = diag((pri_sd/2)^2,length(par.init)),
                               n     = runs[2] + runs[3], 
                               adapt = runs[2],    
                               acc.rate = 0.3, # maybe bit smaller (0.27)           use gamma ?
                               n.start = 100 # maybe larger
                             )
  
  
  # Time taken so far
  end_time=proc.time()
  time_taken=end_time-start_time
  time_taken
  
  return(list(pr.dis=prior.pbdis, Opt.precal=Opt.precal, par.optim.1=par.optim.1, RAM=RAM,
              par.tr=par.tr, par.fix=par.fix))
}




###################################################################################################################################




Pre <- function(prodata, model, transf, runs, RAM, par.tr, par.fix) {
  
  dataCa  <- prodata$Ca
  dataPre <- prodata$Pre
  
  ##############################################################################
  ## Predictions...   (with uncertainty propagation)
  
  if (runs[4] < runs[3] ) {
    MCMC.propa <- RAM$samples[ sample( seq((nrow(RAM$samples) - runs[3] + 1), nrow(RAM$samples), by=1), runs[4] ),  ]
  } 
  if (runs[4] == runs[3]) {
    MCMC.propa <- RAM$samples[seq((nrow(RAM$samples) - runs[3] + 1), nrow(RAM$samples), 1), ]
  }
  
  # Calibration Phase
  
  res.swmm.LCa <- list()
  for (i in 1 : length(dataCa)) {
    res.swmm.LCa[[i]] <- CaPre.predict.Ca(evdata = dataCa[[i]], model = model, MCMC.propa = MCMC.propa, par.tr = par.tr, par.fix = par.fix)
  }
  
  
  # Validation Phase (predicts for Pre events)
  
  res.swmm.LPre <- list()
  for (i in 1 : length(dataPre)) {
    res.swmm.LPre[[i]] <- CaPre.predict.Pre(evdata = dataPre[[i]], model = model, MCMC.propa = MCMC.propa, par.tr = par.tr, par.fix = par.fix)
  }
  
  
  
  
  #########################################################
  ## Back Transform

  #  Ca events
  
  bTr.Ca <- list()
  for (i in 1 : length(dataCa)) {
    bTr.Ca[[i]] <- CaPre.bTr.Ca(transf = transf, res.swmm.LCa = res.swmm.LCa[[i]], L.Ca = dataCa[[i]][[1]])
  }
  
  
  #  Pre events
  
  bTr.Pre <- list()
  for (i in 1 : length(dataPre)) {
    bTr.Pre[[i]] <- CaPre.bTr.Pre(transf = transf, res.swmm.LPre = res.swmm.LPre[[i]], L.Pre = dataPre[[i]][[1]])
  }

  
 
  
  return(list(MCMC.propa=MCMC.propa, bTr.Ca=bTr.Ca, bTr.Pre=bTr.Pre))

}

###################################################################################################################################

