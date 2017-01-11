######################################
## J. Pastorek, JAN 2017
######################################



#######################################
## name of this package
Package <- "noSWMM3"  # stays only a global variable...
pack.dir <- substr(system.file("extdata", "gawk.exe", package = Package), 1,       # path to the package 
                   nchar(system.file("extdata", "gawk.exe", package = Package))-9)

#######################################
## path to the SWMM model
Urquell    <- system.file("extdata", "inpfile.inp", package = Package) 

#######################################
## reads FG overview, FG data and RG overview
uni.data <- read.uni.data(FG.ov.path  = system.file("extdata", "statistics_(I+mwl1+2+rg1+2)+FG.csv", package = Package),
                          FG.dat.path = system.file("extdata", "(I+mwl1+2+rg1+2)+FG_MP1_reg_2min.csv", package = Package),
                          RG.ov.2mm.path  = system.file("extdata", "info_(I+mwl1+2+rg1+2)+FG_2mm.csv", package = Package),
                          RG.ov.path      = system.file("extdata", "info_(I+mwl1+2+rg1+2)+FG.csv", package = Package) )

#######################################
## selects events to work with by their number; you can decide based on       View(uni.data$RG.overview.2mm[good.events,])
good.events <- setdiff(68:86, c(76, 81, 84)) # usable events
eventsCa  <-  c(69, 73, 79, 80, 86) # desired events for CALIBRATION 
set.seed(42)
eventsPre <-  setdiff(good.events, eventsCa)    # desired events for PREDICTION  

#######################################
## reads and prepares data

      # from all RGs
prodata <- read.RG.data(RG.dat.path = system.file("extdata", "data_(I+mwl1+2+rg1+2)+FG.csv", package = Package),
                        uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell)
#       # from the RG #1
# prodata <- read.RG.data(RG.dat.path = system.file("extdata", "data_(I+mwl1+2+rg1+2)+FG.csv", package = Package),
#                         uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell, whichRGs = 1)
# 
#       # from MWs
# prodata <- read.MW.data(MW.dat.path = system.file("extdata", "MWL_rain_JA.RData", package = Package),
#                         uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell)

  ### stage II ###
# data.sources <- c("mwl1","mwl2", "rg1", "rg2")
#       # mwl1 
# prodata <- read.MW.dataII(MW.dat.path = system.file("extdata", "mwl_data_Jaro_JA_rounded.RData", package = Package), 
#                           data.source = data.sources[1],
#                           uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell)
#       # mwl2 
# prodata <- read.MW.dataII(MW.dat.path = system.file("extdata", "mwl_data_Jaro_JA_rounded.RData", package = Package), 
#                           data.source = data.sources[2],
#                           uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell)
#       # rg1 
# prodata <- read.MW.dataII(MW.dat.path = system.file("extdata", "mwl_data_Jaro_JA_rounded.RData", package = Package), 
#                           data.source = data.sources[3],
#                           uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell)
#       # rg2 
# prodata <- read.MW.dataII(MW.dat.path = system.file("extdata", "mwl_data_Jaro_JA_rounded.RData", package = Package), 
#                           data.source = data.sources[4],
#                           uni.data = uni.data, eventsCa = eventsCa, eventsPre = eventsPre, Urquell = Urquell)



#######################################
## defines model and parameters (incliding their prior) to be used  

par.2res    <- data.frame(mult.A  = c("NormalTrunc", 1, 1, 0.1, 2 ),
                          #mult.S1_0 = c("NormalTrunc", 1, 1, 0.1, 5 ),
                          #mult.S2_0 = c("NormalTrunc", 1, 1, 0.1, 5 ),
                          mult.k1  = c("NormalTrunc", 1, 1, 0.1, 5 ),
                          mult.k2  = c("NormalTrunc", 1, 1, 0.1, 5 ) )


par.1res    <- data.frame(mult.A    = c("NormalTrunc", 1, 1, 0.1, 2 ),
                          #mult.S1_0 = c("NormalTrunc", 1, 1, 0.1, 5 ),
                          mult.k    = c("NormalTrunc", 1, 1, 0.1, 5 ) )

                           # when modifying, do not forget to change also parameters in modelSWMM and CaPre and .awk file !
par.SWMM     <- data.frame(#mult.imp = c("NormalTrunc", 1, 1, 0.8, 1.2 ),
                           #mult.wid = c("NormalTrunc", 1, 1, 0.3, 1.7 ),
                           #mult.slo = c("NormalTrunc", 1, 1, 0.3, 1.7 ),
                           mult.Nim = c("NormalTrunc", 1, 1, 0.3, 1.7 ),
                           mult.Sim = c("NormalTrunc", 1, 1, 0.3, 1.7 ),
                           #mult.Spe = c("NormalTrunc", 1, 1, 0.3, 1.7 ),
                           #mult.Pze = c("NormalTrunc", 1, 1, 0.3, 1.7 ),
                           mult.rou = c("NormalTrunc", 1, 1, 0.3, 1.7 ) )

#par.hydr <- par.SWMM
#par.hydr <- par.1res
par.hydr <- par.2res

par.hydr <- t(par.hydr); colnames(par.hydr) <- c("distr", "EV", "sd", "min", "max")


#model <- model.SWMM 
#model <- model.1res
model <- model.2res 


#######################################
## checks model functioning
group.run.plot(dir=pack.dir, name="2res knlkn", par=par.hydr[,"EV"], ev.data=prodata, model=model)

# hlp <- list()
# for (j in 1 : length(prodata$Ca)) {
#   hlp[[j]] <-  model.swmm(par = par, inp.file = prodata$Ca[[j]][[3]], 
#                           out.data = prodata$Ca[[j]][[2]], L = prodata$Ca[[j]][[1]])
# }
# hlp <- list()
# for (j in 1 : length(prodata$Pre)) {
#   hlp[[j]] <-  model.swmm(par = par, inp.file = prodata$Pre[[j]][[3]], 
#                           out.data = prodata$Pre[[j]][[2]], L = prodata$Pre[[j]][[1]])
# }

#######################################
## sets up output transformation
transf <- list( transf="BC", par.tr=c(l1=0.45, l2=1) )
#transf <- list( transf="BC",  par.tr=c(l1=0.35, l2=0) ) # DelGiudice (2013)
#transf <- list( transf="LogSinh",  par.tr=c(alpha=5, beta=100) ) # DelGiudice (2013)

#######################################
## calibration and predictions
runs <- 10*c(1, 90, 50, 50)   
seed <- 42

start_time <- proc.time()
  Ca.res  <- Ca (prodata = prodata, model = model, par.hydr = par.hydr, input = input, transf = transf, seed = seed, runs = runs)
end_time=proc.time(); time_taken=end_time-start_time; time_taken


check <- FALSE
check <- plot.Ca.res(Ca.res = Ca.res, pack.dir = pack.dir) # plots calibration results
if (check==FALSE) {dev.off()} # closes graphic device if plotting fails

runs[4] <- 100*2
Pre.res <- Pre(prodata = prodata, model = model, transf = transf, runs = runs, RAM = Ca.res$RAM, par.tr = Ca.res$par.tr, par.fix = Ca.res$par.fix) 

#######################################
## creates statistics of the inference results
subsets <- list(wo.ev.2 =  c(1, 3:11), 
                weak    =  c(1, 4, 6, 10, 11), 
                strong  =  c(3, 5, 7, 8, 9) 
               )
statistics <- list()
for (i in 1 : length(subsets)) {
  statistics[[ names(subsets)[[i]] ]] <- statist.CaPre.res(Pre.res = Pre.res, prodata = prodata, skip = setdiff(1:11, subsets[[i]]) )
} 
all.RGs <- list(Pre.res = Pre.res, statistics = statistics, transf=transf) # a list for plotting


#######################################           
## plots prediction results
to.plot.list   <- list(all.RGs = all.RGs, l2 = all.RGs, l3 = all.RGs)
check <- FALSE
check <- plot.Pre.res(pack.dir = pack.dir, prodata = prodata, to.plot.list = to.plot.list, Ca.res = Ca.res)
if (check==FALSE) {dev.off()} # closes graphic device if plotting fails





save.image(file = paste(pack.dir, "/all.RGs.Rdata", sep=""))
