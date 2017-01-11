######################################
# J. Pastorek, FEB 2016
######################################


read.uni.data <- function(FG.ov.path, RG.ov.path, RG.ov.2mm.path, FG.dat.path) {
  
  # reads data
  FG.overview <- read.csv(FG.ov.path, sep=";", header=T)
  for (i in c(2,3,4,7,10,13)) {
    FG.overview[,i] <- as.POSIXct( FG.overview[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
  
  RG.overview.2mm <- read.csv(RG.ov.2mm.path, sep=";", header=T)
  for (i in c(1,2)) {
    RG.overview.2mm[,i] <- as.POSIXct( RG.overview.2mm[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
  
  RG.overview <- read.csv(RG.ov.path, sep=";", header=T)
  for (i in c(1,2)) {
    RG.overview[,i] <- as.POSIXct( RG.overview[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
  
  FG.data <- read.csv(FG.dat.path, sep=";", header=T)
  for (i in c(1,4)) {
    FG.data[,i] <- as.POSIXct( FG.data[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
 
  # selects FG events according to RG events  (more RG than FG data available)
  matchRGFG <- match(RG.overview.2mm$st, FG.overview$id)
  if ( length(which(is.na(matchRGFG))) != 0 ) {
    matchRGFG <- matchRGFG[ - which(is.na(matchRGFG)) ]
  }
  FG.overview$id[tail(matchRGFG, 1) + 1] <- FG.overview$en[tail(matchRGFG, 1)] +1
  FG.overview <- FG.overview[ - (tail(matchRGFG, 1) + 2 : length(FG.overview[,1])),]
  FG.sel.overview <- FG.overview[matchRGFG,]
  if (length(FG.sel.overview[,1]) != length(match(FG.sel.overview$id, RG.overview.2mm$st))) 
  { stop("RG FG data timestamp mismatch") }
  
  
  return(list(FG.sel.overview = FG.sel.overview, RG.overview.2mm = RG.overview.2mm, RG.overview = RG.overview, FG.data = FG.data))
}

######################################################################################
######################################################################################

read.RG.data <- function(RG.dat.path, uni.data, eventsCa, eventsPre, Urquell, whichRGs=3) {
 
  RG.data <- read.csv(RG.dat.path, sep=";", header=T)
  if (whichRGs == 1) {
    RG.data <- RG.data[,c(1, 2,3,4, 10)]
  }
  for (i in c(1,length(RG.data[1,]))) {
    RG.data[,i] <- as.POSIXct( RG.data[,i], origin="1970-01-01 00:00:00 UTC", tz="UTC" )
  }
  
  FG.sel.overview <- uni.data$FG.sel.overview; RG.overview.2mm <- uni.data$RG.overview.2mm; FG.data <- uni.data$FG.data
  
  prodata <- list(); prodata$Ca <- list(); prodata$Pre <- list()
  for (j in eventsCa) {
    prodata$Ca[[j]]  <- setupSWMM.RG(j=j, fg.sel.overview=FG.sel.overview, rg.overview=RG.overview.2mm, 
                                    fg.data=FG.data, rg.data=RG.data, urquell=Urquell)
  }
  for (j in eventsPre) {
    prodata$Pre[[j]] <- setupSWMM.RG(j=j, fg.sel.overview=FG.sel.overview, rg.overview=RG.overview.2mm, 
                                     fg.data=FG.data, rg.data=RG.data, urquell=Urquell)
  }
  
  j <- 0; dataCa <- list()
  for (i in 1 : length(prodata$Ca)) {
    if ( !is.null(prodata$Ca[[i]]) ) {
      j <- j+1
      dataCa[[j]] <- prodata$Ca[[i]]
    }
  }
  
  j <- 0; dataPre <- list()
  for (i in 1 : length(prodata$Pre)) {
    if ( !is.null(prodata$Pre[[i]]) ) {
      j <- j+1
      dataPre[[j]] <- prodata$Pre[[i]]
    }
  }
  
  prodata$Ca  <- dataCa
  prodata$Pre <- dataPre
  
  return(prodata)
}



######################################################################################
######################################################################################

read.MW.data <- function(MW.dat.path, uni.data, eventsCa, eventsPre, Urquell) {
  
  FG.sel.overview <- uni.data$FG.sel.overview; FG.data <- uni.data$FG.data
  RG.overview.2mm <- uni.data$RG.overview.2mm; RG.overview <- uni.data$RG.overview
  
  # reads MW data
  load(MW.dat.path)
  
  MW.data <- data.frame( time = tim1 )
  #MW.data <- cbind(MW.data, r2)
  MW.data <- cbind(MW.data, Rmwl.all)
  
  # matches with periods of interest and labels data with id
  st_en_ev <- data.frame( "st_ev" = RG.overview.2mm$st, "en_ev" = RG.overview.2mm$en)
  st_pos <- match(RG.overview.2mm$st, MW.data$time)  # where in the MW data periods of interest are
  NA_pos <- which(is.na(st_pos))
  st_en_ev <- st_en_ev[- NA_pos, ]  # deletes periods of interest when there is no MW data available
  st_pos   <- st_pos  [- NA_pos ]
  
  ev_lengths <- c()  
  step <- MW.data$time[2] - MW.data$time[1]
  MW.data$id <- NA
  
  for (i in 2:length(st_en_ev[,1])) {
    help <- seq(st_en_ev[i,1],st_en_ev[i,2], step)
    ev_lengths[i] <- length(help)  # remembers length of each event     
    MW.data$id[ st_pos[i] : ( st_pos[i] + ev_lengths[i] -1) ] <- rep( st_en_ev[i,1] , ev_lengths[i] )
    #print(i)
  }
  MW.data$id <- as.POSIXct( MW.data$id, origin="1970-01-01 00:00:00 UTC", tz = "UTC")
  
  
  prodata <- list(); prodata$Ca <- list(); prodata$Pre <- list()
  for (j in eventsCa) {
    prodata$Ca[[j]]  <- setupSWMM.MW(j=j, fg.sel.overview=FG.sel.overview, rg.overview=RG.overview.2mm, 
                                     fg.data=FG.data, rg.data=MW.data, urquell=Urquell)
    if (prodata$Ca[[j]][[1]] == 42) {stop(paste("invalid data in: prodata$Ca[",j, "]", sep=""))}
  }
  for (j in eventsPre) {
    prodata$Pre[[j]] <- setupSWMM.MW(j=j, fg.sel.overview=FG.sel.overview, rg.overview=RG.overview.2mm, 
                                     fg.data=FG.data, rg.data=MW.data, urquell=Urquell)
    if (prodata$Pre[[j]][[1]] == 42) {stop(paste("invalid data in: prodata$Pre[",j, "]", sep=""))}
  }
  
  j <- 0; dataCa <- list()
  for (i in 1 : length(prodata$Ca)) {
    if ( !is.null(prodata$Ca[[i]]) ) {
      j <- j+1
      dataCa[[j]] <- prodata$Ca[[i]]
    }
  }
  
  j <- 0; dataPre <- list()
  for (i in 1 : length(prodata$Pre)) {
    if ( !is.null(prodata$Pre[[i]]) ) {
      j <- j+1
      dataPre[[j]] <- prodata$Pre[[i]]
    }
  }
  
  prodata$Ca  <- dataCa
  prodata$Pre <- dataPre
  
  return(prodata)
}

######################################################################################
######################################################################################

read.MW.dataII <- function(MW.dat.path, data.source, uni.data, eventsCa, eventsPre, Urquell) {
  
  FG.sel.overview <- uni.data$FG.sel.overview; FG.data <- uni.data$FG.data
  RG.overview.2mm <- uni.data$RG.overview.2mm; RG.overview <- uni.data$RG.overview
  
  # reads MW data
  load(MW.dat.path)
  
  MW.data <- data.frame( time = t0 )
  MW.data <- cbind(MW.data, mwl1, mwl2, rg1, rg2)
  hlp <- c(which(names(MW.data) == "time"), which(names(MW.data) == data.source))
  MW.data <- MW.data[, hlp]
  
  # matches with periods of interest and labels data with id
  st_en_ev <- data.frame( "st_ev" = RG.overview.2mm$st, "en_ev" = RG.overview.2mm$en)
  st_pos <- match(RG.overview.2mm$st, MW.data$time)  # where periods of interest are in the MW data 
  NA_pos <- which(is.na(st_pos))
  st_en_ev <- st_en_ev[- NA_pos, ]  # deletes periods of interest when there is no MW data available
  st_pos   <- st_pos  [- NA_pos ]
  
  ev_lengths <- c()  
  step <- MW.data$time[2] - MW.data$time[1]
  MW.data$id <- NA
  
  for (i in 1:length(st_en_ev[,1])) {
    help <- seq(st_en_ev[i,1],st_en_ev[i,2], step)
    ev_lengths[i] <- length(help)  # remembers length of each event     
    MW.data$id[ st_pos[i] : ( st_pos[i] + ev_lengths[i] -1) ] <- rep( st_en_ev[i,1] , ev_lengths[i] )
    #print(i)
  }
  MW.data$id <- as.POSIXct( MW.data$id, origin="1970-01-01 00:00:00 UTC", tz = "UTC")
  
  
  prodata <- list(); prodata$Ca <- list(); prodata$Pre <- list()
  for (j in eventsCa) {
    prodata$Ca[[j]]  <- setupSWMM.MW(j=j, fg.sel.overview=FG.sel.overview, rg.overview=RG.overview.2mm, 
                                     fg.data=FG.data, rg.data=MW.data, urquell=Urquell)
    if (prodata$Ca[[j]][[1]] == 42) {stop(paste("invalid data in: prodata$Ca[",j, "]", sep=""))}
  }
  for (j in eventsPre) {
    prodata$Pre[[j]] <- setupSWMM.MW(j=j, fg.sel.overview=FG.sel.overview, rg.overview=RG.overview.2mm, 
                                     fg.data=FG.data, rg.data=MW.data, urquell=Urquell)
    if (prodata$Pre[[j]][[1]] == 42) {stop(paste("invalid data in: prodata$Pre[",j, "]", sep=""))}
  }
  
  j <- 0; dataCa <- list()
  for (i in 1 : length(prodata$Ca)) {
    if ( !is.null(prodata$Ca[[i]]) ) {
      j <- j+1
      dataCa[[j]] <- prodata$Ca[[i]]
    }
  }
  
  j <- 0; dataPre <- list()
  for (i in 1 : length(prodata$Pre)) {
    if ( !is.null(prodata$Pre[[i]]) ) {
      j <- j+1
      dataPre[[j]] <- prodata$Pre[[i]]
    }
  }
  
  prodata$Ca  <- dataCa
  prodata$Pre <- dataPre
  
  return(prodata)
}



######################################################################################
######################################################################################

setupSWMM.RG <- function(j, fg.sel.overview, rg.overview, fg.data, rg.data, urquell)  {
  
  # chooses event specified as j
  remain <- j %% length(fg.sel.overview$st)
  if (remain==0) {remain <- length(fg.sel.overview$st)}
  Ev.nmr <- remain
  Ev.id <- rg.overview$st[Ev.nmr]
  
  # selects data for the given event
  Ev.FG.pos  <- which( fg.data$id == Ev.id )
  Ev.FG.data <- fg.data[Ev.FG.pos,]
  
  Ev.RG.pos  <- which( rg.data$id == Ev.id )
  Ev.RG.data <- rg.data[Ev.RG.pos,]
  
  if (length(Ev.FG.data[,1]) < 3) {stop("event too short")}
  
  pack.dir <- substr(system.file("extdata", "gawk.exe", package = Package),1 , 
                     nchar(system.file("extdata", "gawk.exe", package = Package))-9)
  
  # process FG data
  # -------------------------
  Ev.FG.data  = Ev.FG.data[,-c(2,4)]
  Ev.FG.data[,2] = Ev.FG.data[,2]*1000   # m^3/s  ---> l/s
  
  Ev.FG.data[,1] = format(Ev.FG.data[,1], "%d/%m/%Y %H:%M")
  
  #Ev.FG.data  =Ev.FG.data[seq(1,nrow(Ev.FG.data),by=2),] # 4 minutes time step
  Ev.FG.data  = Ev.FG.data                                # 2 minutes time step now !
  colnames(Ev.FG.data) = c("Date.Time","Flow.L.s.")
  
  t = unlist(strsplit(as.character(Ev.FG.data[,1]), ":"))
  t = matrix(t, ncol = 2, byrow = T)
  if ( as.numeric(t[2,2])-as.numeric(t[1,2]) > 0 ) {
    dt.mod = (as.numeric(t[2,2])-as.numeric(t[1,2]))/60 # calculates time step of the flow data in hours
  } else {
    dt.mod = (as.numeric(t[3,2])-as.numeric(t[2,2]))/60 # calculates time step of the flow data in hours
  }
  
  
  
  # process RG data
  # -------------------------
  Ev.RG.data[,1] = format(Ev.RG.data[,1], "%m/%d/%Y %H:%M:%S")
  
  if (length(Ev.RG.data[1,]) == 10) { # data from all RGs
    Ev.RG.data[3] <- apply(Ev.RG.data[,c(3,6)], 1, mean, na.rm=T)  # avrg from 2 RGs at the same place  
    Ev.RG.data[4] <- apply(Ev.RG.data[,c(4,7)], 1, mean, na.rm=T)  # avrg from 2 RGs at the same place
    Ev.RG.data <- Ev.RG.data[,1:4]
    for (i in 2:4) {
      if ( length((which(is.na(Ev.RG.data[,i]) == TRUE))) != 0 ) {return(list(42, "dummy list entry", "dummy list entry"))}
    }
  }
  
  if (length(Ev.RG.data[1,]) == 5) { # data from 1 RG 
    Ev.RG.data[3] <- Ev.RG.data[2]              # data from the RG #1 for the whole catchment 
    Ev.RG.data[4] <- Ev.RG.data[2]              # data from the RG #1 for the whole catchment 
    Ev.RG.data <- Ev.RG.data[,1:4]
    for (i in 2:4) {
      if ( length((which(is.na(Ev.RG.data[,i]) == TRUE))) != 0 ) {return(list(42, "dummy list entry", "dummy list entry"))}
    }
  }
  
  # saves rain time series files and remembers the path
  out.name <- Ev.RG.data$time[1]
  hlp2 <- substr(out.name, 7, 10)
  out.name = sub( substr(out.name, 7, 10), substr(out.name, 1, 5), out.name)
  out.name = sub( substr(out.name, 1, 5), hlp2, out.name)
  out.name <- gsub("/", ".", out.name); out.name <- gsub(":", ".", out.name); out.name <- gsub(" ", "_", out.name)
  write.table(Ev.RG.data, paste(pack.dir, "/data/", out.name, "_rain.dat", sep=""), col.names=T, row.names=F, sep=";",quote=F)
  rain.path <- c()
  for (i in 1:3) {
    rain.path[i] <- paste(pack.dir, "/data/", out.name, "_rainRG", i, ".dat", sep="")
    write.table(Ev.RG.data[,c(1, 1+i)], rain.path[i],  col.names=F, row.names=F, sep=" ",quote=F) # which rain data to write
  }
  RG.path <- substr(rain.path[1], 1, nchar(rain.path[1])-5)
  
  
  # modifies inp files (sets proper time details and adds path to rain time series) 
  # and remembers the path to the inp file
  
  date = format(as.Date(Ev.FG.data[,1], format = "%d/%m/%Y"), format = "%m/%d/%Y")  
  time = format(as.POSIXct(Ev.FG.data[,1], format="%d/%m/%Y %H:%M",tz=""), format = "%H:%M:%S")
  
  START_DATE <- head(date,n=1)
  END_DATE   <- date[length(date)]
  
  START_TIME <- head(time,n=1)    
  START_TIME <- as.POSIXct(START_TIME, format = "%H:%M:%S", tz="UTC")
  START_TIME <- as.numeric(START_TIME) - 2*2*60 # 2*2 mins before FG data starts     
  START_TIME <- format(as.POSIXct(START_TIME, origin="1970-01-01 00:00:00 UTC", tz="UTC"), format = "%H:%M:%S")
  
  END_TIME  <- time[length(time)]   
  END_TIME  <- as.POSIXct(END_TIME, format = "%H:%M:%S", tz="UTC")
  END_TIME1 <- END_TIME + 2*60    # 2 mins after FG data ends
  daydif    <- (  as.numeric(format(END_TIME1, format = "%Y%m%d"))  
                  - as.numeric(format(END_TIME, format = "%Y%m%d")) 
  )
  
  if (daydif != 0)    # if the event lasts over more days
  {
    END_DATE <- as.POSIXct(END_DATE, format = "%m/%d/%Y", tz="UTC")
    END_DATE <- as.numeric(END_DATE) + daydif * 60*60*24
    END_DATE <- format(as.POSIXct(END_DATE, origin="1970-01-01 00:00:00 UTC", tz="UTC"), format = "%m/%d/%Y")
  }
  
  END_TIME <- format(END_TIME1, format = "%H:%M:%S")
  
  
  t.sim     = c(START_DATE, START_TIME, END_DATE, END_TIME, WET_STEP="00:02:00", REPORT_STEP="00:02:00")  # report step 2 mins
  
  inp.path <- paste(pack.dir, "/", out.name, ".inp", sep="")
  shell( paste(system.file("extdata", "gawk.exe", package = Package),
               " -v START_DATE=", t.sim[1],
               " -v START_TIME=", t.sim[2],
               " -v END_DATE=", t.sim[3],
               " -v END_TIME=", t.sim[4], 
               " -v WET_STEP=", t.sim[5],
               " -v REPORT_STEP=", t.sim[6],
               " -v RAIN_PATH=", RG.path,  # name of the time series has to be "RAINDAT" (else change time_modif.awk)
               " -f ", system.file("extdata", "time_modif.awk", package = Package),  
               " ", urquell, " > ",
               inp.path,
               sep="")
  )
  
  # creates time layouts
  origo     = kimisc::hms.to.seconds(format(as.POSIXct(Ev.FG.data[1,1], format="%d/%m/%Y %H:%M",tz=""), format="%H:%M:%S"))/3600 #[h]
  ult.cal   = (nrow(Ev.FG.data)-1)*dt.mod + origo # end of event 1 [h]
  t.grid = format(seq(origo,ult.cal,dt.mod), nsmall=6)
  L      = paste("Q",t.grid,sep="_")
  
  #4.2.16#
  if (length(which(is.na(Ev.FG.data[,2]))) > 0) { 
    L <- L[- which(is.na(Ev.FG.data[,2]))] 
    Ev.FG.data <- Ev.FG.data[- which(is.na(Ev.FG.data[,2])),]
  }
  #4.2.16#
  
  return( list(L, Ev.FG.data, inp.path) )
}


######################################################################################
######################################################################################

setupSWMM.MW <- function(j, fg.sel.overview, rg.overview, fg.data, rg.data, urquell)  {
  
  # chooses event specified as j
  remain <- j %% length(fg.sel.overview$st)
  if (remain==0) {remain <- length(fg.sel.overview$st)}
  Ev.nmr <- remain
  Ev.id <- rg.overview$st[Ev.nmr]
  
  # selects data for the given event
  Ev.FG.pos  <- which( fg.data$id == Ev.id )
  Ev.FG.data <- fg.data[Ev.FG.pos,]
  
  Ev.RG.pos  <- which( rg.data$id == Ev.id )
  Ev.RG.data <- rg.data[Ev.RG.pos,]
  
  if (length(Ev.FG.data[,1]) < 3) {stop("event too short")}
  
  pack.dir <- substr(system.file("extdata", "gawk.exe", package = Package),1 , 
                     nchar(system.file("extdata", "gawk.exe", package = Package))-9)
  
  # process FG data
  # -------------------------
  Ev.FG.data  = Ev.FG.data[,-c(2,4)]
  Ev.FG.data[,2] = Ev.FG.data[,2]*1000   # m^3/s  ---> l/s
  
  Ev.FG.data[,1] = format(Ev.FG.data[,1], "%d/%m/%Y %H:%M")
  
  #Ev.FG.data  =Ev.FG.data[seq(1,nrow(Ev.FG.data),by=2),] # 4 minutes time step
  Ev.FG.data  = Ev.FG.data                                # 2 minutes time step now !
  colnames(Ev.FG.data) = c("Date.Time","Flow.L.s.")
  
  t = unlist(strsplit(as.character(Ev.FG.data[,1]), ":"))
  t = matrix(t, ncol = 2, byrow = T)
  if ( as.numeric(t[2,2])-as.numeric(t[1,2]) > 0 ) {
    dt.mod = (as.numeric(t[2,2])-as.numeric(t[1,2]))/60 # calculates time step of the flow data in hours
  } else {
    dt.mod = (as.numeric(t[3,2])-as.numeric(t[2,2]))/60 # calculates time step of the flow data in hours
  }
  
  
  
  # process RG data
  # -------------------------
  Ev.RG.data[,1] = format(Ev.RG.data[,1], "%m/%d/%Y %H:%M:%S")
  colnames(Ev.RG.data) <- c("time", "MW.avrg") # mean over all MWLs
  Ev.RG.data <- Ev.RG.data[,1:2]
  for (i in 2:2) {
    if ( length((which(is.na(Ev.RG.data[,i]) == TRUE))) != 0 ) {return(list(42, "dummy list entry", "dummy list entry"))}
  }
  
  # saves rain time series files and remembers the path
  out.name <- Ev.RG.data$time[1]
  hlp2 <- substr(out.name, 7, 10)
  out.name = sub( substr(out.name, 7, 10), substr(out.name, 1, 5), out.name)
  out.name = sub( substr(out.name, 1, 5), hlp2, out.name)
  out.name <- gsub("/", ".", out.name); out.name <- gsub(":", ".", out.name); out.name <- gsub(" ", "_", out.name)
  write.table(Ev.RG.data, paste(pack.dir, "/data/", out.name, "_rain.dat", sep=""), col.names=T, row.names=F, sep=";",quote=F)
  rain.path <- c()
  for (i in 1:3) {
    rain.path[i] <- paste(pack.dir, "/data/", out.name, "_rainRG", i, ".dat", sep="")
    #write.table(Ev.RG.data[,c(1, 1+i)], rain.path[i],  col.names=F, row.names=F, sep=" ",quote=F) # which rain data to write
    write.table(Ev.RG.data[,c(1, 2  )], rain.path[i],  col.names=F, row.names=F, sep=" ",quote=F) # which rain data to write
  }
  RG.path <- substr(rain.path[1], 1, nchar(rain.path[1])-5)
  
  
  # modifies inp files (sets proper time details and adds path to rain time series) 
  # and remembers the path to the inp file
  
  date = format(as.Date(Ev.FG.data[,1], format = "%d/%m/%Y"), format = "%m/%d/%Y")  
  time = format(as.POSIXct(Ev.FG.data[,1], format="%d/%m/%Y %H:%M",tz=""), format = "%H:%M:%S")
  
  START_DATE <- head(date,n=1)
  END_DATE   <- date[length(date)]
  
  START_TIME <- head(time,n=1)    
  START_TIME <- as.POSIXct(START_TIME, format = "%H:%M:%S", tz="UTC")
  START_TIME <- as.numeric(START_TIME) - 2*2*60 # 2*2 mins before FG data starts     
  START_TIME <- format(as.POSIXct(START_TIME, origin="1970-01-01 00:00:00 UTC", tz="UTC"), format = "%H:%M:%S")
  
  END_TIME  <- time[length(time)]   
  END_TIME  <- as.POSIXct(END_TIME, format = "%H:%M:%S", tz="UTC")
  END_TIME1 <- END_TIME + 2*60    # 2 mins after FG data ends
  daydif    <- (  as.numeric(format(END_TIME1, format = "%Y%m%d"))  
                  - as.numeric(format(END_TIME, format = "%Y%m%d")) 
  )
  
  if (daydif != 0)    # if the event lasts over more days
  {
    END_DATE <- as.POSIXct(END_DATE, format = "%m/%d/%Y", tz="UTC")
    END_DATE <- as.numeric(END_DATE) + daydif * 60*60*24
    END_DATE <- format(as.POSIXct(END_DATE, origin="1970-01-01 00:00:00 UTC", tz="UTC"), format = "%m/%d/%Y")
  }
  
  END_TIME <- format(END_TIME1, format = "%H:%M:%S")
  
  
  t.sim     = c(START_DATE, START_TIME, END_DATE, END_TIME, WET_STEP="00:02:00", REPORT_STEP="00:02:00")  # report step 2 mins
  
  inp.path <- paste(pack.dir, "/", out.name, ".inp", sep="")
  shell( paste(system.file("extdata", "gawk.exe", package = Package),
               " -v START_DATE=", t.sim[1],
               " -v START_TIME=", t.sim[2],
               " -v END_DATE=", t.sim[3],
               " -v END_TIME=", t.sim[4], 
               " -v WET_STEP=", t.sim[5],
               " -v REPORT_STEP=", t.sim[6],
               " -v RAIN_PATH=", RG.path,  # name of the time series has to be "RAINDAT" (else change time_modif.awk)
               " -f ", system.file("extdata", "time_modif.awk", package = Package),  
               " ", urquell, " > ",
               inp.path,
               sep="")
  )
  
  # creates time layouts
  origo     = kimisc::hms.to.seconds(format(as.POSIXct(Ev.FG.data[1,1], format="%d/%m/%Y %H:%M",tz=""), format="%H:%M:%S"))/3600 #[h]
  ult.cal   = (nrow(Ev.FG.data)-1)*dt.mod + origo # end of event 1 [h]
  t.grid = format(seq(origo,ult.cal,dt.mod), nsmall=6)
  L      = paste("Q",t.grid,sep="_")
  
  #4.2.16#
  if (length(which(is.na(Ev.FG.data[,2]))) > 0) { 
    L <- L[- which(is.na(Ev.FG.data[,2]))] 
    Ev.FG.data <- Ev.FG.data[- which(is.na(Ev.FG.data[,2])),]
  }
  #4.2.16#
  
  return( list(L, Ev.FG.data, inp.path) )
}

