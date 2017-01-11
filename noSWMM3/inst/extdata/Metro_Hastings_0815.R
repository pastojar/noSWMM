Metro_Hastings_0815 <-
function(li_func,pars,prop_sigma=NULL,iterations=50000,burn_in=1,adapt_par=c(100,100,0.5,0.45),quiet=FALSE,...)
{     
    if (!is.finite(li_func(pars, ...))) 
#       print(!is.finite(li_func(pars, ...)))
#       print(pars)
#       print(str(li_func(pars, ...)))
        stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")

    if(!is.null(dim(prop_sigma)))
    {
        if( ( dim(prop_sigma)[1] != length(pars) ||  dim(prop_sigma)[2] != length(pars) ) && !is.null(prop_sigma) )
            stop("prop_sigma not of dimension length(pars) x length(pars)")
    }

    if(is.null(prop_sigma)) #if no proposal matrix given, estimate the Fisher information to use as the diagonal (start in the right variance scale)
    {
        if(length(pars)!=1)
        {
            fit<-optim(pars,li_func,control=list("fnscale"=-1),method = "Nelder-Mead", hessian=TRUE,...)
            fisher_info<-solve(-fit$hessian)
            prop_sigma<-sqrt(diag(fisher_info))
            prop_sigma<-diag(prop_sigma)
        }else{
            prop_sigma<-1+pars/2
        }
    }

    prop_sigma <- MHadaptive::makePositiveDefinite(prop_sigma)
	mu<-pars
	pi_X<-li_func(pars,...)			        # initial likelihood evaluation
	k_X<-pars
	trace          <-array(dim=c(iterations,length(pars)))
  # -*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-  
  colnames(trace)    = names(pars)
  log.post           = array(dim=c(iterations,1))
  colnames(log.post) = "log.post"
   acc.opt           = .3
  # -*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*-  
    deviance<-array(dim=iterations) 
    announce<-floor(seq(iterations/10,iterations,length.out=10))
   nacc  <- 0
   npro  <- 0
#    nacc.t  <- 0
#    npro.t  <- 0
   
	for(i in 1:iterations)
	{
	  
	  npro  <- npro+1
# 	  npro.t  <- npro.t+1
	        k_Y<-MASS::mvrnorm(1,mu=k_X,Sigma=prop_sigma)	# Draw proposal point
        				
        	pi_Y<-li_func(k_Y,...)		        # evaluate likelihood at proposal point
        		
        	a_X_Y = (pi_Y)-(pi_X)		        # Compare relative likelihoods
        	if(!is.finite(a_X_Y))           # -*--*--*--*--*--*--*--*--*-
                    a_X_Y<--Inf                          # never jump outside the range of the parameter space
        	if( log(runif(1,0,1)) <= a_X_Y)	        # Make the jump according to MH probability
              	{
        		k_X = k_Y
        		pi_X = pi_Y
        		nacc  <- nacc+1
#         		nacc.t  <- nacc.t+1
        	}	
        
        	trace[i,]<-k_X	
          log.post[i,] <-(pi_X)   # -*--*--*--*--*--*--*--*--*-
#           print(paste("parameters", trace[i,], "posterior",log.post[i]))
                deviance[i]<-(-2*pi_X)                      # Store the deviance for calculating DIC
		if(i > adapt_par[1] && i %% adapt_par[2] == 0 && i < (adapt_par[4]*iterations) )	        # adapt the proposal covariance structure
		{   
#       browser()
                    acc.rr  <- nacc/npro
                    len<-floor(i*adapt_par[3]):i
                    x<-trace[len,]
                    N<-length(len)
#                     print("points considered"); print(x)                  
#                     acc.tot <-length(unique(trace[,1]))/(i-burn_in)
# #                     acc     <-length(unique(tail(trace[,1],n=N+length(iterations)-i) ))/(N)
#                     acc     <-length(unique(x))/(N); 
#                 
#                     print(paste("length(unique(x[,1])): ", length(unique(x[,1]))))
#                     print(paste("length unique(x): ", length(unique(x))))
#                     print(paste("length(unique(trace[,1])): ", length(unique(trace[,1]))))
#                     print(paste("i: ", i))
#                     print(paste("acc tot: ",acc.tot ))
#                     print(paste("acc considered period: ", format(acc, digits=2,scientific=F)))
#                     print(paste("real acc considered period: ", format(acc.rr, digits=2,scientific=F)))
#                             
                    c = acc.rr/acc.opt
#                     c = 1
#                     if (acc >1) stop("there is some problem")
#                     print("c");print(c)  
                    p_sigma <- (N-1) * var(x)/N
                    p_sigma <-c*MHadaptive::makePositiveDefinite(p_sigma)   # To deal with rounding problems that can de-symmetrize
                    c     <- 1
                    nacc  <- 0
                    npro  <- 0
                    if(!(0 %in% p_sigma) ) 
                        prop_sigma<-p_sigma
         }
         if(!quiet && i %in% announce)
            print(paste('updating: ------------------------------------------',i/iterations*100,'%',sep=''))
	}
    trace<-trace[burn_in:iterations,]
    DIC<-NULL
    ## Calculate the DIC
    D_bar<-mean(deviance[burn_in:iterations])
    if(length(pars)>1)
    {
        theta_bar<-sapply(1:length(pars),function(x){mean( trace[,x] )})
    }else
        theta_bar<-mean( trace )

    D_hat<-li_func(theta_bar,...)
    pD<-D_bar-D_hat
    DIC<-D_hat + 2*pD

        accept_rate<-length(unique(trace[,1]))/(iterations-burn_in) 
   
# print(paste("real acc tot: ", format(nacc.t/npro.t, digits=2,scientific=F)))
   
   
    val<-list("trace"=trace,"prop_sigma"=prop_sigma,"DIC"=DIC,'acceptance_rate'=accept_rate,
              "unnor.log.post"=log.post)
    class(val)<-"MHposterior"
	return(val)
}
