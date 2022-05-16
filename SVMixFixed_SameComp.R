# General Function for Mixture to eliminate excess code

# Inputs: data product as string (e.g. "aphro" or "trmm")
#         number of mixture components, including point mass
#         value to use for vconst in fixing means/vars
#         number of total iterations to keep, burnin, and thin
#         write plots for sample distribution fits?
#         write trace plots?
SVMix <- function(product,ncomp,vconst,keep.its,burn.its,thin.its,plot.fits=TRUE,plot.trace=TRUE,num.cores=10L) {
  # load relevant libraries
  suppressMessages(library(invgamma))
  suppressMessages(library(ggplot2))
  suppressMessages(library(truncnorm))
  suppressMessages(library(tidyverse))
  suppressMessages(library(parallel))
  source("AMCMCUpdate.R")

  # Read in the data and pull off appropriate months
  # load all the data files to get overall quantiles
  load("./aphro.Rdata")
  load("./era5.Rdata")
  load("./merra2.Rdata")
  load("./trmm.Rdata")
  a.dat <- eval(as.name("aphro"))
  e.dat <- eval(as.name("era5"))
  m.dat <- eval(as.name("merra2"))
  t.dat <- eval(as.name("trmm"))
  dates.a <- dimnames(a.dat)[3] %>% unlist()
  dates.e <- dimnames(e.dat)[3] %>% unlist()
  dates.m <- dimnames(m.dat)[3] %>% unlist()
  dates.t <- dimnames(t.dat)[3] %>% unlist()
  summer.a <- substr(dates.a, nchar(dates.a)-1, nchar(dates.a)) %in% c("04","05","06","07","08","09")
  summer.e <- substr(dates.e, nchar(dates.e)-1, nchar(dates.e)) %in% c("04","05","06","07","08","09")
  summer.m <- substr(dates.m, nchar(dates.m)-1, nchar(dates.m)) %in% c("04","05","06","07","08","09")
  summer.t <- substr(dates.t, nchar(dates.t)-1, nchar(dates.t)) %in% c("04","05","06","07","08","09")
  a.dat <- a.dat[,,which(summer.a)]
  e.dat <- e.dat[,,which(summer.e)]
  m.dat <- m.dat[,,which(summer.m)]
  t.dat <- t.dat[,,which(summer.t)]

  # then just do this data that we are fitting to
  load(paste0("./",product,".Rdata"))
  data <- eval(as.name(product))
  dates <- dimnames(data)[3] %>% unlist()
  summer.months <- substr(dates, nchar(dates)-1, nchar(dates)) %in% c("04","05","06","07","08","09")
  data <- data[,,which(summer.months)]
  logdata <- log(data)
  n.T <- dim(data)[3]

  # Set number of mixture components to ncomp specified
  K <- ncomp

  ## Prior for mu.u
  u.m <- 0
  u.s2 <- 1
  # Priors for means and variances
  #mu.k.m <- rep(0, K-1)
  #mu.k.s2 <- rep(2, K-1)^2
  #sigma2.k.a <- rep(2.1, K-1)
  #sigma2.k.b <- 1*(sigma2.k.a-1)

  ## MCMC settings
  burn <- burn.its # want 400,000 draws this time (approx 11 days on rencher)
  n.it <- keep.its # iterations to keep
  thin <- thin.its
  kp <- 0
  tot.it <- n.it*thin + burn
  kp.seq <- seq(burn+thin, tot.it, by=thin)

  ## Matrices for MCMC Draws
  #mu.k.draws <- matrix(NA, nrow=n.it, ncol=K-1)
  #sigma2.k.draws <- matrix(NA, nrow=n.it, ncol=K-1)
  #loglike <- NULL
  umean.draws <- matrix(NA,nrow=n.it, ncol=prod(dim(data)[1:2]))
  cutpoints.draws <- matrix(NA,nrow=n.it, ncol=K-1)
  
  ## Proposal Variance for Cut Points
  cut.prop.sd <- 0.001#0.0005
  cut.amcmc <- list(mn=matrix(0, nrow=K-2, ncol=1),
                    var=matrix(0, nrow=K-2, ncol=K-2))
  amcmc.it <- 250

  ###
  ## Initialize Parameters
  ###

  ## Initialize component means and variances (on log scale since log-normal)
  #mu.k <- aggregate(logdata[-zero.locs], by=list(z=z[-zero.locs]), FUN=mean)[,2] #get mean and var of logdata belonging to each nonzero component
  #sigma2.k <- aggregate(logdata[-zero.locs], by=list(z=z[-zero.locs]), FUN=var)[,2]
  
  ## Instead of initializing, fixing means and variances:
  # mn.raw <- seq(min(data[data>0]),max(data[data>0]), length=K)
  mn.raw <- quantile(c(a.dat[a.dat>0],e.dat[e.dat>0],m.dat[m.dat>0],t.dat[t.dat>0]),probs=seq(0,1,length=K))
  sig2.raw <- ((diff(mn.raw)/2)/vconst)^2
  mn.raw <- mn.raw[-length(mn.raw)] + 0.5*diff(mn.raw)
  mu.k <- log((mn.raw^2)/sqrt(sig2.raw+mn.raw^2))
  sigma2.k <- log(1+sig2.raw/(mn.raw^2))
  
  ## Plot the Set of Log-Normal Mixture Components
  # xseq <- seq(min(data[data>0]), max(data[data>0]), length=250)
  # df <- expand.grid(x=xseq, grp=1:(K-1))
  # df <- within(df, {
  #   dln <- dlnorm(x, meanlog=mu.k[grp], sdlog=sqrt(sigma2.k[grp]))
  # })
  # ggplot(data=df, aes(x=x, y=dln, color=factor(grp))) + geom_line()

  # initialize Zs as whichever is most likely (if equals 0, is 1)
  z <- rep(1,prod(dim(data))) # Initialize all zs at 1
  z.nzmat <- mapply(function(mu, sig){dlnorm(data[data>0], mu, sig)},
              mu=mu.k, sig=sqrt(sigma2.k))
  z[data>0] <- apply(z.nzmat,1,which.max)+1 # because zero comp is 1 right now

  ## Initialize cut points
  cutpts <- qnorm(c(0, cumsum(prop.table(table(z)))), mean=0, sd=1) #Pick cutpoints based on proportions in each z
  cutpts.trans <- log(diff(cutpts[2:K]))
  
  ###
  ## Make list with item for each location containing loc specific params
  ###
  info <- vector('list',length=dim(data)[1]*(dim(data)[2]))
  locs <- expand.grid(dimnames(data)[[1]],dimnames(data)[[2]])
  for(i in 1:nrow(locs)) {
    this.data <- data[as.character(locs[i,1]),as.character(locs[i,2]),]
    info[[i]]$coord <- locs[i,]
    info[[i]]$umu <- 0
    z <- rep(1, length(this.data)) # Initialize all zs at 1
    zero.locs <- (this.data==0) #then replace all the zs not from the zero component
    z[!zero.locs] <- matrix(mapply(function(mns, vars){
      dlnorm(this.data[!zero.locs], meanlog=mns, sdlog=sqrt(vars), log=TRUE)
    }, mns=mu.k, vars=sigma2.k), ncol=(K-1)) %>% apply(., 1, which.max)+1
    info[[i]]$z <- z
    info[[i]]$u <- rep(NA,length=dim(data)[3])
    info[[i]]$precip <- this.data
    info[[i]]$ztab <- table(factor(z, levels=1:K))
    om <- diff(pnorm(cutpts, mean=info[[i]]$umu, sd=1))
    # lke <- mapply(function(mns, vars, mixprob){
    #   log(mixprob) + dlnorm(this.data[!zero.locs], 
    #                         meanlog=mns, sdlog=sqrt(vars),
    #                         log=TRUE)
    # }, mns=mu.k, vars=sigma2.k, mixprob=om[-1])
    # lke <- max(lke) + log(sum(exp(lke-max(lke))))
    # lke <- lke + sum(this.data==0)*log(om[1])
    # info[[i]]$llike <- lke
    info[[i]]$llike <- sum(info[[i]]$ztab[om>0]*log(om[om>0]))
  }

  ###
  ## Moran Basis Function Expansion for mu_u
  ###
  D <- fields::rdist(locs)
  A <- 1/D #Inverse-weighted Neighborhood Matrix (rather than 0/1)
  diag(A) <- 0
  P <- diag(nrow(A)) - matrix(1/nrow(A), nrow=nrow(A), ncol=ncol(A))
  E <- (P%*%A%*%P)  %>% eigen(., symmetric=TRUE)
  E$vectors <- E$vectors[,E$values>0]
  E$values <- E$values[E$values>0]
  # qplot(1:length(E$values), cumsum(E$values)/sum(E$values), geom="line")
  # ggplot() + geom_raster(mapping=aes(x=locs[,1], y=locs[,2], fill=E$vectors[,255])) +
  #   scale_fill_distiller(palette="Spectral")
  n.B <- which((cumsum(E$values)/sum(E$values))>0.50)[1]
  B <- cbind(1,E$vectors[,1:n.B])
  BtBinv <- (1/diag(crossprod(B)))
  BtBinvBt <- scale(B, center=FALSE, scale=1/BtBinv) %>% t(.)

  # smooth.mu.u <- B%*%BtBinvBt%*%mu.u
  # plt.mu.u <- ggplot() + geom_raster(mapping=aes(x=locs[,1], y=locs[,2], fill=mu.u)) +
  #   scale_fill_distiller(palette="Spectral")
  # plt.smooth.mu.u <- ggplot() + geom_raster(mapping=aes(x=locs[,1], y=locs[,2], fill=smooth.mu.u)) +
  #   scale_fill_distiller(palette="Spectral")
  # gridExtra::grid.arrange(plt.mu.u, plt.smooth.mu.u)
  # cor(mu.u, smooth.mu.u)^2

  ###
  ## Functions for Gibbs and MH Sampling
  ###

  #Function to sample zs for each location
  get.zs <- function(lst) {
    om <- diff(pnorm(cutpts, mean=lst$umu, sd=1))
    tweight <- matrix(mapply(function(meanlog, sdlog, wgt){
      dlnorm(lst$precip[lst$precip!=0], meanlog, sdlog)*wgt
      }, meanlog=mu.k, sdlog=sqrt(sigma2.k), wgt=om[-1]),ncol=(K-1))
    if (length(lst$precip[lst$precip!=0])!=0) {
      lst$z[lst$precip!=0] <- apply(tweight, 1, function(x){sample(2:K, 1, prob=x/sum(x))})}
    lst$ztab <- table(factor(lst$z, levels=1:K))
    lst$llike <- sum(lst$ztab[om>0]*log(om[om>0]))
    return(lst)
  }

  update.us <- function(lst){
    lst$u <- rtruncnorm(length(lst$z),a=cutpts[lst$z],b=cutpts[lst$z+1],
                        mean=lst$umu, sd=1)
    return(lst)
  }

  recalculate.umu <- function(ind){
   info[[ind]]$umu <- B[ind,]%*%umu.coef
   return(info[[ind]])
  }

  cps.llike <- function(lst, cps=cutpts) { 
    om <- diff(pnorm(cps, mean=lst$umu, sd=1))
    if(sum(lst$ztab[om==0])>0){
      return(-Inf)
    } else {
      num <- sum(lst$ztab[om>0]*log(om[om>0]))
      return(num-lst$llike)
    }
    # om <- diff(pnorm(cps, mean=lst$umu, sd=1))
    # lke <- mapply(function(mns, vars, mixprob){
    #   log(mixprob) + dlnorm(lst$precip[lst$precip>0], 
    #                         meanlog=mns, sdlog=sqrt(vars),
    #                         log=TRUE)
    # }, mns=mu.k, vars=sigma2.k, mixprob=om[-1])
    # lke <- max(lke) + log(sum(exp(lke-max(lke))))
    # lke <- lke + sum(lst$precip==0)*log(om[1])
    # return(lke-lst$llike)
  }

  ###
  ## Run MCMC Loop
  ###

  ## Progress Bar
  #pb <- txtProgressBar(min = 0, max = tot.it, style = 3)
  system.time(
  for(it in 1:tot.it){ # 
    
    ## Update cutpoints
    if(it>amcmc.it){
      prop.var <- (2.4^2/(K-2))*((cut.prop.sd^2)*diag(K-2)+cut.amcmc$var)
    } else {
      prop.var <- (cut.prop.sd^2)*diag(K-2)
    }
    a.prop <- cutpts.trans+t(chol(prop.var))%*%rnorm(length(3:K))
    c.prop <- c(cutpts[1:2],
                cutpts[2]+cumsum(exp(a.prop)),
                cutpts[K+1])
    mh.ratio <- sapply(info, cps.llike, cps=c.prop)# assumes uniform prior
    sum(mh.ratio)
    if(log(runif(1))<sum(mh.ratio)){
        cutpts <- c.prop
        cutpts.trans <- a.prop
    }
    cut.amcmc <- AMCMC.update(cutpts.trans, cut.amcmc$mn, cut.amcmc$var, it)
    
    # Sample zs
    info <- mclapply(info, FUN=get.zs,mc.cores=num.cores)

    #draw mu1-muk and sigma1-sigmak
    # ztab <- rowSums(sapply(info, function(lst){return(lst$ztab)}))
    #tpar <- sapply(2:K, function(x){
      #data.k <- sapply(info, function(lst){
        #return(log(lst$precip[lst$z==x]))
      #}) %>% unlist()
    
      ## Draw mu.k
      #mu.k.post.s2 <- (mu.k.s2[x-1]*sigma2.k[x-1])/(ztab[x]*mu.k.s2[x-1]+sigma2.k[x-1])
      #mu.k.post.m <- (sum(data.k)*(mu.k.s2[x-1])+mu.k.m[x-1]*sigma2.k[x-1])/
      #  (ztab[x]*mu.k.s2[x-1]+sigma2.k[x-1])
      #mdraw <- rnorm(1, mean=mu.k.post.m, sd=sqrt(mu.k.post.s2))
    
      ## Draw sigma2.k
      #sigma2.k.astar <- sigma2.k.a[x-1]+ztab[x]/2
      #sigma2.k.bstar <- sigma2.k.b[x-1]+0.5*sum((data.k-mdraw)^2)
      #s2draw <- rinvgamma(1, shape=sigma2.k.astar, 
      #                     rate=sigma2.k.bstar)
      #return(c(mdraw, s2draw))
    #})
    #mu.k <- tpar[1,]
    #sigma2.k <- tpar[2,]
  
    ## Update u's
    info <- lapply(info,FUN=update.us)
  
    ## Update mu.u
    ubar <- sapply(info, FUN=function(lst){
      return(mean(lst$u))
    })
    umu.coef <- BtBinvBt%*%ubar + rnorm(n.B+1, 0, sqrt(BtBinv/n.T))
    smooth.mu.u <- B%*%umu.coef
    info <- lapply(1:length(info), FUN=recalculate.umu)

    #calculate -2 log likelihood # NOT WORKING YET!!!
    # Numberzero*w1 + sum(numberklog(wk))+sumwherez!=0(dlnorm(dat,mu,sig,log=TRUE))

    #save values, sorting to keep them constant
    if(it%in%kp.seq){
      kp <- kp + 1
      #mu.k.draws[kp,] <- sort(mu.k)
      #sigma2.k.draws[kp,] <- sigma2.k[order(mu.k)]
      umean.draws[kp,] <- sapply(info, FUN=function(lst){return(lst$umu)})
      cutpoints.draws[kp,] <- cutpts[-c(1,(K+1))]
      #loglike[i] <- -2*sum(sapply(info,FUN=calc.lik)) #calculate likelihoods
    }
  
    # Save all of the values to a file to try to catch errors!
   # all.values <- list(info,a.prop,c.prop,mh.ratio,cutpts,cutpts.trans,cut.amcmc,ubar,umu.coef,smooth.mu.u)
   # save(all.values,file="AllCurrentValues.Rdata")

    ## Increment progress bar
    #setTxtProgressBar(pb, it)
  }
  #close(pb)
  )

  ####
  ## Posterior Convergence Assessment
  ####

  # Trace Plots
  if (plot.trace==TRUE) {
    pdf(paste0('Sample Trace ',product,' Fixed2.pdf'))
    #par(ask=TRUE)
    #for (i in 1:ncol(mu.k.draws)){
    #  plot(mu.k.draws[,i],type='l',main=paste("Mean",i))
    #}
    #for (i in 1:ncol(sigma2.k.draws)){
     # plot(sigma2.k.draws[,i],type='l',main=paste("Variance",i))
    #}
    for (i in 1:ncol(cutpoints.draws)){
      plot(cutpoints.draws[,i],type='l',main=paste("Cutpoint",i))
    }
    pick.locs <- sample(1:nrow(locs),size=K-1)
    for (i in pick.locs){
      plot(umean.draws[,i],type='l',main=paste("U mean",i))
    }
    dev.off()
  }

  if (plot.fits==TRUE) {
    # Distribution Check (Random Locations)
    pick.locs <- sample(1:nrow(locs),size=6)
    mn <- mu.k #apply(mu.k.draws,2,mean)
    sig2 <- sigma2.k #apply(sigma2.k.draws,2,mean)
    meancp <- apply(cutpoints.draws,2,mean)

    pdf(paste0('Sample Dist Fits ',product,'2.pdf'))
    for(i in pick.locs) {
      this.umu <- mean(umean.draws[,i],na.rm=TRUE)
      this.data <- info[[i]]$precip
      omega <- diff(pnorm(c(-Inf,meancp,Inf),mean=this.umu,sd=1))
      ztst <- sample(1:K, size=1000, replace=TRUE, prob=omega)
      vals <- rep(0, 1000)
      vals[ztst!=1] <- rlnorm(sum(ztst!=1), meanlog=mn[ztst[ztst!=1]-1],sdlog=sqrt(sig2[ztst[ztst!=1]-1]))
      plot(density(this.data,from=0))
      lines(density(vals,from=0),col='red')
    }
    dev.off()
  }

  #results.list <- list(mu=mu.k.draws,sigma2=sigma2.k.draws,umean=umean.draws,cutpoints=cutpoints.draws)
  results.list <- list(mu=matrix(rep(mu.k,keep.its),byrow=TRUE,nrow=keep.its),sigma2=matrix(rep(sigma2.k,keep.its),byrow=TRUE,nrow=keep.its),umean=umean.draws,cutpoints=cutpoints.draws)
  save(results.list,file=paste0('SameComp',K,product,'2.Rdata'))

  #DIC <- means['n2LogLike']+(1/2)*var(results$n2LogLike
}
