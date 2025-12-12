basket_aLasso<-function(y,n,lambda,gamma){
  if(length(y)<2){
    stop("y should be a vector of length 2 or more")
  } 
  
  if(length(y)!=length(n)){
    stop("The vectors y and n should have the same length")
  }
  
  if(lambda<0 | gamma<0){
    stop("lambda and gamma should be non-negative")
  }
  
  K<- length(y)
  I<-diag(rep(1,K))
  
  # Fit an adaptive lasso model for each of the K baskets
  pEst<-rep(NA,K)
  res.penFactor<-diag(K)
  res.penFactor[res.penFactor==0]<-NA
  cd<-rep(NA,K)

  for(l in 1:K){
    
    # When the number of responses is 0 or n, we cannot use the MLE to 
    # provide an initial estimate of the parameters for the adaptive lasso model 
    # For these cases, we will just set the estimate to be the MLE (no borrowing)
    if(y[l]==0){
      pEst[l]<-0
      # EPS No.31
      # Store weight 0 in places other than NA,
      # so that weights are not overwritten unless the number of responses is 0 or n
      res.penFactor[l, -l][is.na(res.penFactor[l, -l])]<-0
      cd[l]<-1
    } else if(y[l]==n[l]) {
      pEst[l]<-1
      res.penFactor[l, -l][is.na(res.penFactor[l, -l])]<-0 # EPS No.31
      cd[l]<-1
    } else {
      
      # Use the MLE as an initial estimate of the parameters 
      penFactor<- abs(logit((y[-l]/n[-l]))-logit((y[l]/n[l])))
      penFactor[abs(penFactor)==Inf]<-999
      # If the MLE of the parameter is 0, drop that parameter from the adaptive lasso model
      include<-which(penFactor!=0)
      penFactor<-penFactor[include]
      penFactor<-(penFactor^(-gamma))
      
      XNew<-I[,-l]
      XNew<-XNew[,include]
      
      if(length(include)==0){
        fit<-glm(cbind(y,n-y)~1,family=binomial())
      } else if(length(include)==1){
        # glmnet rescales the penaltyFactor argument internally, so we will 
        # scale the lambda to counter this
        lambdaNew<- lambda*sum(penFactor)/length(include)
        # glmnet requires at least 2 variables. To get around this, include a 
        # second dummy intercept. This dummy variable will always have a 
        # coefficient of zero. 
        fit <- glmnet(cbind(1,XNew), cbind(n-y,y), family = "binomial", 
                      intercept=T, 
                      standardize=F, 
                      lambda=lambdaNew, 
                      penalty.factor = c(0,penFactor))
      } else {
        # glmnet rescales the penaltyFactor argument internally, so we will 
        # scale the lambda to counter this
        lambdaNew<- lambda*sum(penFactor)/length(include)
        fit <- glmnet(XNew, cbind(n-y,y), family = "binomial", intercept=T, 
                      standardize=F, 
                      lambda=lambdaNew, 
                      relax = TRUE,
                      penalty.factor = penFactor)
      }
      
      if (fit$npasses <= 1e+05) {
        pEst[l]<-inv.logit(coef(fit)[1])
        cd[l]<-1
      } else {
        pEst[l]<-y[l]/n[l]
        cd[l]<-0
      }

      # EPS No.31
      res.penFactor[l, -l][include]<-penFactor
      res.penFactor[-l, l][include]<-penFactor
    }
  }
  
  list(pEst=pEst, penFactor=res.penFactor, cd=cd)
}

BCHM_distance<-function (nDat, xDat, mu = 0.2, sigma02 = 10, sigmaD2 = 0.001, 
                         alpha = 1e-60, d0 = 0.05, alpha1 = 50, beta1 = 10, tau2 = 0.1, 
                         phi1 = 0.1, deltaT = 0.05, thetaT = 0.6, burnIn = 10000, 
                         MCIter = 20000, MCNum = 20000, seed = 1000) 
{
  numArm <- length(nDat)
  if (numArm != length(xDat)) {
    stop("Numbers of subgroups in nDat and xDat are not equal.")
  }
  if (numArm > 20) {
    stop("Numbers of subgroups is more than 20.")
  }
  set.seed(seed)
  weight <- nDat
  delta <- deltaT
  alphaP <- alpha
  alpha <- alpha1
  beta <- beta1
  posi <- xDat
  phi2 <- phi1 + delta
  rate <- c()
  priorMean <- mu
  priorVar <- sigma02
  res <- posi/weight
  x <- as.matrix(res)
  estVarGroup <- sigmaD2
  inSD <- sd(res)
  result <- gibbsSampler(x, alphaP, priorMean, priorVar, estVarGroup, 
                         weight, burnIn, MCIter)
  tables <- result$tables
  sm <- matrix(0, numArm, numArm)
  tSize <- dim(tables)[1]
  for (i in 1:tSize) {
    rr <- tables[i, ]
    for (j in 1:numArm) {
      for (k in 1:numArm) {
        if (rr[j] == rr[k]) {
          sm[j, k] <- sm[j, k] + 1
        }
      }
    }
  }
  sm <- sm/tSize
  sm=sm
}

bin2dec<-function(bits){
  dec<-rep(0,nrow(bits))
  for(i in 1:nrow(bits)){
    for(j in 1:ncol(bits)){
      dec[i]<-dec[i]+bits[i,j]*2^(j-1)
    }
  }
  dec
}

create_data<-function(Nbask,Npat,true.pi,base,seed,cntn,anal.time){
  set.seed(seed)#different seed for each scenario
  #not using this function, but is applicable promptly
  
  if(sum(cntn)>0){
    Ncntn<-sum(cntn)
    
    tp<-true.pi*cntn
    tp<-tp[tp>0]
    if(Npat$method=="random"||Npat$method=="realloc"){
      ratio<-Npat$ratio*cntn
      ratio<-ratio[ratio>0]
      
      if(anal.time==1)Nmax<-Npat$max[1]
      else Nmax<-Npat$max[anal.time]-Npat$max[anal.time-1]
      
      for(i in 1:100){
        burnin<-rmultinom(n=1,size=Nmax,prob=ratio)
      }
      n<-c(rmultinom(n=1,size=Nmax-Ncntn,prob=ratio))+rep(1,Ncntn)
    }else if(Npat$method=="fix"){
      n<-Npat$nfix
    }
    
    for(i in 1:100){
      burnin<-rbinom(n=Ncntn,size=n,prob=tp)
    }
    x<-rbinom(n=Ncntn,size=n,prob=tp)
  }else{
    n<-0
    x<-0
  }
  
  k<-1
  aftr<-data.frame(x=rep(0,Nbask),n=rep(0,Nbask))
  for(i in 1:Nbask){
    if(cntn[i]){
      aftr$n[i]<-n[k]
      aftr$x[i]<-x[k]
      k<-k+1
    }
  }
  
  data<-data.frame(i=1:Nbask,x=0,n=0)
  data$x<-base$x+aftr$x
  data$n<-base$n+aftr$n
  
  data
}

create_lump_data<-function(Nbask,Npat,true.pi,Nsim,Nanal){
  if(length(true.pi)!=Nbask)stop("Different number of cohort and Nbask\n")
  
  n.tab<-list()
  length(n.tab)<-Nanal
  x.tab<-list()
  length(x.tab)<-Nanal
  
  for(j in 1:Nanal){
    
    if(Npat$method=="fix"){
      if(j==1) n.tab[[j]]<-matrix(Npat$nfix[j,],ncol=Nbask,nrow=Nsim,byrow=T)
      else{
        n.tab[[j]]<-matrix(Npat$nfix[j,],ncol=Nbask,nrow=Nsim,byrow=T)-n.tab[[j-1]]
      }
    }
    else if(Npat$method=="random"|Npat$method=="poisson"){
      if (Npat$method=="random") {
        if(j==1)Nmax<-Npat$max[1]
        else Nmax<-Npat$max[j]-Npat$max[j-1]
      
        n.tab[[j]]<-matrix(NaN,ncol=Nbask,nrow=Nsim)
        for(i in 1:Nsim){
          n.tab[[j]][i,]<-c(1+rmultinom(n=1,size=Nmax-Nbask,prob=Npat$ratio))
        }
      } else if (Npat$method=="poisson") {
        po.tab <- list()
        length(po.tab) <- Nsim

        for (i in 1:Nsim) {
          po.tab[[i]] <- matrix(NaN, ncol = 24, nrow = Nbask)
          for (k in 1:Nbask) {
            po.tab[[i]][k, ] <- rpois(24, 0.5)
          }
        }
        n.tab[[j]] <- matrix(NaN,ncol=Nbask,nrow=Nsim)

        if (j == 1) {
          for (i in 1:Nsim) {
            if (Nanal == 1) {
              n.tab[[j]][i, ] <- rowSums(po.tab[[i]][, 1:24])
            } else {
              n.tab[[j]][i, ] <- rowSums(po.tab[[i]][, 1:12])
            }
          }
        } else {
          for (i in 1:Nsim) {
            n.tab[[j]][i, ] <- rowSums(po.tab[[i]][, 1:24]) - rowSums(po.tab[[i]][, 1:12])
          }
        }
      }
    }
    x.tab[[j]]<-matrix(NaN,ncol=Nbask,nrow=Nsim)
    for(i in 1:Nsim){
      x.tab[[j]][i,]<-rbinom(size=n.tab[[j]][i,],prob=true.pi,n=Nbask)
    }
    
  }
  #n.tab[[2]] is a matrix to add newly n to monitor in stage 2
  
  list(n=n.tab,x=x.tab)
}

create_nk<-function(seed,now,Npat.rand){
  set.seed(seed+now)#
  nk<-rep(0,length(Npat.rand$ratio))
  while(any(nk==0)){
    nk<-rmultinom(n=1,size=Npat.rand$max-length(Npat.rand$ratio),prob=Npat.rand$ratio)
  }
  c(nk)#at least nk=1
}

create_pi<-function(Nbask,scen,pi.null,pi.alt){
  rem1<-as.numeric(Nbask%%3==1)
  rem2<-as.numeric(Nbask%%3==2)
  remn0<-as.numeric(Nbask%%3!=0)
  div3f<-floor(Nbask/3)
  div3c<-ceiling(Nbask/3)
  div2f<-floor(Nbask/2)
  div2c<-ceiling(Nbask/2)
  div3.fl2<-floor(floor(Nbask/3)/2)
  div3.ce2<-ceiling(floor(Nbask/3)/2)
  if(scen<=Nbask) pi<-c(rep(pi.null,Nbask-scen),rep(pi.alt,scen))
  else if(scen==7){
    if(pi.null==0.05&&pi.alt==0.20){
      pi<-c(rep(pi.null,div3f+rem2),rep(pi.null+0.10,div3f+rem1),rep(pi.alt,div3f+rem2))
    }
    else  pi<-c(rep(pi.null,div3f+rem2),rep((pi.null+pi.alt)/2,div3f+rem1),rep(pi.alt,div3f+rem2))
  }
  else if(scen==8)pi<-c(rep(pi.null,div3f+rem2),
                        rep(pi.alt+(pi.alt-pi.null)/2,div3f+rem1),
                        rep(pi.alt,div3f+rem2))
  else if(scen==9)pi<-c(rep(pi.null,div3f+rem2),
                        rep(pi.alt+(pi.alt-pi.null),div3f+rem1),
                        rep(pi.alt,div3f+rem2))
  else if(scen==10)pi<-c(rep(pi.null/2,div3.fl2),rep(pi.null,div3.ce2),
                         rep(pi.alt+(pi.alt-pi.null)/2,div3f+remn0),
                         rep(pi.alt,div3f+rem2))
  else if(scen==11)pi<-c(rep(pi.null,div3f+rem2),
                         rep(pi.null/2,div3f+rem1),
                         rep(pi.alt,div3f+rem2))
  else if(scen==12){
    if(pi.null==0.05&&pi.alt==0.20&&Nbask==4){
      pi<-c(pi.null,pi.null+0.10,pi.alt,pi.alt+0.10)
    }
    else  pi<-c(rep(pi.null/2,div3.fl2),rep(pi.null,div3.ce2+remn0),
                rep((pi.null+pi.alt)/2,div3.fl2),rep((pi.null+3*pi.alt)/4,div3.ce2),
                rep(pi.alt,div3.fl2+rem2),rep(pi.alt+(pi.alt-pi.null)/2,div3.ce2))
  }
  else if(scen==13)pi<-c(rep(pi.null/2,div3f+rem2),rep(pi.null/2,div3f-1),
                         rep(pi.null,div3f+remn0),
                         pi.alt)
  else if(scen==14)pi<-c(rep(pi.null,div3f+div3.fl2+rem2),
                         rep(pi.alt+(pi.alt-pi.null)/2,div3.ce2+rem1),
                         rep(pi.alt,div3f+rem2))
  else if(scen==15)pi<-c(rep(pi.null,Nbask-1),
                         rep((pi.alt+pi.null)/2,1))
  else if(scen==16)pi<-c(rep(pi.null,div2c),
                         rep((pi.alt+pi.null)/2,div2f))
  else if(scen==17)pi<-c(rep(pi.null,Nbask-1),
                         rep(pi.alt+(pi.alt-pi.null),1))
  else if(scen==18)pi<-c(rep(pi.null,div2c),
                         rep(pi.alt+(pi.alt-pi.null),div2f))
  else if(scen==19)pi<-c(rep(pi.null,div2c),
                         rep(pi.alt+(pi.alt-pi.null)/2,div2f))
  else if(scen==20){
    if(pi.null==0.05&&pi.alt==0.20){
      pi<-c(rep(pi.null,Nbask-1),rep(pi.alt+0.10,1))
    }
    else  pi<-c(rep(pi.null,Nbask-1),rep(pi.alt+(pi.alt-pi.null)/2,1))
  }
  else if(scen==21)pi<-c(rep(pi.null,div2f),
                         rep(pi.null+0.10,div2c))
  else if(scen==22)pi<-c(rep(pi.alt-0.05,Nbask-3),pi.alt,pi.alt+0.05,pi.alt+0.10)
  else if(scen==23)pi<-c(rep(pi.null,Nbask-2),pi.alt+0.10,pi.alt)
  else if(scen==24)pi<-c(pi.null,rep(pi.alt-0.05,Nbask-1))
  else if(scen==25)pi<-c(rep(pi.null,Nbask-2),pi.null+0.10,pi.alt)
  else if(scen==26)pi<-c(pi.null,rep(pi.alt+0.10,Nbask-1))
  else if(scen==27)pi<-c(rep(pi.null,ceiling((Nbask-2)/2)),rep(pi.null+0.15,floor((Nbask-2)/2)),rep(pi.alt,2))
  else if(scen==28)pi<-seq(pi.null-0.05,pi.alt+0.05,by=0.05)
  else if(scen==29)pi<-c(rep(pi.null,ceiling((Nbask)/2)),rep(pi.null+0.10,3),rep(pi.alt+0.05,1))
  else if(scen==29)pi<-c(rep(pi.null,ceiling((Nbask)/2)),rep(pi.null+0.10,3),rep(pi.alt+0.05,1))
  else if(scen==30)pi<-c(rep(pi.null,ceiling((Nbask)/2)),rep(pi.null+0.05,2),pi.null+0.15,pi.alt+0.10)
  else if(scen==31)pi<-c(rep(pi.null,ceiling((Nbask)/2)-1),rep(pi.null+0.10,ceiling((Nbask)/2)),pi.alt-0.05)
  else if(scen==32)pi<-c(rep(pi.null,2),rep(pi.null+0.15,3),rep(pi.alt,2),pi.alt+0.10)
  else if(scen==33)pi<-c(pi.null,rep(pi.null+0.10,floor((Nbask)/2)-1),rep(pi.null+0.15,floor((Nbask)/2)-1),pi.alt)
  else if(scen==34)pi<-c(rep(pi.null-0.05,Nbask-1),pi.alt)
  else if(scen==35)pi<-c(rep(pi.null,floor((Nbask)/2)-1),rep(pi.null+0.05,2),pi.null+0.10,rep(pi.alt,2))
  else if(scen==36)pi<-c(rep(pi.null-0.05,2),rep(pi.null,2),rep(pi.alt-0.05,2),rep(pi.alt+0.10,2))
  else if(scen==37)pi<-c(pi.null,rep(pi.null+0.05,floor((Nbask)/2)-1),rep(pi.null+0.10,floor((Nbask)/2)-1),pi.alt)
  else if(scen==38)pi<-c(rep(pi.null-0.05,floor(Nbask/3)),rep(pi.null,floor(Nbask/2)),pi.alt+0.05)
  else if(scen==39)pi<-c(pi.null-0.05,rep(pi.alt+0.10,floor(Nbask/3)),rep(pi.alt,floor(Nbask/2)))
  else if(scen==40)pi<-c(rep(pi.alt-0.10,floor(Nbask/3)),rep(pi.alt,floor(Nbask/3)),rep(pi.alt+0.10,floor(Nbask/3)))
  else if(scen==41)pi<-c(rep(pi.null-0.05,floor(Nbask/3)),rep(pi.null,floor(Nbask/3)),pi.alt,pi.alt+0.05)
  else if(scen==42)pi<-c(rep(pi.null-0.05,floor(Nbask/2)),rep(pi.null,floor(Nbask/3)),pi.alt-0.05)
  else if(scen==43)pi<-c(pi.null-0.05,pi.null,pi.null+0.10,pi.alt-0.10,pi.alt,pi.alt+0.10)
  else if(scen==44)pi<-c(pi.null-0.05,rep(pi.null,2),pi.alt-0.10,pi.alt,pi.alt+0.10)
  else if(scen==45)pi<-c(rep(pi.null,div3f+rem2),rep(pi.alt-0.10,div3f+rem1),rep(pi.alt,div3f+rem2))
  else if(scen==46)pi<-c(rep(pi.null-0.05,div3f+rem2),rep(pi.alt+(pi.alt-pi.null)/2,div3f+rem1),rep(pi.alt,div3f+rem2))
  else if(scen==47)pi<-c(rep(pi.null-0.05,div3f+rem2),rep(pi.alt+0.10,div3f+rem1),rep(pi.alt,div3f+rem2))
  else if(scen==48)pi<-c(rep(pi.null-0.05,div3f+rem2),rep(pi.alt+0.10,div3f+rem1),rep(pi.alt-0.05,div3f+rem2))
  else if(scen==49)pi<-c(rep(pi.null-0.05,div3f+rem2),pi.alt-0.10,pi.alt,rep(pi.alt+0.10,div3f+rem2))
  else if(scen==50)pi<-c(pi.null-0.05,pi.null,(pi.null+pi.alt)/2,pi.alt-0.10,pi.alt,pi.alt+0.10)
  else if(scen==51)pi<-c(pi.null-0.05,pi.null,pi.alt-0.10,pi.alt-0.05,pi.alt,pi.alt+0.10)
  else if(scen==101)pi<-c(rep(pi.null,1),
                          rep(pi.alt+(pi.alt-pi.null),Nbask-1))
  else if(scen==102)pi<-c(rep(pi.null,floor((Nbask-2)/2)),
                          rep(pi.alt,ceiling((Nbask-2)/2)),
                          rep(pi.alt+(pi.alt-pi.null),2))
  else if(scen==103)pi<-c(rep(pi.null,Nbask-2),
                          rep(pi.alt,1),
                          rep(pi.alt+(pi.alt-pi.null),1))
  else if(scen>=1000){#considering different ns
    scen.bin<-dec2bin(num=scen-1000,digit=Nbask)
    pi<-rep(pi.alt,Nbask)*scen.bin+rep(pi.null,Nbask)*!scen.bin
  }
  pi
}

create_pi_random<-function(Nbask,seed,now,pi.null,pi.alt){
  set.seed(seed+now)
  Neff<-sample(0:(Nbask-1),1)
  Aeff<-c(rep(FALSE,Neff),rep(TRUE,Nbask-Neff))
  pi<-c()
  for(i in 1:Nbask){
    if(Aeff[i])pi[i]<-runif(n=1,min=(pi.null+pi.alt)/2,max=0.50)
    else pi[i]<-pi.null
  }
  list(Aeff=Aeff,pi=pi)
}

dec2bin <- function(num, digit=10){
  bin.tab<-matrix(0,nrow=length(num),ncol=digit)
  for(j in 1:length(num)){
    bin<-c()
    tmp<-num[j]
    for(i in digit:1){
      bin[i]<-tmp%/%2^(i-1)
      tmp<-tmp-bin[i]*2^(i-1)
    }
    bin.tab[j,]<-bin
  }
  bin.tab
}

dir_name<-function(sim,prefix="none"){
  if(sim$Npat$method=="random"){
    
    Npat.name<-as.character(paste(c("_N",sim$Npat$max[1]),collapse=""))
    if(length(sim$Npat$max)>1){
      for(i in 2:length(sim$Npat$max)){
        Npat.name<-as.character(paste(c(Npat.name,"-",sim$Npat$max[i]),collapse=""))
      }
    }
  }else if(sim$Npat$method=="realloc"){
    
    Npat.name<-as.character(paste(c("_N",sim$Npat$max[1]),collapse=""))
    if(length(sim$Npat$max)>1){
      for(i in 2:length(sim$Npat$max)){
        Npat.name<-as.character(paste(c(Npat.name,"-",sim$Npat$max[i]),collapse=""))
      }
    }
    Npat.name<-as.character(paste(c(Npat.name,"ra"),collapse=""))
  }else if(sim$Npat$method=="fix"){
    
    Npat.name<-as.character(paste(c("_n",sim$Npat$nfix[1,1]),collapse=""))
    if(nrow(sim$Npat$nfix)>1){
      for(i in 2:nrow(sim$Npat$nfix)){
        Npat.name<-as.character(paste(c(Npat.name,"-",sim$Npat$nfix[i,1]),collapse=""))
      }
    }
  }else if (sim$Npat$method=="poisson"){
    Npat.name<-"_Nanal1_po"
    if(length(sim$Npat$max)>1){
      Npat.name<-as.character(paste(c("_Nanal",length(sim$Npat$max), "_po"),collapse=""))
    }
  }
  
  pn.name<-as.character(paste(c("0.",formatC(sim$pi.null*100,width=2,flag="0")),collapse=""))
  pa.name<-as.character(paste(c("0.",formatC(sim$pi.alt*100,width=2,flag="0")),collapse=""))
  
  if(prefix!="none")pf.name<-as.character(paste(c("result_eps/",prefix,"_K"),collapse=""))
  else pf.name<-"result_eps/K"
  
  dir<-as.character(paste(c(pf.name,sim$Nbask,Npat.name,"_p",pn.name,"-",pa.name,"_s",sim$seed),collapse=""))
  dir
}

drop_data<-function(Nbask,data,stop){
  if(length(stop)==0){
    data<-cbind(data,data.frame(j=1:Nbask))
    data
  }
  else{
    data<-cbind(data[-stop,],data.frame(j=1:(Nbask-length(stop)))) 
    data
  }
}

find_dir<-function(prt.dir,dir){
  #find prt.dir
  tmp<-c("./../",prt.dir)
  name<-paste(tmp,collapse="")
  if(file.exists(name)){
    cat("The parent directory exists.")
  }
  else{
    cat("There is no parent directories. They will being made...")
    dir.create(name)
  }
  
  #find dir
  tmp<-c("./../",prt.dir,"/",dir)
  name<-paste(tmp,collapse="")
  if(file.exists(name)){
    cat("The directory exists.", "\n")
  }
  else{
    cat("There is no directories. They will being made...", "\n")
    dir.create(name)
  }
}

get_theta0 = function(pi0, J){
  
  if(length(pi0) != 1 & length(pi0) != J){
    stop("Length of pi0 must either be 1 or match with the length of n and r.")
  }
  
  theta0 = logit(pi0)
  theta0[theta0 < -10] = -10
  theta0[theta0 > 10] = 10

  if(length(pi0) == 1){
    theta0 = rep(theta0, J)
  }

  return(theta0)
}

get.part <- function(R, max_cl){
  ## generate all the possible partitions 
  ## and store them in a matrix 
  part_mat <- t(setparts(R))
  part_mat <- part_mat[apply(part_mat, 1, function(x){
    length(unique(x))<=max_cl
  }), ]
  part_mat <- data.frame(part_mat)
  names(part_mat) <- LETTERS[1:R]
  #part_mat[-1, ]
  part_mat
}

H_distance<-function(K,nik,rik,q1,q0){
  D=matrix(NA,K,K)
  epsilon = 3*(q1-q0)/K
  for (j in 1:(K-1))
  {
    for (k in (j+1):K)
    {
      if ((rik[j]/nik[j])==(rik[k]/nik[k]))
      {
        rik[k]=rik[k] + epsilon
        epsilon<-epsilon+3*(q1-q0)/K
      }
    }
  }
  for (j in 1:K)
  {
    for (k in j:K)
    {
      if (j == k)
      {
        D[j,k] = 0
      }
      if (j != k)
      {
        a1=1+rik[j]
        a2=1+rik[k]
        b1=1+nik[j]-rik[j]
        b2=1+nik[k]-rik[k]
        D[j,k] = D[k,j] = sqrt(1-beta((a1+a2)/2,(b1+b2)/2)/sqrt(beta(a1,b1)*beta(a2,b2)))
      }
    }
  }
  D
}

JS_distance<-function(K,a,b,nk,xk){
  p <- seq(0, 1, length=1000)
  
  for(i in 1:K){
    if(i==1)beta<-matrix(dbeta(p, a+xk[i], b+(nk[i]-xk[i]))/sum(dbeta(p, a+xk[i], b+(nk[i]-xk[i]))),nrow=1)
    else beta<-rbind(beta,matrix(dbeta(p, a+xk[i], b+(nk[i]-xk[i]))/sum(dbeta(p, a+xk[i], b+(nk[i]-xk[i]))),nrow=1))
  }
  
  distance<-matrix(ncol=K,nrow=K)
  for(i in 1:K){
    for(j in i:K){
      if(i==j)distance[i,j]<-0
      else{
        distance[i,j]<-suppressMessages(JSD(rbind(c(beta[i,]),c(beta[j,])),unit = "log"))
        distance[j,i]<-distance[i,j]
      }
    }
  }
  distance
}

KL_distance2<-function(K,nik,rik,q1,q0,a,b){
  D=matrix(NA,K,K)
  for (j in 1:(K-1))
  {
    for (k in (j+1):K)
    {
      if ((rik[j]/nik[j])==(rik[k]/nik[k]))
      {
        rik[k]=rik[k]
      }
    }
  }
  for (j in 1:K)
  {
    for (k in j:K)
    {
      if (j == k)
      {
        D[j,k] = 0
      }
      if (j != k)
      {
        cdf1 = pbeta(seq(0,1,1/100),a+rik[j],b+nik[j]-rik[j])
        freqs1 = sapply(1:(length(cdf1)-1),FUN=function(x){cdf1[x+1]-cdf1[x]}) + 0.0001
        cdf2 = pbeta(seq(0,1,1/100),a+rik[k],b+nik[k]-rik[k])
        freqs2 = sapply(1:(length(cdf2)-1),FUN=function(x){cdf2[x+1]-cdf2[x]}) + 0.0001
        D[j,k] = D[k,j] = (KL.plugin(freqs1, freqs2)+KL.plugin(freqs2, freqs1))/2
      }
    }
  }
  D
}

MEM_distance<-function(data,Nbask,a,b,mem.pri){
  N.omega.j<-2^(Nbask-1)
  num.elem<-(Nbask^2-Nbask)/2
  omega.all<-dec2bin(1:2^num.elem-1,digit=num.elem)
  omega.eachk<-dec2bin(1:N.omega.j-1,digit=Nbask-1)
  
  beta.sub<-beta(a+data$x,b+data$n-data$x)
  
  omega.j<-data.frame()
  k<-1
  for(i in 1:Nbask){
    for(j in 1:N.omega.j){
      omega.obs<-rep(0,Nbask)
      omega.obs[i]<-1
      omega.obs[-i]<-as.integer(intToBits(j-1))[1:(Nbask-1)]
      if(k==1){
        omega.j<-omega.obs
      }else{
        omega.j<-rbind(omega.j,omega.obs)
      }
      k<-k+1
    }
  }
  sum.tmp<-as.data.frame(as.matrix(omega.j)%*%as.matrix(data))
  sum.df<-cbind(sum.tmp,data.frame(ny=sum.tmp$n-sum.tmp$x))
  
  
  weighted.post.m<-c()
  post.m<-c()
  for(i in 1:nrow(omega.j)){
    prod.beta<-prod((beta.sub/beta(a,b))^(1-omega.j[i,]))
    post.m[i]<-(beta(a+sum.df$x[i],b+sum.df$ny[i])/beta(a,b))*prod.beta
  }
  
  
  omega.matrix<-list()
  weighted.m.matrix<-c()
  for(i in 1:2^num.elem){
    matrix.tmp<-diag(Nbask)
    matrix.tmp[upper.tri(matrix.tmp)]<-omega.all[i,]
    matrix.tmp[lower.tri(matrix.tmp)]<-t(matrix.tmp)[lower.tri(matrix.tmp)]
    omega.matrix[[i]]<-matrix.tmp
    
    eachk<-matrix(matrix.tmp[row(matrix.tmp)!=col(matrix.tmp)],ncol=Nbask-1,nrow=Nbask,byrow=T)
    ref.k<-bin2dec(eachk)
    
    post.m.matrix<-post.m[N.omega.j*(1:Nbask-1)+ref.k+1]
    
    weighted.m.matrix[i]<-prod(post.m.matrix)*prod(sapply(matrix.tmp[upper.tri(matrix.tmp)], function(x){x*mem.pri+(1-x)*(1-mem.pri)}))
    
    if(i==1){
      g.matrix<-matrix(ref.k,nrow=1)
    }else g.matrix<-rbind(g.matrix,ref.k)
  }
  
  sum.wmm<-sum(weighted.m.matrix)
  post.omega<-weighted.m.matrix/sum.wmm
  
  
  omega.similarity<-matrix(0,nrow=Nbask,ncol=Nbask)
  for(i in 1:2^num.elem){
    omega.similarity<-omega.similarity+omega.matrix[[i]]*post.omega[i]
  }
  
  list(Nomega.j=N.omega.j,sum.df=sum.df,post.omega=post.omega,similarity=omega.similarity,
       g.matrix=g.matrix,maximizer=omega.matrix[[which.max(weighted.m.matrix)]])
  
}

pi_name<-function(Nbask){
  for(i in 1:Nbask){
    if(i==1){
      name<-"pi[1]"
    }else name<-cbind(name,paste(c("pi[",i,"]"),collapse=""))
  }
  name
}

pi_success<-function(fit,pn)colMeans(as.data.frame(fit>pn))

read_result<-function(prt.dir,dir,title="noname",subtitle="output",sep){
  OC.tmp<-c("./../",prt.dir,"/",dir,"/",title,"_",subtitle,".csv")
  OC.name<-paste(OC.tmp,collapse="")
  OC<-read.table(OC.name,sep=sep)
  OC
}

ref_matrix<-function(Nbask){
  ref<-matrix(0,nrow=Nbask,ncol=Nbask)
  ref[upper.tri(ref)]<-1:choose(Nbask,2)
  ref<-ref+t(ref)
  
  ref
}

result_atleast1<-function(ps,Nbask,sim,q,true.clus){
  over.q<-ps>q
  clus.tab<-matrix(true.clus,ncol=Nbask,nrow=sim,byrow=T)
  match<-over.q==clus.tab
  
  if(sum(true.clus)==Nbask){
    fwer<-0
    poc<-mean(rowSums(match)==Nbask)
  }else if(sum(true.clus)==0){
    fwer<-mean(rowSums(match)<Nbask)
    poc<-0
  }else{
    fwer<-mean(rowSums(matrix(match[,!true.clus],ncol=sum(!true.clus)))<sum(!true.clus))
    poc<-mean(rowSums(matrix(match[,true.clus],ncol=sum(true.clus)))==sum(true.clus))
  }
  data.frame(fwer=fwer,all=poc)
}

result_bias<-function(pmean,true.pi,sim){
  mse<-colMeans((pmean-matrix(true.pi,ncol=length(true.pi),nrow=sim,byrow=T))^2,na.rm=T)
  mae<-colMeans(abs(pmean-matrix(true.pi,ncol=length(true.pi),nrow=sim,byrow=T)),na.rm=T)
  bias<-colMeans(pmean-matrix(true.pi,ncol=length(true.pi),nrow=sim,byrow=T),na.rm=T)
  list(mse=mse,mae=mae,bias=bias)
}

result_terminate<-function(ps,Nbask,sim,q.int){
  terminate<-c()
  for (i in 1:Nbask) {
    terminate[i]<-sum(as.numeric(sort(ps[,i])<=q.int))/sim
  }
  terminate
}

result_perfect<-function(ps,Nbask,sim,q,true.clus){
  over.q<-ps>q
  clus.tab<-matrix(true.clus,ncol=Nbask,nrow=sim,byrow=T)
  match<-over.q==clus.tab
  mean(rowSums(match)==Nbask)
}

result_power<-function(ps,Nbask,sim,q){
  power<-c()
  for (i in 1:Nbask) {
    power[i]<-sum(as.numeric(sort(ps[,i])>q))/sim
  }
  power
}

result_quantile<-function(ps,TIE){
  ps.mtx<-as.matrix(ps)
  q<-quantile(ps.mtx,1-TIE,na.rm=T)
  q
}

RoBoT_MCMC = function(n, r, pi0, niter, burnin = 5000, thin = 5, 
  return_MCMC_spls = FALSE, pi_prior_shape01 = 0.01, pi_prior_shape02 = 0.01,
  pi_prior_shape11 = 0.01, pi_prior_shape12 = 0.01, alpha = 2,
  mu_prior_mean = 0, mu_prior_sd = 2, tau_prior_loc = 0, tau_prior_scale = 1){
  
  #################################################################
  # INPUT: 
  # n: J length vector, total number of patients in each group
  # r: J length vector, number of responses in each group
  # niter: number of posterior samples we want 
  #################################################################

  J = length(n)
  
  if(length(pi0) != 1 & length(pi0) != J){
    stop("Length of pi0 must either be 1 or match with the length of n and r.")
  }

  if(length(pi0) == 1){
    pi0 = rep(pi0, J)
  }

  
  #################################################################
  # Set hyperparameters
  #################################################################
  
  # alpha = 2
  # pi_prior_shape1 = 0.01
  # pi_prior_shape2 = 0.01
  # mu_prior_mean = 0
  # mu_prior_sd = 2
  # tau_prior_loc = 0
  # tau_prior_scale = 1


  #################################################################
  # Initialize Markov chains
  #################################################################
  
  pi_spls = matrix(0, J, niter)
  theta_spls = matrix(0, J, niter)
  mu_spls = matrix(0, J, niter)
  tau_spls = matrix(0, J, niter)
  
  
  pi_init = r/n
  pi_init[pi_init < 0.01] = 0.01
  pi_init[pi_init > 0.99] = 0.99
  pi_spls[ , 1] = pi_init
  theta_spls[ , 1] = rep(0, J)
  mu_spls[ , 1] = rep(0, J)
  tau_spls[ , 1] = rep(1, J)

  #################################################################
  # MCMC burnin
  #################################################################
  
  output_C = .C("RoBoT_MCMC", pi = as.double(pi_spls[ , 1]),
                theta = as.double(theta_spls[ , 1]),
                mu = as.double(mu_spls[ , 1]),
                tau = as.double(tau_spls[ , 1]),
                alpha = as.double(alpha),
                J = as.integer(J), 
                n = as.double(n), 
                r = as.double(r), 
                pi0 = as.double(pi0),
                pi_prior_shape01 = as.double(pi_prior_shape01), 
                pi_prior_shape02 = as.double(pi_prior_shape02),
                pi_prior_shape11 = as.double(pi_prior_shape11), 
                pi_prior_shape12 = as.double(pi_prior_shape12),
                mu_prior_mean = as.double(mu_prior_mean),
                mu_prior_sd = as.double(mu_prior_sd),
                tau_prior_loc = as.double(tau_prior_loc),
                tau_prior_scale = as.double(tau_prior_scale),
                m = as.integer(5),
                niter = as.integer(burnin))
  
  pi_spls[ , 1] = output_C$pi
  theta_spls[ , 1] = output_C$theta
  mu_spls[ , 1] = output_C$mu
  tau_spls[ , 1] = output_C$tau
  
  #################################################################
  # MCMC iterations
  #################################################################

  for(i in 2:niter){
    
    pi_spls[ , i] = pi_spls[ , i-1]
    theta_spls[ , i] = theta_spls[ , i-1]
    mu_spls[ , i] = mu_spls[ , i-1]
    tau_spls[ , i] = tau_spls[ , i-1]

    output_C = .C("RoBoT_MCMC", pi = as.double(pi_spls[ , i]),
                  theta = as.double(theta_spls[ , i]),
                  mu = as.double(mu_spls[ , i]),
                  tau = as.double(tau_spls[ , i]),
                  alpha = as.double(alpha),
                  J = as.integer(J), 
                  n = as.double(n), 
                  r = as.double(r), 
                  pi0 = as.double(pi0),
                  pi_prior_shape01 = as.double(pi_prior_shape01), 
                  pi_prior_shape02 = as.double(pi_prior_shape02),
                  pi_prior_shape11 = as.double(pi_prior_shape11), 
                  pi_prior_shape12 = as.double(pi_prior_shape12),
                  mu_prior_mean = as.double(mu_prior_mean),
                  mu_prior_sd = as.double(mu_prior_sd),
                  tau_prior_loc = as.double(tau_prior_loc),
                  tau_prior_scale = as.double(tau_prior_scale),
                  m = as.integer(5),
                  niter = as.integer(thin))

    pi_spls[ , i] = output_C$pi
    theta_spls[ , i] = output_C$theta
    mu_spls[ , i] = output_C$mu
    tau_spls[ , i] = output_C$tau
  
  } 
  
  if(return_MCMC_spls){
    
    MCMC_spls = list()
    MCMC_spls$pi_spls = pi_spls
    MCMC_spls$theta_spls = theta_spls
    MCMC_spls$mu_spls = mu_spls
    MCMC_spls$tau_spls = tau_spls
    return(MCMC_spls)

  } else {

    return(pi_spls)
    
  }
}

update.part <- function(x, n, prior_part, part, a0 = 1, b0 = 1){
  R <- length(x)
  K <- nrow(part)
  
  p <- foreach(k = 1:K,.combine = "c")%do%{
    grp <- unlist(part[k, ])
    S <- aggregate(x, by=list(grp), sum)$x
    N <- aggregate(n, by=list(grp), sum)$x
    #calculate marginal probs m(s_j)
    prod((beta(a0+S, b0+N-S)/beta(a0,b0)))*prior_part[k]
  }
  #marginal probs m(s_j)
  mp <- sum(p)
  ## calculate posterior prob of each grouping structure
  p/mp
}

write_custom<-function(dir.name,Atitl,OC.save,sub1=NULL,sub2="output"){
  for(i in 1:length(OC.save)){
    if(!is.null(sub1)){
      name<-paste(c(Atitl[i],"_",sub1),collapse="")
    }else name<-Atitl[i]
    
    write_result(prt.dir=dir.name$prt.dir,
               dir=dir.name$dir,
               title=name,
               data=OC.save[[i]],
               subtitle=sub2)
  }
}

write_result<-function(prt.dir,dir,title="noname",data,subtitle="output"){
  OC.tmp<-c("./../",prt.dir,"/",dir,"/",title,"_",subtitle,".csv")
  OC.name<-paste(OC.tmp,collapse="")
  write.table(data, OC.name,row.names = FALSE,col.names = FALSE,sep = ",")
}

write_set<-function(prt.dir,dir,title="noname",data,subtitle="output"){
  OC.tmp<-c("./../",prt.dir,"/",dir,"/",title,"_",subtitle,".RData")
  OC.name<-paste(OC.tmp,collapse="")
  save(data, file=OC.name)
}





