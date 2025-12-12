analysis<-function(data,Nbask,pi.null,pi.alt,design,mcmc.set){
  ps <- rep(NA, Nbask)
  pu <- rep(NA, Nbask)
  mn <- rep(NA, Nbask)
  sd <- rep(NA, Nbask)
  cd <- rep(NA, Nbask)
  lc <- rep(NA, Nbask)
  md <- rep(NA, Nbask)
  uc <- rep(NA, Nbask)
  non0 <- (data$n != 0)

  if(strsplit(design,"_")[[1]][1]=="EXNEX"){
    start_time <- Sys.time()

    tau<-as.numeric(strsplit(strsplit(design,"_t")[[1]][2],"_")[[1]][1])
    pnex<-as.numeric(strsplit(strsplit(design,"_pn")[[1]][2],"_")[[1]][1])
    pinex<-as.numeric(strsplit(strsplit(design,"_pin")[[1]][2],"_")[[1]][1])
    
    data.mcmc <- list(
      "Nexch"=2, "Nmix"=3,
      "pMix" = c(rep((1-pnex)/2,2),pnex),
      # alternative weights for EX and EXNEX-2
      "Nstrata"=sum(non0),
      "n" = data$n[which(non0)],
      # original data
      "r" = data$x[which(non0)],
      # prior means and precisions for EX parameter mu
      "mu.mean"=c(logit(pi.null),logit(pi.alt)), "mu.prec"=c(1/(1/pi.null+1/(1-pi.null)-1),1/(1/pi.alt+1/(1-pi.alt)-1)),#correct(2022/05/14)
      # scale parameter of Half-Normal prior for tau
      "tau.HN.scale"=c(tau,tau),
      # NEX priors; make them strata-specific if needed
      "nex.mean"=logit(pinex), "nex.prec"=1/(1/pinex+1/(1-pinex)),
      "p.cut" = pi.null
    )
    
    var.name<-c("tau","p","mu","exch")
    
    fit <- jags(model.file = "exnex.txt",data = data.mcmc,parameters.to.save = var.name,
                n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)
    
    ps[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=pi.null)
    pu[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=(pi.null+pi.alt)/2)
    
    sum.fit<-fit[["BUGSoutput"]][["summary"]]
    mn[c(which(non0), seq(Nbask+1, Nbask+4))]<-sum.fit[c((1:sum(non0)+2),1,2,(1:2+sum(non0)+2))+sum(non0)*3,"mean"]
    sd[c(which(non0), seq(Nbask+1, Nbask+4))]<-sum.fit[c((1:sum(non0)+2),1,2,(1:2+sum(non0)+2))+sum(non0)*3,"sd"]
    cd[c(which(non0), seq(Nbask+1, Nbask+4))]<-sum.fit[c((1:sum(non0)+2),1,2,(1:2+sum(non0)+2))+sum(non0)*3,"Rhat"]
    lc[c(which(non0), seq(Nbask+1, Nbask+4))]<-sum.fit[c((1:sum(non0)+2),1,2,(1:2+sum(non0)+2))+sum(non0)*3,"2.5%"]
    md[c(which(non0), seq(Nbask+1, Nbask+4))]<-sum.fit[c((1:sum(non0)+2),1,2,(1:2+sum(non0)+2))+sum(non0)*3,"50%"]
    uc[c(which(non0), seq(Nbask+1, Nbask+4))]<-sum.fit[c((1:sum(non0)+2),1,2,(1:2+sum(non0)+2))+sum(non0)*3,"97.5%"]

    names(mn)<-c(pi_name(Nbask),"mu[1]","mu[2]","tau[1]","tau[2]")

    R<-as.numeric(sum.fit[1:(Nbask*3),"mean"])
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="CBHM"){
    start_time <- Sys.time()

    tau1<-as.numeric(strsplit(strsplit(design,"_tau1")[[1]][2],"_")[[1]][1])
    tau2<-as.numeric(strsplit(strsplit(design,"_tau2")[[1]][2],"_")[[1]][1])
    
    var.name<-c("theta0","phi","tausq","tausq2","tausq3","p","R")
    
    if(Nbask>1)distance<-H_distance(K=sum(non0),nik=data$n[which(non0)],rik=data$x[which(non0)],q1=pi.alt,q0=pi.null)
    else distance<-1
    data.mcmc <- list("n"=data$n[which(non0)], "Y"=data$x[which(non0)], "D"=distance, "K"=sum(non0), "zero"=rep(0,sum(non0)),
                      'mu0'=logit((pi.alt+pi.null)/2), tau1=tau1, tau2=tau2)
    fit <- jags(model.file = "cbhm_h.txt",data = data.mcmc,parameters.to.save = var.name,
                n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)
    
    ps[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=pi.null)
    pu[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=(pi.null+pi.alt)/2)
    
    sum.fit<-fit[["BUGSoutput"]][["summary"]]
    mn[c(which(non0), seq(Nbask+1, Nbask+5))]<-sum.fit[-(1:sum(non0)^2),"mean"]
    sd[c(which(non0), seq(Nbask+1, Nbask+5))]<-sum.fit[-(1:sum(non0)^2),"sd"]
    cd[c(which(non0), seq(Nbask+1, Nbask+5))]<-sum.fit[-(1:sum(non0)^2),"Rhat"]
    lc[c(which(non0), seq(Nbask+1, Nbask+5))]<-sum.fit[-(1:sum(non0)^2),"2.5%"]
    md[c(which(non0), seq(Nbask+1, Nbask+5))]<-sum.fit[-(1:sum(non0)^2),"50%"]
    uc[c(which(non0), seq(Nbask+1, Nbask+5))]<-sum.fit[-(1:sum(non0)^2),"97.5%"]

    names(mn)<-c(pi_name(Nbask),"phi","tausq","tausq2","tausq3","theta0")

    R<-rep(NA, (Nbask^2-Nbask)/2)
    R[1:((sum(non0)^2-sum(non0))/2)]<-fit[["BUGSoutput"]][["mean"]][["R"]][upper.tri(distance)]
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="BHM"){
    start_time <- Sys.time()
    
    var.name<-c("theta0","tausq","p")
    
    data.mcmc <- list("n"=data$n[which(non0)], "Y"=data$x[which(non0)], "K"=sum(non0),"mu0"=logit((pi.alt+pi.null)/2),
                      "apri"=as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1]),
                      "bpri"=as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1]))
    fit <- jags(model.file = "bhm.txt",data = data.mcmc,parameters.to.save = var.name,
                n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)
    
    ps[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=pi.null)
    pu[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=(pi.null+pi.alt)/2)
    
    sum.fit<-fit[["BUGSoutput"]][["summary"]]
    mn[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[,"mean"]
    sd[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[,"sd"]
    cd[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[,"Rhat"]
    lc[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[,"2.5%"]
    md[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[,"50%"]
    uc[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[,"97.5%"]

    names(mn)<-c(pi_name(Nbask),"tausq","theta0")

    R<-matrix(rep(0,Nbask^2),nrow = Nbask)[upper.tri(matrix(rep(0,Nbask^2),nrow = Nbask))]
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="BBM-JS"){
    start_time <- Sys.time()
    
    epsilon<-as.numeric(strsplit(strsplit(design,"_e")[[1]][2],"_")[[1]][1])
    tau<-as.numeric(strsplit(strsplit(design,"_t")[[1]][2],"_")[[1]][1])
    a<-as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1])
    b<-as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1])
    
    distance<-JS_distance(K=sum(non0),nk=data$n[which(non0)],xk=data$x[which(non0)],a=a,b=b)
    similarity<-(1-distance)^epsilon
    similarity[similarity<tau]<-0
    
    Aa<-similarity%*%(data$x[which(non0)]+a)
    Ab<-similarity%*%(data$n[which(non0)]-data$x[which(non0)]+b)
    
    ps[which(non0)]<-matrix(1-pbeta(pi.null, shape1 = Aa, shape2 = Ab),nrow=1)
    pu[which(non0)]<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = Aa, shape2 = Ab),nrow=1)
    mn<-rep(NA, Nbask*2)
    mn[c(which(non0), which(non0)+Nbask)]<-c(Aa/(Aa+Ab),Aa+Ab-data$n[which(non0)])
    sd[which(non0)]<-c(sqrt(Aa*Ab/((Aa+Ab)^2*(Aa+Ab+1))))
    cd[which(non0)]<-rep(1,sum(non0))
    lc[which(non0)]<-c(qbeta(0.025,shape1 = Aa, shape2 = Ab))
    md[which(non0)]<-c(qbeta(0.50,shape1 = Aa, shape2 = Ab))
    uc[which(non0)]<-c(qbeta(0.975,shape1 = Aa, shape2 = Ab))
    
    names(mn)<-c(pi_name(Nbask),sapply(1:Nbask,function(x){paste(c("pri_ESS[",x,"]"),collapse="")}))
    names(sd)<-pi_name(Nbask)
    names(cd)<-pi_name(Nbask)
    names(lc)<-pi_name(Nbask)
    names(md)<-pi_name(Nbask)
    names(uc)<-pi_name(Nbask)
    
    R<-rep(NA, (Nbask^2-Nbask)/2)
    R[1:((sum(non0)^2-sum(non0))/2)]<-similarity[upper.tri(similarity)]
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="MEM"){
    start_time <- Sys.time()
    
    a<-as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1])
    b<-as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1])
    pri<-as.numeric(strsplit(strsplit(design,"_p")[[1]][2],"_")[[1]][1])
    
    omega.list<-MEM_distance(data[data$n>0, ],sum(non0),a,b,mem.pri=pri)
    
    k<-1
    poj.R<-c()
    for(i in 1:sum(non0)){
      for(j in 1:omega.list$Nomega.j){
        poj.tmp<-sum(omega.list$post.omega[omega.list$g.matrix[,i]==(j-1)])
        poj.R[k]<-poj.tmp
        sumdf.tmp<-omega.list$sum.df[(i-1)*2^(sum(non0)-1)+j,]
        wpp.tmp<-poj.tmp*(1-pbeta(q=pi.null,shape1=a+sumdf.tmp$x,shape2=b+sumdf.tmp$ny))
        wpp.tmp.pu<-poj.tmp*(1-pbeta(q=(pi.null+pi.alt)/2,shape1=a+sumdf.tmp$x,shape2=b+sumdf.tmp$ny))
        wp.tmp<-poj.tmp*((a+sumdf.tmp$x)/(a+sumdf.tmp$x+b+sumdf.tmp$ny))
        wlci.tmp<-poj.tmp*qbeta(p=0.025,shape1=a+sumdf.tmp$x,shape2=b+sumdf.tmp$ny)
        wuci.tmp<-poj.tmp*qbeta(p=0.975,shape1=a+sumdf.tmp$x,shape2=b+sumdf.tmp$ny)
        if(k==1){
          post.omega.j<-c(i,j,poj.tmp)
          weighted.post.Pi<-c(i,j,wpp.tmp)
          weighted.post.pu<-c(i,j,wpp.tmp.pu)
          weighted.Pi<-c(i,j,wp.tmp)
          weighted.lowerCI<-c(i,j,wlci.tmp)
          weighted.upperCI<-c(i,j,wuci.tmp)
        }else{
          post.omega.j<-rbind(post.omega.j,c(i,j,poj.tmp))
          weighted.post.Pi<-rbind(weighted.post.Pi,c(i,j,wpp.tmp))
          weighted.post.pu<-rbind(weighted.post.pu,c(i,j,wpp.tmp.pu))
          weighted.Pi<-rbind(weighted.Pi,c(i,j,wp.tmp))
          weighted.lowerCI<-rbind(weighted.lowerCI,c(i,j,wlci.tmp))
          weighted.upperCI<-rbind(weighted.upperCI,c(i,j,wuci.tmp))
        }
        k<-k+1
      }
    }
    
    for(i in 1:sum(non0)){
      res.tmp<-c(sum(weighted.post.Pi[weighted.post.Pi[,1]==i,3]),
                 sum(weighted.post.pu[weighted.post.pu[,1]==i,3]),
                    sum(weighted.Pi[weighted.Pi[,1]==i,3]),
                    sum(weighted.lowerCI[weighted.lowerCI[,1]==i,3]),
                    sum(weighted.upperCI[weighted.upperCI[,1]==i,3]))
      if(i==1){
        res.estimate<-t(data.frame(res.tmp))
        colnames(res.estimate)<-c("weighted_post_Pi",
                                  "weighted_post_pu",
                                     "weighted_Pi",
                                     "weighted_lowerCI",
                                     "weighted_upperCI")
      }else res.estimate<-rbind(res.estimate,res.tmp)
    }
    
    # rname<-pi_name(Nbask)
    # row.names(res.estimate)<-rname
    
    ps[which(non0)]<-res.estimate[,"weighted_post_Pi"]
    pu[which(non0)]<-res.estimate[,"weighted_post_pu"]

    mn<-rep(NA, Nbask*2)
    mn[c(which(non0), which(non0)+Nbask)]<-c(res.estimate[,"weighted_Pi"],omega.list$similarity%*%data$n[which(non0)]+a+b-data$n[which(non0)])
    sd[which(non0)]<-rep(0,sum(non0))
    cd[which(non0)]<-rep(0,sum(non0))
    lc[which(non0)]<-res.estimate[,"weighted_lowerCI"]
    md[which(non0)]<-rep(0,sum(non0))
    uc[which(non0)]<-res.estimate[,"weighted_upperCI"]

    names(mn)<-c(pi_name(Nbask),sapply(1:Nbask,function(x){paste(c("pri_ESS[",x,"]"),collapse="")}))
    
    R<-rep(NA, Nbask^2-Nbask)
    R[c(1:((sum(non0)^2-sum(non0))/2), ((Nbask^2-Nbask)/2)+1:((sum(non0)^2-sum(non0))/2))]<-c(omega.list$similarity[upper.tri(omega.list$similarity)],omega.list$maximizer[upper.tri(omega.list$maximizer)])
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="independent"){
    start_time <- Sys.time()
    
    a<-as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1])
    b<-as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1])

    Aa<-rep(NA, Nbask)
    Ab<-rep(NA, Nbask)
    Aa[which(non0)]<-data$x[which(non0)]+a
    Ab[which(non0)]<-data$n[which(non0)]-data$x[which(non0)]+b
    
    ps<-matrix(1-pbeta(pi.null, shape1 = Aa, shape2 = Ab),nrow=1)
    pu<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = Aa, shape2 = Ab),nrow=1)
    mn<-c(Aa/(Aa+Ab))
    sd<-c(sqrt(Aa*Ab/((Aa+Ab)^2*(Aa+Ab+1))))
    cd[which(non0)]<-rep(1,sum(non0))
    lc<-c(qbeta(0.025,shape1 = Aa, shape2 = Ab))
    md<-c(qbeta(0.500,shape1 = Aa, shape2 = Ab))
    uc<-c(qbeta(0.975,shape1 = Aa, shape2 = Ab))
    
    names(mn)<-pi_name(Nbask)
    names(sd)<-pi_name(Nbask)
    names(cd)<-pi_name(Nbask)
    names(lc)<-pi_name(Nbask)
    names(md)<-pi_name(Nbask)
    names(uc)<-pi_name(Nbask)
    
    R<-matrix(rep(0,Nbask^2),nrow = Nbask)[upper.tri(matrix(rep(0,Nbask^2),nrow = Nbask))]
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="UPSiDe-D"){
    start_time <- Sys.time()
    
    Ncomb<-sum(non0)*(sum(non0)-1)/2
    z<-as.numeric(strsplit(strsplit(design,"_z")[[1]][2],"_")[[1]][1])
    ref<-ref_matrix(sum(non0))
    zvct<-rep(z,Ncomb)
    
    
    imtx<-matrix(0,ncol=sum(non0)-1,nrow=sum(non0))
    for(i in 1:sum(non0)){
      imtx[i,]<-(1:sum(non0))[1:sum(non0)!=i]
    }
    
    Means <- data$x[which(non0)]/data$n[which(non0)]
    Means[Means<0.001] <- 0.001
    Means[Means>0.999] <- 0.999
    upM <-as.numeric(strsplit(strsplit(design,"_M")[[1]][2],"_")[[1]][1])
    
    jagsData <- list(x = data$x[which(non0)], n = data$n[which(non0)], I = sum(non0), imtx = imtx, zvct = zvct,
                     ref = ref, upM = upM, Means = Means)
    jagsParam <- c("u", "pi", "wmtx", "Mwmtx")
    model.name<-paste(c("UPSiDe-D",strsplit(design,"_")[[1]][2],".txt"),collapse="")
    fit <- jags(model.file = model.name,data = jagsData, parameters.to.save = jagsParam,
                n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)
    
    ps[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$pi,pn=pi.null)
    pu[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$pi,pn=(pi.null+pi.alt)/2)
    
    sum.fit<-fit[["BUGSoutput"]][["summary"]]
    mn[c(which(non0), Nbask+1)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+1),"mean"]
    sd[c(which(non0), Nbask+1)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+1),"sd"]
    cd[c(which(non0), Nbask+1)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+1),"Rhat"]
    lc[c(which(non0), Nbask+1)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+1),"2.5%"]
    md[c(which(non0), Nbask+1)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+1),"50%"]
    uc[c(which(non0), Nbask+1)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+1),"97.5%"]

    names(mn)<-c(pi_name(Nbask),"u")
    
    R<-rep(NA, length(c(1:(Nbask*(Nbask-1)),Nbask*(Nbask-1)+Nbask+1+1:(Nbask*(Nbask-1)))))
    R[c(1:(sum(non0)*(sum(non0)-1)), (Nbask*(Nbask-1))+1:(sum(non0)*(sum(non0)-1)))]<-sum.fit[c(1:(sum(non0)*(sum(non0)-1)),sum(non0)*(sum(non0)-1)+sum(non0)+1+1:(sum(non0)*(sum(non0)-1))),"mean"]
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
    
  }
  else if(strsplit(design,"_")[[1]][1]=="UPSiDe-FIX"){
    start_time <- Sys.time()
    
    Ncomb<-sum(non0)*(sum(non0)-1)/2
    distance<-KL_distance2(K=sum(non0),nik=data$n[which(non0)],rik=data$x[which(non0)],q1=pi.alt,q0=pi.null,a=1,b=1)
    s<-as.numeric(strsplit(strsplit(design,"_s")[[1]][2],"_")[[1]][1])
    dsmtx<-exp(-distance/s)
    diag(dsmtx)<-0
    
    ref<-ref_matrix(sum(non0))
    imtx<-matrix(0,ncol=sum(non0)-1,nrow=sum(non0))
    for(i in 1:sum(non0)){
      imtx[i,]<-(1:sum(non0))[1:sum(non0)!=i]
    }
    
    sumw<-sum(dsmtx[upper.tri(dsmtx)])
    w<-dsmtx[upper.tri(dsmtx)]/sumw
    
    Means <- data$x[which(non0)]/data$n[which(non0)]
    Means[Means<0.001] <- 0.001
    Means[Means>0.999] <- 0.999
    upM <-as.numeric(strsplit(strsplit(design,"_M")[[1]][2],"_")[[1]][1])
    
    avg<-c()
    prec<-c()
    alp<-c()
    bet<-c()
    wmtx<-matrix(0,nrow=sum(non0),ncol=sum(non0)-1)
    Mwmtx<-matrix(0,nrow=sum(non0),ncol=sum(non0)-1)
    MwUImtx<-matrix(0,nrow=sum(non0),ncol=sum(non0)-1)
    for(i in 1:sum(non0)){
      for(j in 1:(sum(non0)-1)){
        wmtx[i,j] <- w[ref[i,imtx[i,j]]] / 2
        Mwmtx[i,j] <- upM * wmtx[i,j]
        MwUImtx[i,j] <- Mwmtx[i,j] * min(1/(0.05*0.95), 1 / (Means[imtx[i,j]] * (1 - Means[imtx[i,j]])))
      }
      avg[i]<-sum(wmtx[i,1:(sum(non0)-1)] * Means[imtx[i,1:(sum(non0)-1)]]) / sum(wmtx[i,1:(sum(non0)-1)])
      prec[i]<-sum(MwUImtx[i,1:(sum(non0)-1)])
      alp[i] <- max(avg[i] * (avg[i] * (1 - avg[i]) * prec[i] - 1),0.5)
      bet[i] <- max((1 - avg[i]) * (avg[i] * (1 - avg[i]) * prec[i] - 1),0.5)
    }
    Aa<-rep(NA, Nbask)
    Ab<-rep(NA, Nbask)
    Aa[which(non0)]<-alp+data$x[which(non0)]
    Ab[which(non0)]<-bet+(data$n[which(non0)]-data$x[which(non0)])
    
    ps<-matrix(1-pbeta(pi.null, shape1 = Aa, shape2 = Ab),nrow=1)
    pu<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = Aa, shape2 = Ab),nrow=1)
    mn<-c(Aa/(Aa+Ab),Aa+Ab-data$n)
    sd<-c(sqrt(Aa*Ab/((Aa+Ab)^2*(Aa+Ab+1))))
    cd[which(non0)]<-rep(1,sum(non0))
    lc<-c(qbeta(0.025,shape1 = Aa, shape2 = Ab))
    md<-c(qbeta(0.500,shape1 = Aa, shape2 = Ab))
    uc<-c(qbeta(0.975,shape1 = Aa, shape2 = Ab))
    
    names(mn)<-c(pi_name(Nbask),sapply(1:Nbask,function(x){paste(c("pri_ESS[",x,"]"),collapse="")}))
    names(sd)<-pi_name(Nbask)
    names(cd)<-pi_name(Nbask)
    names(lc)<-pi_name(Nbask)
    names(md)<-pi_name(Nbask)
    names(uc)<-pi_name(Nbask)
    
    R<-rep(NA, length(c(1:(Nbask*(Nbask-1)),1:(Nbask*(Nbask-1)))))
    R[c(1:(sum(non0)*(sum(non0)-1)), (Nbask*(Nbask-1))+1:(sum(non0)*(sum(non0)-1)))]<-c(Mwmtx,wmtx)
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
    
  }
  else if(strsplit(design,"_")[[1]][1]=="UPSiDe-JS"){
    start_time <- Sys.time()
    
    Ncomb<-sum(non0)*(sum(non0)-1)/2
    distance<-KL_distance2(K=sum(non0),nik=data$n[which(non0)],rik=data$x[which(non0)],q1=pi.alt,q0=pi.null,a=1,b=1)
    dvct<-c(distance[upper.tri(distance)])
    r<-as.numeric(strsplit(strsplit(design,"_r")[[1]][2],"_")[[1]][1])
    
    ref<-ref_matrix(sum(non0))
    imtx<-matrix(0,ncol=sum(non0)-1,nrow=sum(non0))
    for(i in 1:sum(non0)){
      imtx[i,]<-(1:sum(non0))[1:sum(non0)!=i]
    }
    
    Means <- data$x[which(non0)]/data$n[which(non0)]
    Means[Means<0.001] <- 0.001
    Means[Means>0.999] <- 0.999
    upM <-as.numeric(strsplit(strsplit(design,"_M")[[1]][2],"_")[[1]][1])
    
    jagsData <- list(x = data$x[which(non0)], n = data$n[which(non0)], I = sum(non0), imtx = imtx, dvct = dvct, Ncomb = Ncomb,
                     spri = r, ref = ref, upM = upM, Means = Means)
    jagsParam <- c("u", "pi", "s", "wmtx", "Mwmtx")
    model.name<-paste(c("UPSiDe-JS",strsplit(design,"_")[[1]][2],".txt"),collapse="")
    fit <- jags(model.file = model.name,data = jagsData, parameters.to.save = jagsParam,
                n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)
    
    ps[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$pi,pn=pi.null)
    pu[which(non0)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$pi,pn=(pi.null+pi.alt)/2)
    
    sum.fit<-fit[["BUGSoutput"]][["summary"]]
    mn[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+2),"mean"]
    sd[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+2),"sd"]
    cd[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+2),"Rhat"]
    lc[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+2),"2.5%"]
    md[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+2),"50%"]
    uc[c(which(non0), Nbask+1, Nbask+2)]<-sum.fit[sum(non0)*(sum(non0)-1)+1:(sum(non0)+2),"97.5%"]
    
    names(mn)<-c(pi_name(Nbask),"s","u")

    R<-rep(NA, length(c(1:(Nbask*(Nbask-1)),Nbask*(Nbask-1)+Nbask+2+1:(Nbask*(Nbask-1)))))
    R[c(1:(sum(non0)*(sum(non0)-1)), (Nbask*(Nbask-1))+1:(sum(non0)*(sum(non0)-1)))]<-sum.fit[c(1:(sum(non0)*(sum(non0)-1)),sum(non0)*(sum(non0)-1)+sum(non0)+2+1:(sum(non0)*(sum(non0)-1))),"mean"]
    
    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
    
  }
  else if(strsplit(design,"_")[[1]][1]=="RoBoT"){
    start_time <- Sys.time()

    a0<-as.numeric(strsplit(strsplit(design,"_a0")[[1]][2],"_")[[1]][1])
    b0<-as.numeric(strsplit(strsplit(design,"_b0")[[1]][2],"_")[[1]][1])
    a1<-as.numeric(strsplit(strsplit(design,"_a1")[[1]][2],"_")[[1]][1])
    b1<-as.numeric(strsplit(strsplit(design,"_b1")[[1]][2],"_")[[1]][1])
    alpha<-as.numeric(strsplit(strsplit(design,"_alpha")[[1]][2],"_")[[1]][1])
    mu<-as.numeric(strsplit(strsplit(design,"_mu")[[1]][2],"_")[[1]][1])
    sigma<-as.numeric(strsplit(strsplit(design,"_sigma")[[1]][2],"_")[[1]][1])
    tau<-as.numeric(strsplit(strsplit(design,"_tau")[[1]][2],"_")[[1]][1])
    xi<-as.numeric(strsplit(strsplit(design,"_xi")[[1]][2],"_")[[1]][1])

    MCMC_spls = RoBoT_MCMC(
      data$n[which(non0)], data$x[which(non0)], pi.alt, niter = 10000, burnin = 50000, thin = 5,
      return_MCMC_spls = TRUE, pi_prior_shape01 = a0, pi_prior_shape02 = b0,
      pi_prior_shape11 = a1, pi_prior_shape12 = b1, alpha = alpha,
      mu_prior_mean = mu, mu_prior_sd = sigma, tau_prior_loc = tau, tau_prior_scale = xi
    )

    pi_spls = MCMC_spls$pi_spls
    theta_spls = MCMC_spls$theta_spls
    mu_spls = MCMC_spls$mu_spls
    tau_spls = MCMC_spls$tau_spls

    ps[which(non0)]<-rowMeans(pi_spls > pi.null)
    pu[which(non0)]<-rowMeans(pi_spls > (pi.null+pi.alt)/2)
    
    mn[c(which(non0))]<-apply(pi_spls, 1, mean)
    sd[c(which(non0))]<-apply(pi_spls, 1, stats::sd)
    cd[c(which(non0))]<-rep(1, sum(non0))
    lc[c(which(non0))]<-apply(pi_spls, 1, function(x){quantile(x, probs = 0.025)})
    md[c(which(non0))]<-apply(pi_spls, 1, median)
    uc[c(which(non0))]<-apply(pi_spls, 1, function(x){quantile(x, probs = 0.975)})
    
    names(mn)<-pi_name(Nbask)

    comb <- combn(sum(non0), 2)
    R <- rep(NA, dim(combn(Nbask, 2))[2])
    for (i in 1:dim(comb)[2]) {
      R[i] <- cor(pi_spls[comb[1, i], ], pi_spls[comb[2, i], ])
    }

    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="Cluster-BHM"){
    start_time <- Sys.time()
    
    psi<-as.numeric(strsplit(strsplit(design,"_psi")[[1]][2],"_")[[1]][1])
    a<-as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1])
    b<-as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1])
    
    Aa<-rep(NA, Nbask)
    Ab<-rep(NA, Nbask)
    Aa[which(non0)]<-data$x[which(non0)]+a
    Ab[which(non0)]<-data$n[which(non0)]-data$x[which(non0)]+b
    
    pu<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = Aa, shape2 = Ab),nrow=1)
    cond <- (pu > psi)

    ps <- rep(NA, Nbask)
    pu <- rep(NA, Nbask)
    mn <- rep(NA, Nbask+4)
    sd <- rep(NA, Nbask+4)
    cd <- rep(NA, Nbask+4)
    lc <- rep(NA, Nbask+4)
    md <- rep(NA, Nbask+4)
    uc <- rep(NA, Nbask+4)

    names(mn)<-cbind(pi_name(Nbask), "tausq[1]", "tausq[2]", "theta0[1]", "theta0[2]")
    names(sd)<-cbind(pi_name(Nbask), "tausq[1]", "tausq[2]", "theta0[1]", "theta0[2]")
    names(cd)<-cbind(pi_name(Nbask), "tausq[1]", "tausq[2]", "theta0[1]", "theta0[2]")
    names(lc)<-cbind(pi_name(Nbask), "tausq[1]", "tausq[2]", "theta0[1]", "theta0[2]")
    names(md)<-cbind(pi_name(Nbask), "tausq[1]", "tausq[2]", "theta0[1]", "theta0[2]")
    names(uc)<-cbind(pi_name(Nbask), "tausq[1]", "tausq[2]", "theta0[1]", "theta0[2]")

    # valid case -------------
    if (sum(cond, na.rm=T) > 1) {
      var.name<-c("theta0","tausq","p")

      data.mcmc <- list("n"=data$n[which(cond)], "Y"=data$x[which(cond)], "K"=sum(cond, na.rm=T),"mu0"=logit(pi.alt),
                        "apri"=as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1]),
                        "bpri"=as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1]))
      fit <- jags(model.file = "bhm.txt",data = data.mcmc,parameters.to.save = var.name,
                  n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                  n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)

      ps[which(cond)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=pi.null)
      pu[which(cond)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=(pi.null+pi.alt)/2)

      sum.fit<-fit[["BUGSoutput"]][["summary"]]
      mn[c(which(cond), Nbask+1, Nbask+3)]<-sum.fit[,"mean"]
      sd[c(which(cond), Nbask+1, Nbask+3)]<-sum.fit[,"sd"]
      cd[c(which(cond), Nbask+1, Nbask+3)]<-sum.fit[,"Rhat"]
      lc[c(which(cond), Nbask+1, Nbask+3)]<-sum.fit[,"2.5%"]
      md[c(which(cond), Nbask+1, Nbask+3)]<-sum.fit[,"50%"]
      uc[c(which(cond), Nbask+1, Nbask+3)]<-sum.fit[,"97.5%"]

    } else if (sum(cond, na.rm=T) == 1) {
      # only one member
      ps[which(cond)]<-matrix(1-pbeta(pi.null, shape1 = Aa, shape2 = Ab),nrow=1)[, which(cond)]
      pu[which(cond)]<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = Aa, shape2 = Ab),nrow=1)[, which(cond)]
      mn[which(cond)]<-c(Aa/(Aa+Ab))[which(cond)]
      sd[which(cond)]<-c(sqrt(Aa*Ab/((Aa+Ab)^2*(Aa+Ab+1))))[which(cond)]
      cd[which(cond)]<-rep(1,Nbask)[which(cond)]
      lc[which(cond)]<-c(qbeta(0.025,shape1 = Aa, shape2 = Ab))[which(cond)]
      md[which(cond)]<-c(qbeta(0.500,shape1 = Aa, shape2 = Ab))[which(cond)]
      uc[which(cond)]<-c(qbeta(0.975,shape1 = Aa, shape2 = Ab))[which(cond)]
    }

    # invalid case -------------
    if (sum(!cond, na.rm=T) > 1) {
      var.name<-c("theta0","tausq","p")

      data.mcmc <- list("n"=data$n[which(!cond)], "Y"=data$x[which(!cond)], "K"=sum(!cond, na.rm=T),"mu0"=logit(pi.null),
                        "apri"=as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1]),
                        "bpri"=as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1]))
      fit <- jags(model.file = "bhm.txt",data = data.mcmc,parameters.to.save = var.name,
                  n.chains = mcmc.set$Nchain,n.iter=mcmc.set$Niter,
                  n.thin=mcmc.set$Nthin,n.burnin = mcmc.set$Nburn,quiet=T,DIC = F)

      ps[which(!cond)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=pi.null)
      pu[which(!cond)]<-pi_success(fit=fit[["BUGSoutput"]][["sims.list"]]$p,pn=(pi.null+pi.alt)/2)

      sum.fit<-fit[["BUGSoutput"]][["summary"]]
      mn[c(which(!cond), Nbask+2, Nbask+4)]<-sum.fit[,"mean"]
      sd[c(which(!cond), Nbask+2, Nbask+4)]<-sum.fit[,"sd"]
      cd[c(which(!cond), Nbask+2, Nbask+4)]<-sum.fit[,"Rhat"]
      lc[c(which(!cond), Nbask+2, Nbask+4)]<-sum.fit[,"2.5%"]
      md[c(which(!cond), Nbask+2, Nbask+4)]<-sum.fit[,"50%"]
      uc[c(which(!cond), Nbask+2, Nbask+4)]<-sum.fit[,"97.5%"]

    } else if (sum(!cond, na.rm=T) == 1) {
      # only one member
      ps[which(!cond)]<-matrix(1-pbeta(pi.null, shape1 = Aa, shape2 = Ab),nrow=1)[, which(!cond)]
      pu[which(!cond)]<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = Aa, shape2 = Ab),nrow=1)[, which(!cond)]
      mn[which(!cond)]<-c(Aa/(Aa+Ab))[which(!cond)]
      sd[which(!cond)]<-c(sqrt(Aa*Ab/((Aa+Ab)^2*(Aa+Ab+1))))[which(!cond)]
      cd[which(!cond)]<-rep(1,Nbask)[which(!cond)]
      lc[which(!cond)]<-c(qbeta(0.025,shape1 = Aa, shape2 = Ab))[which(!cond)]
      md[which(!cond)]<-c(qbeta(0.500,shape1 = Aa, shape2 = Ab))[which(!cond)]
      uc[which(!cond)]<-c(qbeta(0.975,shape1 = Aa, shape2 = Ab))[which(!cond)]

    }

    R<-as.numeric(cond)

    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="local-MEM"){
    start_time <- Sys.time()

    delta <- as.numeric(strsplit(strsplit(design,"_delta")[[1]][2],"_")[[1]][1])
    a0 <- as.numeric(strsplit(strsplit(design,"_a")[[1]][2],"_")[[1]][1])
    b0 <- as.numeric(strsplit(strsplit(design,"_b")[[1]][2],"_")[[1]][1])

    part <- get.part(R = sum(non0), max_cl = sum(non0))
    K <- nrow(part)
    n_bk <- apply(part, 1, function(x){
      length(unique(x))
    })

    prior_mat <- n_bk^delta / sum(n_bk^delta)
    prior_part <- prior_mat

    post_part <- update.part(
      x = data$x[which(non0)],
      n = data$n[which(non0)],
      prior_part = prior_part,
      part = part
    )

    # partition with the maximum pp
    part_hat <- unlist(part[which.max(post_part), ])
    p_hat <- post_part[which.max(post_part)]
    
    a1 <- b1 <- NULL
    for(b in 1:sum(non0)){
      # exchangeability of B baskets
      exc <- rep(0, sum(non0))
      exc[part_hat==part_hat[b]] <- p_hat
      exc[b] <- 1
      a1[b] <- a0 + sum(data$x[which(non0)]*exc)
      b1[b] <- b0 + sum(data$n[which(non0)]*exc)
    }
    
    ps[which(non0)]<-matrix(1-pbeta(pi.null, shape1 = a1, shape2 = b1),nrow=1)
    pu[which(non0)]<-matrix(1-pbeta((pi.null+pi.alt)/2, shape1 = a1, shape2 = b1),nrow=1)
    mn[which(non0)]<-c(a1/(a1+b1))
    sd[which(non0)]<-c(sqrt(a1*b1/((a1+b1)^2*(a1+b1+1))))
    cd[which(non0)]<-rep(1,sum(non0))
    lc[which(non0)]<-c(qbeta(0.025,shape1 = a1, shape2 = b1))
    md[which(non0)]<-c(qbeta(0.500,shape1 = a1, shape2 = b1))
    uc[which(non0)]<-c(qbeta(0.975,shape1 = a1, shape2 = b1))
    
    names(mn)<-pi_name(Nbask)
    names(sd)<-pi_name(Nbask)
    names(cd)<-pi_name(Nbask)
    names(lc)<-pi_name(Nbask)
    names(md)<-pi_name(Nbask)
    names(uc)<-pi_name(Nbask)
    
    R<-rep(NA, Nbask+1)
    R[c(1, which(non0)+1)]<-c(p_hat, part_hat)

    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
  else if(strsplit(design,"_")[[1]][1]=="Adaptive-Lasso"){
    start_time <- Sys.time()

    lambda<-as.numeric(strsplit(strsplit(design,"_lambda")[[1]][2],"_")[[1]][1])
    gamma<-as.numeric(strsplit(strsplit(design,"_gamma")[[1]][2],"_")[[1]][1])

    options(warn = 0)
    res.estimate<-basket_aLasso(data$x[which(non0)], data$n[which(non0)], lambda, gamma)
    options(warn = 2)

    pEst = res.estimate$pEst
    penFactor = res.estimate$penFactor
    cd.tmp = res.estimate$cd

    ps[which(non0)]<-rep(0, sum(non0))
    pu[which(non0)]<-rep(0, sum(non0))
    
    mn[c(which(non0))]<-pEst
    sd[c(which(non0))]<-rep(0, sum(non0))
    cd[c(which(non0))]<-cd.tmp
    lc[c(which(non0))]<-rep(0, sum(non0))
    md[c(which(non0))]<-rep(0, sum(non0))
    uc[c(which(non0))]<-rep(0, sum(non0))
    
    names(mn)<-pi_name(Nbask)

    R<-rep(NA, (Nbask^2-Nbask)/2)
    R[1:((sum(non0)^2-sum(non0))/2)]<-penFactor[upper.tri(penFactor)]

    diff_time <- Sys.time() - start_time
    list(ps=ps,pu=pu,mn=mn,sd=sd,cd=cd,R=R,lc=lc,md=md,uc=uc,et=diff_time)
  }
}
