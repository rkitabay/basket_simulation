trial<-function(sim,design,true.pi=true.pi){
  
  mcmc.set<-data.frame(Nchain=2,Niter=10000,Nthin=2,Nburn=1000)
  orres<-list()
  
  if((sim$Npat$method=="random"||sim$Npat$method=="realloc"||sim$Npat$method=="poisson")&&length(sim$Npat$max)>1) Nanal<-length(sim$Npat$max)
  else if(sim$Npat$method=="fix"&&nrow(sim$Npat$nfix)>1) Nanal<-nrow(sim$Npat$nfix)
  else Nanal<-1
  length(orres)<-Nanal
  
  #To applicate UPSiDe to interim analyses
  design.tmp<-design
  
  set.seed(sim$seed+sum(round(exp(true.pi)*10000)))
  if(sim$Npat$method=="realloc"){
    seed.tab<-round(runif(n=sim$Nsim*Nanal)*100000)
  }else{
    monitor.lump<-create_lump_data(Nbask=sim$Nbask,Npat=sim$Npat,true.pi,Nsim=sim$Nsim,Nanal)
  }
    
  for(s in 1:sim$Nsim){
    cntn<-rep(1,sim$Nbask)
    monitor<-data.frame(i=1:sim$Nbask,x=0,n=0)

    for(t in 1:Nanal){
      if(sim$Npat$method=="realloc"){
        monitor<-create_data(Nbask=sim$Nbask,Npat=sim$Npat,true.pi=true.pi,base=monitor,seed=seed.tab[t+Nanal*(s-1)],cntn=cntn,anal.time=t)
      }else{
        monitor$x<-monitor$x+cntn*monitor.lump$x[[t]][s,]
        monitor$n<-monitor$n+cntn*monitor.lump$n[[t]][s,]
      }

      #To applicate UPSiDe to interim analyses
       if(strsplit(design,"-")[[1]][1]=="UPSiDe"&&Nanal>1){
         design.div<-strsplit(design,"_")[[1]]
         design.div[3]<-paste(c("M", strsplit(strsplit(design,"_")[[1]][3],"[M-]")[[1]][t+1]), collapse = "")
         design.tmp<-paste(design.div,collapse="_")
         output<-analysis(data=monitor,Nbask=sim$Nbask,pi.null=sim$pi.null,
                          pi.alt=sim$pi.alt,design=design.tmp,mcmc.set)
       } else {
         output<-analysis(data=monitor,Nbask=sim$Nbask,pi.null=sim$pi.null,
                          pi.alt=sim$pi.alt,design,mcmc.set)
       }

      #----append----
      if(s==1){
        orres[[t]]<-list()
        length(orres[[t]])<-length(output)+2
        for(i in 1:length(output)){
          orres[[t]][[i]]<-output[[i]]
        }
        orres[[t]][[length(output)+1]]<-monitor$n
        orres[[t]][[length(output)+2]]<-monitor$x
        names(orres[[t]])<-c(names(output),"n","x")
      }else{
        for(i in 1:length(output)){
          orres[[t]][[i]]<-rbind(orres[[t]][[i]],output[[i]])
        }
        orres[[t]][[length(output)+1]]<-rbind(orres[[t]][[length(output)+1]],monitor$n)
        orres[[t]][[length(output)+2]]<-rbind(orres[[t]][[length(output)+2]],monitor$x)
      }
      if (strsplit(design,"_")[[1]][1]!="Adaptive-Lasso") {
      cntn<-c(output$pu>sim$q.int)
      } else {
        cntn<-c(output$mn>sim$q.int)
      }
      cntn[is.na(cntn)] <- FALSE
    }
    cat(s,":")
    if(!is.null(sim$q)&&s>1)cat(result_power(ps=as.data.frame(orres[[Nanal]][c(names(output)=="ps",F,F)]),
                                             sim$Nbask,sim=s,sim$q))
    cat("\n")
  }
  orres
}
