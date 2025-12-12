simulation_q<-function(sim,design){
  
  true.pi<-rep(sim$pi.null,sim$Nbask)
  #----calculate q----
  q.OC<-trial(sim,design,true.pi=true.pi)

  prt.dir<-sim$dir.name$prt.dir
  dir<-sim$dir.name$dir
  #prepare to output
  if (strsplit(design,"_")[[1]][1]!="Adaptive-Lasso") {
  q<-result_quantile(ps=q.OC[[length(q.OC)]]$ps,sim$TIE)
  } else {
    q<-result_quantile(ps=q.OC[[length(q.OC)]]$mn,sim$TIE)
  }
  
  write_result(prt.dir,dir,title="q",data=q)
  
  
  for(i in 1:length(q.OC)){
    #Storage the tables
    Atitl<-sapply(c("power.OC","p.upper","mean","sd","cd","R","lc","md","uc","exectime","n","x"),function(x){paste(c(x,i),collapse="")})
    OC.save<-list(q.OC[[i]]$ps,q.OC[[i]]$pu,q.OC[[i]]$mn,
                  q.OC[[i]]$sd,q.OC[[i]]$cd,q.OC[[i]]$R,
                  q.OC[[i]]$lc,q.OC[[i]]$md,q.OC[[i]]$uc,
                  q.OC[[i]]$et,q.OC[[i]]$n,q.OC[[i]]$x)#1-2,1-3,2-3,1-4,2-4,...
    write_custom(sim$dir.name,Atitl,OC.save,sub2="q")
    
    #summaries
    if (strsplit(design,"_")[[1]][1]!="Adaptive-Lasso") {
    power<-result_power(ps=q.OC[[i]]$ps,sim$Nbask,sim$Nsim,q)
    terminate<-result_terminate(ps=q.OC[[i]]$pu,sim$Nbask,sim$Nsim,sim$q.int)
    } else {
      power<-result_power(ps=q.OC[[i]]$mn,sim$Nbask,sim$Nsim,q)
      terminate<-result_terminate(ps=q.OC[[i]]$mn,sim$Nbask,sim$Nsim,sim$q.int)
    }
    estimate<-list(name=colnames(q.OC[[i]]$mn),mean=colMeans(q.OC[[i]]$mn, na.rm=T),sd=colMeans(q.OC[[i]]$sd, na.rm=T))
    bias<-result_bias(pmean=q.OC[[i]]$mn[,1:sim$Nbask],true.pi,sim$Nsim)
    mean_et<-mean(q.OC[[i]]$et)
    
    #Storage the summaries
    Atitl<-sapply(c("power","terminate","estimate","bias","mean_exectime"),function(x){paste(c(x,i),collapse="")})
    OC.save<-list(power,terminate,estimate,bias,mean_et)
    write_custom(sim$dir.name,Atitl,OC.save,sub2="q")
  }
  
  q
}
simulation_q_strg<-function(sim,design){
  
  sim["seed"]<-sim$seed+sim$Nsim*10
  #----calculate q----
  pi.strg<-c(sim$pi.null,rep(sim$pi.alt,sim$Nbask-1))
  q.OC<-trial(sim,design,true.pi=pi.strg)
  
  prt.dir<-sim$dir.name$prt.dir
  dir<-sim$dir.name$dir
  #Storage the result
  write_result(prt.dir,dir,title="q.OC",data=q.OC$ps)
  #prepare to output
  q<-result_quantile(ps=q.OC$ps[,1],TIE)
  write_result(prt.dir,dir,title="q",subtitle="strong",data=q)
  
  q
}
simulation_power<-function(sim,design){
  
  if (sim$scen <= sim$Nbask) {
    true.pi<-create_pi(sim$Nbask,sim$scen,sim$pi.null,sim$pi.alt)
  } else {
    csv.tmp <- c("./../data/true_pi", sim$Nbask, ".csv")
    csv.name <- paste(csv.tmp, collapse = "")
    csv.pi.all <- read.csv(csv.name, header = FALSE)
    csv.pi <- csv.pi.all[csv.pi.all[, 1] == sim$scennm, 2:(sim$Nbask+1)]
    true.pi <- rep(NA, sim$Nbask)
    true.pi[which(csv.pi=="pi.null")] <- sim$pi.null
    true.pi[which(csv.pi=="pi.alt")] <- sim$pi.alt
    true.pi[which(!csv.pi %in% c("pi.null", "pi.alt"))] <- as.numeric(csv.pi[which(!csv.pi %in% c("pi.null", "pi.alt"))])
  }

  #----calculate power----
  OC<-trial(sim,design,true.pi=true.pi)
  
  for(i in 1:length(OC)){
    #Storage the tables
    Atitl<-sapply(c("power.OC","p.upper","mean","sd","cd","R","pi","lc","md","uc","exectime","n","x"),function(x){paste(c(x,i),collapse="")})
    OC.save<-list(OC[[i]]$ps,OC[[i]]$pu,OC[[i]]$mn,
                  OC[[i]]$sd,OC[[i]]$cd,OC[[i]]$R,true.pi,
                  OC[[i]]$lc,OC[[i]]$md,OC[[i]]$uc,
                  OC[[i]]$et,OC[[i]]$n,OC[[i]]$x)#1-2,1-3,2-3,1-4,2-4,...
    write_custom(sim$dir.name,Atitl,OC.save,sub2=sim$scen)
    
    #summaries
    if (strsplit(design,"_")[[1]][1]!="Adaptive-Lasso") {
    power<-result_power(ps=OC[[i]]$ps,sim$Nbask,sim$Nsim,sim$q)
    terminate<-result_terminate(ps=OC[[i]]$pu,sim$Nbask,sim$Nsim,sim$q.int)
    fwerror<-sum(rowSums(as.matrix(OC[[i]]$ps[, true.pi == sim$pi.null]) > sim$q, na.rm = T) > 0) / sim$Nsim
    discovery<-sum(rowSums(as.matrix(OC[[i]]$ps[, true.pi != sim$pi.null]) > sim$q, na.rm = T) > 0) / sim$Nsim
    } else {
      power<-result_power(ps=OC[[i]]$mn,sim$Nbask,sim$Nsim,sim$q)
      terminate<-result_terminate(ps=OC[[i]]$mn,sim$Nbask,sim$Nsim,sim$q.int)
      fwerror<-sum(rowSums(as.matrix(OC[[i]]$mn[, true.pi == sim$pi.null]) > sim$q, na.rm = T) > 0) / sim$Nsim
      discovery<-sum(rowSums(as.matrix(OC[[i]]$mn[, true.pi != sim$pi.null]) > sim$q, na.rm = T) > 0) / sim$Nsim
    }
    estimate<-list(name=colnames(OC[[i]]$mn),mean=colMeans(OC[[i]]$mn, na.rm=T),sd=colMeans(OC[[i]]$sd, na.rm=T))
    bias<-result_bias(pmean=OC[[i]]$mn[,1:sim$Nbask],true.pi,sim$Nsim)
    mean_et<-mean(OC[[i]]$et)

    #Storage the summaries
    Atitl<-sapply(c("power","terminate","fwerror","discovery","estimate","bias","mean_exectime"),function(x){paste(c(x,i),collapse="")})
    OC.save<-list(power,terminate,fwerror,discovery,estimate,bias,mean_et)
    write_custom(sim$dir.name,Atitl,OC.save,sub2=sim$scen)
  }
}



