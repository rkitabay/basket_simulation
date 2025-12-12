scenario_custom1<-function(Ascen,sim,design){
  find_dir(sim$dir.name$prt.dir,sim$dir.name$dir)
  write_set(prt.dir=sim$dir.name$prt.dir,dir=sim$dir.name$dir,title="_setting",data=sim,subtitle="sim")
  
  sim["q"]<-simulation_q(sim,design)
  for(i in 0:(length(Ascen)-1)){
    sim["scen"]<-i
    sim["scennm"]<-Ascen[i+1]
    simulation_power(sim,design)
    invisible(replicate(20, gc()))
  }
}
scenario_custom2<-function(Ascen,sim,design){
  find_dir(sim$dir.name$prt.dir,sim$dir.name$dir)
  write_set(prt.dir=sim$dir.name$prt.dir,dir=sim$dir.name$dir,title="_setting",data=sim,subtitle="sim")
  
  sim["q"]<-as.numeric(read_result(prt.dir=sim$dir.name$prt.dir,dir=sim$dir.name$dir,
                          title="q",subtitle="output",sep=","))
  for(i in 0:(length(Ascen)-1)){
    sim["scen"]<-i
    sim["scennm"]<-Ascen[i+1]
    simulation_power(sim,design)
    invisible(replicate(20, gc()))
  }
}

random_scenario<-function(sim,design){
  find_dir(sim$dir.name$prt.dir,sim$dir.name$dir)
  dir.name2<-data.frame(prt.dir=paste(c(sim$dir.name$prt.dir,"/data"),collapse=""),dir=sim$dir.name$dir)
  find_dir(dir.name2$prt.dir,dir.name2$dir)
  write_set(prt.dir=dir.name2$prt.dir,dir=dir.name2$dir,title="_setting",data=sim,subtitle="sim")
  
  for(s1 in 1:sim$Nscen){
    
    scen.rand<-create_pi_random(sim$Nbask,sim$seed,s1,sim$pi.null,sim$pi.alt)
    ef<-scen.rand$Aeff
    true.pi<-scen.rand$pi
    
    #----calculate power----
    OC<-trial(sim,design,true.pi=true.pi)
    
    for(i in 1:length(OC)){
      #Storage the tables
      Atitl<-sapply(c("power.OC","p.upper","mean","sd","cd","R","pi","lc","md","uc","exectime","n","x"),function(x){paste(c(x,i),collapse="")})
      OC.save<-list(OC[[i]]$ps,OC[[i]]$pu,OC[[i]]$mn,OC[[i]]$sd,OC[[i]]$cd,OC[[i]]$R,scen.rand,
                    OC[[i]]$lc,OC[[i]]$md,OC[[i]]$uc,OC[[i]]$et,OC[[i]]$n,OC[[i]]$x)#1-2,1-3,2-3,1-4,2-4,...
      write_custom(dir.name2,Atitl,OC.save,sub2=s1)
      
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
      write_custom(dir.name2,Atitl,OC.save,sub2=s1)
    }
    
  }
}