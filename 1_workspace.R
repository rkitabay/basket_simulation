rm(list = ls())
cat("\014")

library(R2jags)
library(BCHM)
library(mbend)
library(BayesFactor)
library(data.table)
library(kdensity)
library(cluster)
library(entropy)
library(philentropy)
library(basket)
library(boot)
library(foreach)
library(partitions)
library(svMisc)
library(truncdist)
library(gmp)
library(glmnet)
library(gtools)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"independent_a1_b1"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"BHM_a0.01_b0.01"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"BHM_a2_b2"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"BHM_a2_b20"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"UPSiDe-JS_01_M48-96_r0.01"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"UPSiDe-JS_01_M24-48_r0.01"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"Cluster-BHM_psi0.5_a0.01_b0.01"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"Cluster-BHM_psi0.5_a2_b2"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"Cluster-BHM_psi0.5_a2_b20"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"local-MEM_delta0_a1_b1"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"local-MEM_delta2_a1_b1"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")
dyn.load("fn_RoBoT_MCMC.dll")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"RoBoT_a01_b01_a11_b11_alpha2_mu0_sigma2_tau0_xi1"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")
dyn.load("fn_RoBoT_MCMC.dll")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"RoBoT_a01_b01_a11_b11_alpha0.1_mu0_sigma2_tau0_xi1"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")
dyn.load("fn_RoBoT_MCMC.dll")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)

design<-"RoBoT_a01_b01_a11_b11_alpha10_mu0_sigma2_tau0_xi1"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)
options(nwarnings = 1e+05)

design<-"Adaptive-Lasso_lambda0.0018_gamma1.5"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)
options(nwarnings = 1e+05)

design<-"Adaptive-Lasso_lambda0.003_gamma1.5"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

# -----------------------------------------------
rm(list = ls())
cat("\014")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("2_scenario.R",encoding="UTF-8")
source("3_simulation.R",encoding="UTF-8")
source("4_trial.R",encoding="UTF-8")
source("5_analysis.R",encoding="UTF-8")
source("6_library.R",encoding="UTF-8")

sim<-list(
  seed=123,
  Npat=NULL,
  Nbask=6,
  pi.null=0.10,
  pi.alt=0.35,
  q.int=0.05,
  TIE=0.05,
  Nsim=1000,
  dir.name=NULL,
  exs.list=NULL
)

sim$Npat<-list(method="poisson",max=c(25, 50),ratio=rep(1,sim$Nbask),
               nfix=matrix(c(12, 24),nrow=2,ncol=sim$Nbask,byrow=T))

prt.dir<-dir_name(sim)
options(warn = 2)
options(nwarnings = 1e+05)

design<-"Adaptive-Lasso_lambda0.0018_gamma0.5"
sim$dir.name<-data.frame(dir=design,prt.dir=prt.dir)
system.time(scenario_custom1(Ascen=c(0:sim$Nbask, "qc"),sim,design))

