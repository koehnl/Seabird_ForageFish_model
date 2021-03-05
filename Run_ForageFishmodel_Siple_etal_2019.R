# Running the Forage fish model

library(plyr)
library(here)
# Set directories
setwd(here::here())
subDir <- "Anchovy" # Name of ff type

Nyear = 1000
years.test = Nyear + 2 # use to be +1 but since doing average biomass in a year (so need year1+) and need 1001 years to match with 
# seabird model, moved up to 2. Does change values slightly because messes with set.seed/randomization
#years.test = 250
nsims = 100 # 100 before
tim.params = list(sigma0 = 0.2,tau0 = 0.1)
tau1 = (1/tim.params$tau0^2 + 1/tim.params$sigma0^2)^(-0.5)
sig.s = 0.3 # 
R0.sens = NA #NO DYNAMIC R0 anymore-- ignore


# Load packages
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

# Load rev devs generating fxn, MSE main model, estimator fxns
#toplot=FALSE      # Don't plot examples of recruitement trajectories
source(file.path("ForageFishModel/Recruitment/GenerateDevs.R")) 
source(file.path("ForageFishModel/Estimators/CalcFTrue.R"))
source(file.path("ForageFishModel/HCR_Trajectory_LEK_NEW_2timesteps.R")) 
source(file.path("ForageFishModel/Estimators/Estimators.R"))
source(file.path("ForageFishModel/Run/generate_M.R"))
source(file.path("ForageFishModel/Plots/Megsieggradar.R"))


# Load control files & set parameter values
#     Type     TargetSD
# 1  Sardine 0.5117551
# 2  Anchovy 0.3221512
# 3 Menhaden 0.3009440

# Anchovy/Herring

if(subDir == "Sardine"){
  source(file.path("ForageFishModel/Ctl/Sardine_LHControl.R"))
  source(file.path("ForageFishModel/Ctl/Sardine_FisheryControl.R"))
  # Sardine params
  recruit.sd <- 0.6
  recruit.rho <- 0.9
}
# Anchovy/Herring
if(subDir == "Anchovy"){
  source(file.path("ForageFishModel/Ctl/Anchovy_LHControl.R"))
  source(file.path("ForageFishModel/Ctl/Anchovy_FisheryControl.R"))
  #Anchovy recruitment dev params
  recruit.sd <- 0.6
  recruit.rho <- 0.5
}

# Menhaden
if(subDir == "Menhaden"){
  source(file.path("ForageFishModel/Ctl/Menhaden_LHControl.R"))
  source(file.path("ForageFishModel/Ctl/Menhaden_FisheryControl.R"))
  #Menhaden recruitment dev params
  recruit.sd <- 0.8
  recruit.rho <- 0.2
}

# Create list of matrices with error for the delayed detection scenario. This is important because the random seed will offset otherwise.
set.seed(123)
tim.rands.list <- list() #for all ensuring random values
n.ages = length(lh.test$ages)
tim.inits.vec <- rnorm(n.ages,0,tim.params$sigma0)  # just for initial values
for(sim in 1:nsims){
  tim.mat <- matrix(NA,nrow=n.ages,ncol=years.test)
  for(i in 1:years.test){
    tim.mat[,i] <- rnorm(n.ages,0,tau1)
  }
  tim.rands.list[[sim]] <- tim.mat
}

# Create errors for AC scenario
set.seed(123) #original set.seed was 123
curly.phi.mat <- matrix(NA, nrow = nsims,ncol = years.test)
for(sim in 1:nsims){
  curly.phi.mat[sim,] <- rnorm(years.test,0,sig.s)
}

# Load harvest rules
source(file.path("ForageFishModel/Control Rules/smith_oceana.R"))
source(file.path("ForageFishModel/Control Rules/cfp.R"))
source(file.path("ForageFishModel/Control Rules/hockey-stick.R"))
source(file.path("ForageFishModel/Control Rules/trend-based-rule.R"))


# If a results folder doesn't exist already, create one!
# resultsfolder <- paste(subDir,Sys.Date(),sep="")
# dir.create(file.path(resultsdir, resultsfolder))
# setwd(file.path(resultsdir, resultsfolder))

# Scenarios
h = c(0.6) # should be 0.6, conservative 
obs.error.type = c("AC")

HCR = c("nofish","constF", "constFhi","msycatch" ,"C1", "C2", "C3")#,"constFdelay", "C1delay", "constFhi", "constFhidelay") # Took out trend because it was unrealistic-- but using trend in CPUE as adjustment (data-poor method) might be a good idea!
M.type = c("constant") # took out "regimeshift" and "time-varying" to save time but can be added back in for sensitivity

scenarios <- expand.grid(h,obs.error.type,HCR,recruit.sd,recruit.rho,M.type)
colnames(scenarios) <- c("h","obs.error.type","HCR","recruit.sd","recruit.rho","M.type")
nscenarios <- nrow(scenarios)
scenarios$scenario <- 1:nscenarios #Label them so it's easier to find/index em later

write.table(scenarios,file = "Scenario_Table.txt")

nofish <- msycatch <- constF <- constFhi <- C1 <- C2 <- C3 <- #<-constFdelay  <- C1delay <- constFhi <- constFhidelay <-
  list(popn = array(NA, c(nsims,n.ages,years.test)),
       popbio = array(NA, c(nsims, n.ages, years.test)),
    biomass.oneplus.true=matrix(nrow = nsims,ncol = years.test),
       biomass.total.obs = matrix(nrow = nsims, ncol = years.test),
       total.catch=matrix(nrow = nsims,ncol = years.test),
       fishing= matrix(nrow = nsims,ncol = years.test),
       intended.f=matrix(nrow = nsims,ncol = years.test),
       rec= matrix(nrow = nsims,ncol = years.test),
       biomass.oneplus.obs = matrix(nrow = nsims,ncol = years.test),
       biomass.total.true = matrix(nrow = nsims,ncol = years.test),
       no.fishing.tb = matrix(nrow = nsims,ncol = years.test),
       penguinfood = array(NA, c(2, years.test, nsims)),
    avgbiomass = matrix(nrow = nsims, ncol = years.test-1))

avgbio <- function(biomass, M, Fmort) { 
  biomass_avg = rep(NA, length = years.test-1)
  a = length(ages.test)
  bio_avg_age = matrix(NA, nrow = a, ncol = years.test-1)
  biomass_avg[1] = sum(biomass[2:a,1])
  bio_avg_age[,1] = biomass[1]
  for(i in 2:(years.test-1)) {
    bio_t = biomass[,i]
    bio_t1 = biomass[,i+1]
    
    gamma = log(bio_t1[2:a]/bio_t[1:(a-1)])
    # gamma = log(Ba+1, t+1 / Ba, t) = g - M - F
    #gamma = g - M - Fmort[i]
    print(gamma)
    temp = bio_t[1:(a-1)]*(exp(gamma)-1) / gamma
  #gammaplus = (bio_t1[a] - (bio_t1[a-1]*exp(-M-Fmort[i])))/(bio_t[a])
    # gamma_plus = (B{plus, t+1} - B{plus -1, t+1}* exp(-M - F)) / B{plus, t+1} 
   # gammaplus = -M-Fmort[i]
   #temp2 = bio_t[a]*(exp(gammaplus)-1) / gammaplus
    biomass_avg[i] = sum(temp) + bio_t[a] - temp[1]#

    bio_avg_age[1:(a-1),i] = temp
    #bio_avg_age[a,i] = temp2
  }
  return(biomass_avg)
  #return(bio_avg_age)
}
# if run equation: Bavg = Bt (e(g-m-f)-1/ g-m-f) than Bavg between 10 and 18.22 is 13.7
# using my equations - 
biotest = matrix(c(10,18.22,2,4), nrow = 2, ncol = 2)
g = log(18.22/10)
10*(exp(g)-1) / g # matches what it should be if we knew g so is working for non
# plus age group 


#save model output
output = paste(getwd(), "/", sep="")

#run scenarios
for(s in 1:nscenarios){
  steepness = scenarios$h[s]
  obs.type <- scenarios$obs.error.type[s]
  HCR <- scenarios$HCR[s]
  recruit.sd = scenarios$recruit.sd[s]
  recruit.rho = scenarios$recruit.rho[s]
  M.type = scenarios$M.type[s]
  
  equilib = getEquilibriumConditions(lh = lh.test,fish = seq(0,5,by=.1),years = 150,steepness=steepness) # NO recruitment devs used in the equilibrium calculations, so don't need to embed in the loop
  const.f.rate = 0.5*equilib$Fmsy 
  no.fishing <- matrix(NA, nrow = nsims, ncol = years.test)
  #RNGkind(sample.kind = "Rounding")
  set.seed(123) # Start each round of sims at same random seed
  time.var.m <- NA # Base case: M constant 
  if(M.type == "timevar"){time.var.m <- rw.M(Mbar = lh.test$M, rho.m = 0.6, sigma.m = 0.2,n = years.test)}
  if(M.type == "regimeshift"){time.var.m <- regime.M(Mbar = lh.test$M,cutoff.yr = 201,n = years.test)}
  for (sim in 1:nsims){
    rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
    F0.Type <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                               hcr.type = "constF", const.f.rate = 0, steepness = steepness,obs.type = "noerror",equilib=equilib,R0.traj = R0.sens, 
                               tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
                               tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])$biomass.total.true # Only need to do this 1x for each simulation (not repeat for each CR) because the seed is the same and there is no fishing.
    no.fishing[sim,] <- F0.Type
  }
  
  if(HCR=="nofish"){
    #RNGkind(sample.kind = "Rounding")
    set.seed(123) # same seed
    for (sim in 1:nsims){
      const.f.rate_nofish <- 0
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.nofish <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                     hcr.type = "constF", const.f.rate = const.f.rate_nofish, steepness = steepness,obs.type = "noerror",equilib=equilib,
                                     R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
                                     tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      nofish[["popn"]][sim,,] <- expt.nofish$popn[1,,]
      nofish[["popbio"]][sim,,] <- expt.nofish$popbio[1,,]
      nofish[["biomass.oneplus.true"]][sim,] <- expt.nofish$biomass.oneplus.true
      nofish[["biomass.total.obs"]][sim,] <- expt.nofish$biomass.total.obs
      nofish[["total.catch"]][sim,] <- expt.nofish$total.catch
      nofish[["fishing"]][sim,] <- expt.nofish$fishing
      nofish[["intended.f"]][sim,] <- expt.nofish$intended.f   
      nofish[["rec"]][sim,] <- expt.nofish$rec
      nofish[["biomass.oneplus.obs"]][sim,] <- expt.nofish$biomass.oneplus.obs
      nofish[["biomass.total.true"]][sim,] <- expt.nofish$biomass.total.true
      nofish[["no.fishing.tb"]] <- no.fishing
      temp = nrow(expt.nofish$catch.at.age)
      nofish[["penguinfood"]][1,,sim] <- colSums(expt.nofish$popbio[1,2:temp,])
      nofish[["penguinfood"]][2,,sim] <- colSums(expt.nofish$popbio[2,2:temp,])
      nofish[["avgbiomass"]][sim, ] = avgbio(expt.nofish$popbio[1,,],M = lh.test$M, Fmort = expt.nofish$fishing)
    }
    save(nofish,file=paste(output, subDir,s,"nofish",".RData",sep="_"))
  }
  

  if(HCR=="constF"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      #const.f.rate <- 0.6
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.constF <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                     hcr.type = "constF", const.f.rate =0.25*equilib$Fmsy, steepness = steepness,obs.type = "noerror",equilib=equilib,
                                     R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
                                     tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      constF[["popn"]][sim,,] <- expt.constF$popn[1,,]
      constF[["popbio"]][sim,,] <- expt.constF$popbio[1,,]
      constF[["biomass.oneplus.true"]][sim,] <- expt.constF$biomass.oneplus.true
      constF[["biomass.total.obs"]][sim,] <- expt.constF$biomass.total.obs
      constF[["total.catch"]][sim,] <- expt.constF$total.catch
      constF[["fishing"]][sim,] <- expt.constF$fishing
      constF[["intended.f"]][sim,] <- expt.constF$intended.f   
      constF[["rec"]][sim,] <- expt.constF$rec
      constF[["biomass.oneplus.obs"]][sim,] <- expt.constF$biomass.oneplus.obs
      constF[["biomass.total.true"]][sim,] <- expt.constF$biomass.total.true
      constF[["no.fishing.tb"]] <- no.fishing
      temp = nrow(expt.constF$catch.at.age)
      constF[["penguinfood"]][1,,sim] <- colSums(expt.constF$popbio[1,2:temp,])
      constF[["penguinfood"]][2,,sim] <- colSums(expt.constF$popbio[2,2:temp,])
      constF[["avgbiomass"]][sim, ] = avgbio(expt.constF$popbio[1,,],M = lh.test$M, Fmort = expt.constF$fishing)
      
    }
    save(constF,file=paste(output, subDir,s,"constF",".RData",sep="_"))
  }
  
  if(HCR=="C1"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
      expt.c1 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                 hcr.type = "C1",equilib = equilib,steepness=steepness,obs.type = "noerror",R0.traj = R0.sens, tim.params = tim.params,
                                 time.var.m = time.var.m, sig.s = sig.s, props = props,
                                 tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])

      C1[["popn"]][sim,,] <- expt.c1$popn[1,,]
      C1[["biomass.oneplus.true"]][sim,] <- expt.c1$biomass.oneplus.true
      C1[["biomass.total.obs"]][sim,] <- expt.c1$biomass.total.obs
      C1[["total.catch"]][sim,] <- expt.c1$total.catch
      C1[["fishing"]][sim,] <- expt.c1$fishing
      C1[["intended.f"]][sim,] <- expt.c1$intended.f
      C1[["rec"]][sim,] <- expt.c1$rec
      C1[["biomass.oneplus.obs"]][sim,] <- expt.c1$biomass.oneplus.obs
      C1[["biomass.total.true"]][sim,] <- expt.c1$biomass.total.true
      C1[["no.fishing.tb"]] <- no.fishing
      temp = nrow(expt.nofish$catch.at.age)
      C1[["penguinfood"]][1,,sim] <- colSums(expt.c1$popbio[1,2:temp,])
      C1[["penguinfood"]][2,,sim] <- colSums(expt.c1$popbio[2,2:temp,])
      C1[["avgbiomass"]][sim, ] = avgbio(expt.c1$popbio[1,,],M = lh.test$M, Fmort = expt.c1$fishing)
      
    }
    save(C1,file=paste(output, subDir,s,"C1",".RData",sep="_"))
  }

#   if(HCR=="constFdelay"){
#     set.seed(123) # same seed
#     for (sim in 1:nsims){
#       #const.f.rate <- 0.6
#       rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
#       expt.constFdelay <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
#                                      hcr.type = "constF", const.f.rate = 0.25, steepness = steepness,obs.type = "Tim",equilib=equilib,
#                                      R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
#                                      tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
#       
#       constFdelay[["biomass.oneplus.true"]][sim,] <- expt.constFdelay$biomass.oneplus.true
#       constFdelay[["biomass.total.obs"]][sim,] <- expt.constFdelay$biomass.total.obs
#       constFdelay[["total.catch"]][sim,] <- expt.constFdelay$total.catch
#       constFdelay[["fishing"]][sim,] <- expt.constFdelay$fishing
#       constFdelay[["intended.f"]][sim,] <- expt.constFdelay$intended.f   
#       constFdelay[["rec"]][sim,] <- expt.constFdelay$rec
#       constFdelay[["biomass.oneplus.obs"]][sim,] <- expt.constFdelay$biomass.oneplus.obs
#       constFdelay[["biomass.total.true"]][sim,] <- expt.constFdelay$biomass.total.true
#       constFdelay[["no.fishing.tb"]] <- no.fishing
#       temp = nrow(expt.constFdelay$catch.at.age)
#       constFdelay[["penguinfood"]][1,,sim] <- colSums(expt.constFdelay$popbio[1,2:temp,])
#       constFdelay[["penguinfood"]][2,,sim] <- colSums(expt.constFdelay$popbio[2,2:temp,])
#     }
#     save(constFdelay,file=paste("All",s,"constFdelay",".RData",sep="_"))
#   }  
#   
#   if(HCR=="C1delay"){
#     set.seed(123) # same seed
#     for (sim in 1:nsims){
#       rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
#       expt.c1delay <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
#                                  hcr.type = "C1",equilib = equilib,steepness=steepness,obs.type = "Tim",R0.traj = R0.sens, tim.params = tim.params,
#                                  time.var.m = time.var.m, sig.s = sig.s, props = props, 
#                                  tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
#       
#       C1delay[["biomass.oneplus.true"]][sim,] <- expt.c1delay$biomass.oneplus.true
#       C1delay[["biomass.total.obs"]][sim,] <- expt.c1delay$biomass.total.obs
#       C1delay[["total.catch"]][sim,] <- expt.c1delay$total.catch
#       C1delay[["fishing"]][sim,] <- expt.c1delay$fishing
#       C1delay[["intended.f"]][sim,] <- expt.c1delay$intended.f    
#       C1delay[["rec"]][sim,] <- expt.c1delay$rec
#       C1delay[["biomass.oneplus.obs"]][sim,] <- expt.c1delay$biomass.oneplus.obs
#       C1delay[["biomass.total.true"]][sim,] <- expt.c1delay$biomass.total.true
#       C1delay[["no.fishing.tb"]] <- no.fishing
#       #temp = nrow(expt.nofish$catch.at.age)
#       C1delay[["penguinfood"]][1,,sim] <- colSums(expt.c1delay$popbio[1,2:temp,])
#       C1delay[["penguinfood"]][2,,sim] <- colSums(expt.c1delay$popbio[2,2:temp,])
#     }
#     save(C1delay,file=paste("All",s,"C1delay",".RData",sep="_")) 
#   }
  if(HCR=="C2"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
      expt.C2 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                 hcr.type = "C2",equilib = equilib,steepness=steepness,obs.type = "noerror",R0.traj = R0.sens, tim.params = tim.params,
                                 time.var.m = time.var.m, sig.s = sig.s, props = props,
                                 tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      C2[["popn"]][sim,,] <- expt.C2$popn[1,,]
      C2[["biomass.oneplus.true"]][sim,] <- expt.C2$biomass.oneplus.true
      C2[["biomass.total.obs"]][sim,] <- expt.C2$biomass.total.obs
      C2[["total.catch"]][sim,] <- expt.C2$total.catch
      C2[["fishing"]][sim,] <- expt.C2$fishing
      C2[["intended.f"]][sim,] <- expt.C2$intended.f
      C2[["rec"]][sim,] <- expt.C2$rec
      C2[["biomass.oneplus.obs"]][sim,] <- expt.C2$biomass.oneplus.obs
      C2[["biomass.total.true"]][sim,] <- expt.C2$biomass.total.true
      C2[["no.fishing.tb"]] <- no.fishing
      temp = nrow(expt.nofish$catch.at.age)
      C2[["penguinfood"]][1,,sim] <- colSums(expt.C2$popbio[1,2:temp,])
      C2[["penguinfood"]][2,,sim] <- colSums(expt.C2$popbio[2,2:temp,])
      C2[["avgbiomass"]][sim, ] = avgbio(expt.C2$popbio[1,,],M = lh.test$M, Fmort = expt.C2$fishing)
      
    }
    save(C2,file=paste(output, subDir,s,"C2",".RData",sep="_"))
  }
  
  if(HCR=="C3"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd)
      expt.C3 <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                 hcr.type = "C3",equilib = equilib,steepness=steepness,obs.type = "noerror",R0.traj = R0.sens, tim.params = tim.params,
                                 time.var.m = time.var.m, sig.s = sig.s, props = props,
                                 tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      C3[["popn"]][sim,,] <- expt.C3$popn[1,,]
      C3[["biomass.oneplus.true"]][sim,] <- expt.C3$biomass.oneplus.true
      C3[["biomass.total.obs"]][sim,] <- expt.C3$biomass.total.obs
      C3[["total.catch"]][sim,] <- expt.C3$total.catch
      C3[["fishing"]][sim,] <- expt.C3$fishing
      C3[["intended.f"]][sim,] <- expt.C3$intended.f
      C3[["rec"]][sim,] <- expt.C3$rec
      C3[["biomass.oneplus.obs"]][sim,] <- expt.C3$biomass.oneplus.obs
      C3[["biomass.total.true"]][sim,] <- expt.C3$biomass.total.true
      C3[["no.fishing.tb"]] <- no.fishing
      temp = nrow(expt.nofish$catch.at.age)
      C3[["penguinfood"]][1,,sim] <- colSums(expt.C3$popbio[1,2:temp,])
      C3[["penguinfood"]][2,,sim] <- colSums(expt.C3$popbio[2,2:temp,])
      C3[["avgbiomass"]][sim, ] = avgbio(expt.C3$popbio[1,,],M = lh.test$M, Fmort = expt.C3$fishing)
      
    }
    save(C3,file=paste(output, subDir,s,"C3",".RData",sep="_"))
  }
  
  if(HCR=="constFhi"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      #const.f.rate <- 0.6
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.constFhi <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                     hcr.type = "constF", const.f.rate = 0.5*equilib$Fmsy , steepness = steepness,obs.type = "noerror",equilib=equilib, # constant F hi should normally be 0.5 MSY, testing full MSY on generalist
                                     R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
                                     tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      constFhi[["popn"]][sim,,] <- expt.constFhi$popn[1,,]
      constFhi[["biomass.oneplus.true"]][sim,] <- expt.constFhi$biomass.oneplus.true
      constFhi[["biomass.total.obs"]][sim,] <- expt.constFhi$biomass.total.obs
      constFhi[["total.catch"]][sim,] <- expt.constFhi$total.catch
      constFhi[["fishing"]][sim,] <- expt.constFhi$fishing
      constFhi[["intended.f"]][sim,] <- expt.constFhi$intended.f   
      constFhi[["rec"]][sim,] <- expt.constFhi$rec
      constFhi[["biomass.oneplus.obs"]][sim,] <- expt.constFhi$biomass.oneplus.obs
      constFhi[["biomass.total.true"]][sim,] <- expt.constFhi$biomass.total.true
      constFhi[["no.fishing.tb"]] <- no.fishing
      #temp = nrow(expt.constF$catch.at.age)
      constFhi[["penguinfood"]][1,,sim] <- colSums(expt.constFhi$popbio[1,2:temp,])
      constFhi[["penguinfood"]][2,,sim] <- colSums(expt.constFhi$popbio[2,2:temp,])
      constFhi[["avgbiomass"]][sim, ] = avgbio(expt.constFhi$popbio[1,,],M = lh.test$M, Fmort = expt.constFhi$fishing)
      
    }
    save(constFhi,file=paste(output, subDir,s,"constFhi",".RData",sep="_"))
  }
  
  if(HCR=="msycatch"){
    set.seed(123) # same seed
    for (sim in 1:nsims){
      #const.f.rate <- 0.6
      rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
      expt.msy <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
                                       hcr.type = "constF", const.f.rate = equilib$Fmsy , steepness = steepness,obs.type = "noerror",equilib=equilib, # constant F hi should normally be 0.5 MSY, testing full MSY on generalist
                                       R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
                                       tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
      
      msycatch[["popn"]][sim,,] <- expt.msy$popn[1,,]
      msycatch[["biomass.oneplus.true"]][sim,] <- expt.msy$biomass.oneplus.true
      msycatch[["biomass.total.obs"]][sim,] <- expt.msy$biomass.total.obs
      msycatch[["total.catch"]][sim,] <- expt.msy$total.catch
      msycatch[["fishing"]][sim,] <- expt.msy$fishing
      msycatch[["intended.f"]][sim,] <- expt.msy$intended.f   
      msycatch[["rec"]][sim,] <- expt.msy$rec
      msycatch[["biomass.oneplus.obs"]][sim,] <- expt.msy$biomass.oneplus.obs
      msycatch[["biomass.total.true"]][sim,] <- expt.msy$biomass.total.true
      msycatch[["no.fishing.tb"]] <- no.fishing
      #temp = nrow(expt.constF$catch.at.age)
      msycatch[["penguinfood"]][1,,sim] <- colSums(expt.msy$popbio[1,2:temp,])
      msycatch[["penguinfood"]][2,,sim] <- colSums(expt.msy$popbio[2,2:temp,])
      msycatch[["avgbiomass"]][sim, ] = avgbio(expt.msy$popbio[1,,],M = lh.test$M, Fmort = expt.msy$fishing)
      
    }
    save(msycatch,file=paste(output, subDir,s,"msycatch",".RData",sep="_"))
  }
#   if(HCR=="constFhidelay"){
#     set.seed(123) # same seed
#     for (sim in 1:nsims){
#       #const.f.rate <- 0.6
#       rec.dev.test  <-  generate.devs(N = years.test,rho = recruit.rho,sd.devs = recruit.sd) 
#       expt.constFhidelay <- calc.trajectory(lh = lh.test,obs.cv = 1.2, init = init.test, rec.dev = rec.dev.test, F0 = F0.test, cr = cr.test, years = years.test,
#                                           hcr.type = "constF", const.f.rate = const.f.rate, steepness = steepness,obs.type = "Tim",equilib=equilib,
#                                           R0.traj = R0.sens, tim.params = tim.params,time.var.m = time.var.m, sig.s = sig.s, props = props,
#                                           tim.rand.inits = tim.inits.vec, tim.rands = tim.rands.list[[sim]],curly.phi.vec = curly.phi.mat[sim,])
#       
#       constFhidelay[["biomass.oneplus.true"]][sim,] <- expt.constFhidelay$biomass.oneplus.true
#       constFhidelay[["biomass.total.obs"]][sim,] <- expt.constFhidelay$biomass.total.obs
#       constFhidelay[["total.catch"]][sim,] <- expt.constFhidelay$total.catch
#       constFhidelay[["fishing"]][sim,] <- expt.constFhidelay$fishing
#       constFhidelay[["intended.f"]][sim,] <- expt.constFhidelay$intended.f   
#       constFhidelay[["rec"]][sim,] <- expt.constFhidelay$rec
#       constFhidelay[["biomass.oneplus.obs"]][sim,] <- expt.constFhidelay$biomass.oneplus.obs
#       constFhidelay[["biomass.total.true"]][sim,] <- expt.constFhidelay$biomass.total.true
#       constFhidelay[["no.fishing.tb"]] <- no.fishing
#       #temp = nrow(expt.constFdelay$catch.at.age)
#       constFhidelay[["penguinfood"]][1,,sim] <- colSums(expt.constFhidelay$popbio[1,2:temp,])
#       constFhidelay[["penguinfood"]][2,,sim] <- colSums(expt.constFhidelay$popbio[2,2:temp,])
#     }
#     save(constFhidelay,file=paste("All",s,"constFhidelay",".RData",sep="_"))
#   }  
}


# EXTRA ADDITIONAL CODE NOT NEEDED 
biomass = constF$popn[1,,]

wtuse = lh.test$w.at.age
a= length(lh.test$ages)
g = log(wtuse[2:(a)]/ wtuse[1:(a-1)])
g = c(g,0) # right isn't g in the plus group 0?

example =  constF$popbio#need for each sim
biohalf <- function(biomass, g, M, Fmort) {
  biohalfway = array(NA, c(2, nsims,years.test))
  for(j in 1:nsims) {
   bio = biomass[j,,]
   biohalfway[,,1] = sum(bio[2:a,1])
   for(i in 2:(years.test-1)) {
      change = g - M - Fmort[j,i]
      
      first = bio[ ,i-1] * ((exp(change*0.5)-1)/ change)
      biohalfway[1,j,i]  = sum(first)
      tempsecond = (exp(change) - exp(change*0.5)) / change
      
      second = bio[ ,i-1] * tempsecond
      biohalfway[2,j,i] = sum(second)
    }
  }
  return(biohalfway)
}

attempt = biohalf(example, g = g, M = lh.test$M, Fmort = constF$fishing)

expt.constF$wts[,900:910]
log(expt.constF$wts[2,2]/expt.constF$wts[1,1])

x = avgbio(constF$popn[1,,])
plot(x[1,], type = 'l')
lines(colSums(biomass[1,,]), col = "red", lty = 2)


color = viridis(6, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
nofishrun = rowMeans(nofish$penguinfood[1,,])
constFrun = rowMeans(constF$penguinfood[1,,])
constFhirun = rowMeans(constFhi$penguinfood[1,,])
C1run =rowMeans(C1$penguinfood[1,,])
C2run =rowMeans(C2$penguinfood[1,,])
C3run =rowMeans(C3$penguinfood[1,,])
plot(nofishrun[201:1001]/nofishrun[201:1001], type = 'l', ylim = c(0, 1), col = color[1],
     xlab = "Time", ylab = "Sardine biomass proportion (compared to unfished)")
lines(constFrun[201:1001]/nofishrun[201:1001], col = color[2])
lines(constFhirun[201:1001]/nofishrun[201:1001], col = color[3])
lines(C1run[201:1001]/nofishrun[201:1001], col = color[4])
lines(C2run[201:1001]/nofishrun[201:1001], col = color[5])
lines(C3run[201:1001]/nofishrun[201:1001], col = color[6])
legend("bottomright", legend = c("No F", "Const. F", "Const. Fmsy",
                                 "cut-off 0.5 msy", "cut-off low Blim", "cut-off F msy"),
       col = color[1:6], lty = 1,  ncol = 3, bty = FALSE)

color = viridis(6, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
nofishrun = rowMins(nofish$penguinfood[1,,])
constFrun = rowMins(constF$penguinfood[1,,])
constFhirun = rowMins(constFhi$penguinfood[1,,])
C1run =rowMins(C1$penguinfood[1,,])
C2run =rowMins(C2$penguinfood[1,,])
C3run =rowMins(C3$penguinfood[1,,])
plot(nofishrun[201:1001], type = 'l', ylim = c(0,5000), col = color[1],
     xlab = "Time", ylab = "biomass minimum")
lines(constFrun[201:1001], col = color[2])
lines(constFhirun[201:1001], col = color[3])
lines(C1run[201:1001], col = color[4])
lines(C2run[201:1001], col = color[5])
lines(C3run[201:1001], col = color[6])
legend("bottomright", legend = c("No F", "Const. F", "Const. Fmsy",
                                 "cut-off 0.5 msy", "cut-off low Blim", "cut-off F msy"),
       col = color[1:6], lty = 1,  ncol = 3, bty = FALSE)

color = viridis(6, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
nofishrun = rowMaxs(nofish$penguinfood[1,,])
constFrun = rowMaxs(constF$penguinfood[1,,])
constFhirun = rowMaxs(constFhi$penguinfood[1,,])
C1run =rowMaxs(C1$penguinfood[1,,])
C2run =rowMaxs(C2$penguinfood[1,,])
C3run =rowMaxs(C3$penguinfood[1,,])
plot(nofishrun[201:1001], type = 'l', ylim = c(8000,40000), col = color[1],
     xlab = "Time", ylab = "biomass maximum")
lines(constFrun[201:1001], col = color[2])
lines(constFhirun[201:1001], col = color[3])
lines(C1run[201:1001], col = color[4])
lines(C2run[201:1001], col = color[5])
lines(C3run[201:1001], col = color[6])
legend("bottomright", legend = c("No F", "Const. F", "Const. Fmsy",
                                 "cut-off 0.5 msy", "cut-off low Blim", "cut-off F msy"),
       col = color[1:6], lty = 1,  ncol = 3, bty = FALSE)

nofishrun = rowMedians(nofish$penguinfood[1,,])
constFrun = rowMedians(constF$penguinfood[1,,])
constFhirun = rowMedians(constFhi$penguinfood[1,,])
C1run =rowMedians(C1$penguinfood[1,,])
C2run =rowMedians(C2$penguinfood[1,,])
C3run =rowMedians(C3$penguinfood[1,,])
plot(nofishrun[201:1001]/nofishrun[201:1001], type = 'l', ylim = c(0, 1), col = color[1],
     xlab = "Time", ylab = "Anchovy bio median (proportion)")
lines(constFrun[201:1001]/nofishrun[201:1001], col = color[2])
lines(constFhirun[201:1001]/nofishrun[201:1001], col = color[3])
lines(C1run[201:1001]/nofishrun[201:1001], col = color[4])
lines(C2run[201:1001]/nofishrun[201:1001], col = color[5])
lines(C3run[201:1001]/nofishrun[201:1001], col = color[6])
legend("bottomright", legend = c("No F", "Const. F", "Const. Fmsy",
                                 "cut-off 0.5 msy", "cut-off low Blim", "cut-off F msy"),
       col = color[1:6], lty = 1,  ncol = 3, bty = FALSE)

nofishrun = rowVars(nofish$penguinfood[1,,])
constFrun = rowVars(constF$penguinfood[1,,])
constFhirun = rowVars(constFhi$penguinfood[1,,])
C1run =rowVars(C1$penguinfood[1,,])
C2run =rowVars(C2$penguinfood[1,,])
C3run =rowVars(C3$penguinfood[1,,])
plot(nofishrun[201:1001]/nofishrun[201:1001], type = 'l', ylim = c(-0.2, 1), col = color[1],
     xlab = "Time", ylab = "Anchovy Bio Variance (Proportion of unfished)")
lines(constFrun[201:1001]/nofishrun[201:1001], col = color[2])
lines(constFhirun[201:1001]/nofishrun[201:1001], col = color[3])
lines(C1run[201:1001]/nofishrun[201:1001], col = color[4])
lines(C2run[201:1001]/nofishrun[201:1001], col = color[5])
lines(C3run[201:1001]/nofishrun[201:1001], col = color[6])
legend("bottomright", legend = c("No F", "Const. F", "Const. Fmsy",
                                 "cut-off 0.5 msy", "cut-off low Blim", "cut-off F msy"),
       col = color[1:6], lty = 1,  ncol = 3, bty = FALSE)

# CV
cv <- function(matrix){
  CV = rep(NA, length.out = ncol(matrix))
  for(i in 1:ncol(matrix)) {
  mean = mean(matrix[,i])
  sd = sd(matrix[,i])
  CV[i] <- sd/mean
  }
  return(CV)
}

median(cv(constF$penguinfood[1,,])) # 0.66
median(cv(constFhi$penguinfood[1,,])) # 0.73
median(cv(nofish$penguinfood[1,,])) # 0.61

median(cv(C1$penguinfood[1,,])) # 0.57
median(cv(C2$penguinfood[1,,])) # 0.63
median(cv(C3$penguinfood[1,,])) # 0.54

nofishrun = rowMeans(nofish$penguinfood[1,,])/rowSds(nofish$penguinfood[1,,])
constFrun = rowMeans(constF$penguinfood[1,,])/rowSds(constF$penguinfood[1,,])
constFhirun = rowMeans(constFhi$penguinfood[1,,])/rowSds(constFhi$penguinfood[1,,])
C1run =rowMeans(C1$penguinfood[1,,])/rowSds(C1$penguinfood[1,,])
C2run =rowMeans(C2$penguinfood[1,,])/rowSds(C2$penguinfood[1,,])
C3run =rowMeans(C3$penguinfood[1,,])/rowSds(C3$penguinfood[1,,])
plot(nofishrun[201:1001]/nofishrun[201:1001], type = 'l', ylim = c(0.5, 1.6), col = color[1],
     xlab = "Time", ylab = "anchovy CV proportion")
lines(constFrun[201:1001]/nofishrun[201:1001], col = color[2])
lines(constFhirun[201:1001]/nofishrun[201:1001], col = color[3])
lines(C1run[201:1001]/nofishrun[201:1001], col = color[4])
lines(C2run[201:1001]/nofishrun[201:1001], col = color[5])
lines(C3run[201:1001]/nofishrun[201:1001], col = color[6])
legend("bottomright", legend = c("No F", "Const. F", "Const. Fmsy",
                                 "cut-off 0.5 msy", "cut-off low Blim", "cut-off F msy"),
       col = color[1:6], lty = 1,  ncol = 3, bty = FALSE)

minNo = apply(nofish$penguinfood[1,201:1001,], 2, min)
CIno = median(minNo)+c(-1.96,1.96)*sd(minNo)/sqrt(length(minNo))

minConst = apply(constF$penguinfood[1,201:1001,], 2, min)
CIConst = median(minConst)+c(-1.96,1.96)*sd(minConst)/sqrt(length(minConst))

minC1 = apply(C1$penguinfood[1,201:1001,], 2, min)
CIC1 = median(minC1)+c(-1.96,1.96)*sd(minC1)/sqrt(length(minC1))

minConsthi = apply(constFhi$penguinfood[1,201:1001,], 2, min)
CIConsthi = median(minConsthi)+c(-1.96,1.96)*sd(minConsthi)/sqrt(length(minConsthi))

minC2 = apply(C2$penguinfood[1,201:1001,], 2, min)
CIC2 = median(minC2)+c(-1.96,1.96)*sd(minC2)/sqrt(length(minC2))

minC3 = apply(C3$penguinfood[1,201:1001,], 2, min)
CIC3 = median(minC3)+c(-1.96,1.96)*sd(minC3)/sqrt(length(minC3))


mid = c(median(minC1)/median(minNo), median(minC2)/median(minNo), median(minC3)/median(minNo), median(minConst)/median(minNo), median(minConsthi)/median(minNo) )
high = c(mean(CIC1[1])/mean(CIno[1]), mean(CIC2[1])/mean(CIno[1]), mean(CIC3[1])/mean(CIno[1]), mean(CIConst[1])/mean(CIno[1]), mean(CIConsthi[1])/mean(CIno[1]) )
low = c(mean(CIC1[2])/mean(CIno[2]), mean(CIC2[2])/mean(CIno[2]), mean(CIC3[2])/mean(CIno[2]), mean(CIConst[2])/mean(CIno[2]), mean(CIConsthi[2])/mean(CIno[2]) )

par(mar = c(6,2,2,2))
x = c(1,2,3,4,5)
names = c("Cut-off, 0.5 Fmsy","Cut-off - Low limit", "Cut-off, Fmsy", "Constant: 0.25 Fmsy", "Constant: 0.5 Fmsy")
plot(mid, x, frame = F, xlab = "Minimum Forage Biomass Proportion", ylab = "", 
     xlim = c(0,1), axes = F, pch = 19, col = c("Blue"))
arrows(low, x, 
       high, x, length=0, angle=90, code=3, lwd = 1,col = c("Blue"))
#axis(2, at = c(1,2,3,4,5), labels = c("", "", "", "", ""), las = 2, pos = -0.1, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = 0.8)


# Proportion of years below some value? in a collapse, so <0.2 B0
thres = 0.2*equilib$B0
thresNo = apply(nofish$penguinfood[1,201:1001,], 2, function(x){ length(which(x < thres))/801})
thresConst = apply(constF$penguinfood[1,201:1001,], 2, function(x){ length(which(x < thres))/801})
thresConsthi = apply(constFhi$penguinfood[1,201:1001,], 2, function(x){ length(which(x < thres))/801})
thresC1 = apply(C1$penguinfood[1,201:1001,], 2, function(x){ length(which(x < thres))/801})
thresC2 = apply(C2$penguinfood[1,201:1001,], 2, function(x){ length(which(x < thres))/801})
thresC3 = apply(C3$penguinfood[1,201:1001,], 2, function(x){ length(which(x < thres))/801})
thresmean = c(mean(thresC1), mean(thresC2), mean(thresC3), mean(thresConst), mean(thresConsthi))
par(mar = c(6,2,2,2))
x = c(1,2,3,4,5)
names = c("Cut-off, 0.5 Fmsy","Cut-off - Low limit", "Cut-off, Fmsy", "Constant: 0.25 Fmsy", "Constant: 0.5 Fmsy")
plot(1-thresmean, x, frame = F, xlab = "Average prop years not collapsed", ylab = "", 
     xlim = c(0.8,1), axes = F, pch = 19, col = c("Blue"))
#axis(2, at = c(1,2,3,4,5), labels = c("", "", "", "", ""), las = 2, pos = -0.1, cex.axis = 1)
axis(1, at=c(0.8,0.85,0.9,0.95,1), cex.axis=1, tck=0.01, pos = 0.8)


par(mfrow = c(1,2))
par(mar = c(2,2,1,1))
par(oma = c(3,3,1,1))
plot(nofish$biomass.total.true[1,200:300], type = 'l', xlab= "", ylab = "", 
     ylim = c(0,1500000))
lines(constFhi$biomass.total.true[1,200:300], col = 'grey50', lty = 2)
lines(constF$biomass.total.obs[1,200:300], col = 'grey', lty = 1)
legend(0, 1600000,  legend = c("No fishing", "F = 0.25 Fmsy", "F = 0.5 Fmsy"),
       cex=0.75, bty = 'n', lty = c(1,1,2), col = c('black','grey','grey50'),y.intersp=2)
plot(C1$biomass.total.true[1,200:300], type = 'l', col = 'red', xlab= "", ylab = "", ylim = c(0,1500000))
lines(C2$biomass.total.true[1,200:300], col = 'red', lty = 2) # low Blim
lines(C3$biomass.total.true[1,200:300], col = 'red4') # high Fmax

legend(0, 1600000,  legend = c("Hockey Stick", "Hockey stick low Blim", "Hockey Stick high F max"),
       col = c('red', 'red','red4'), cex=0.75, bty = 'n', lty = c(1,2,1), y.intersp=2)
par(xpd = NULL)
mtext("Time", side = 1, line = 1, outer = TRUE)
mtext("Sardine biomass", side = 2, line = 1, outer = TRUE)

plot(nofish$biomass.oneplus.true[1,200:300], type = 'l', xlab= "", ylab = "", ylim = c(4000,35000))
lines(fish[1,200:300,1]*props_low[1,200:300,1], col = "red")
lines(fish[1,200:300,1]*props[1,200:300,1], col = "red", lty = 2)
legend(0, 40000,  legend = c("No fishing", "Low variance seabird", "High variance seabird"),
       cex=0.5, bty = 'n', lty = c(1,1,2), col = c('black','red','red'),y.intersp=2)


mtext("Time", side = 1, line = 3)
mtext("Anchovy biomass", side = 2, line = 3)

sardine = nofish
anchovy = nofish
par(mfrow = c(1,2))
par(mar = c(3,2,1,1))
par(oma = c(1,4,1,1))
plot(log(sardine$avgbiomass[1,200:300]), type = 'l', xlab= "", ylab = "")
text(x = 32, y = 13.9, labels = "(A) Sardine")
plot(log(anchovy$avgbiomass[1,200:300]), lty = 2, type = 'l', xlab= "", ylab = "")
text(x = 32, y = 10, labels = "(B) Anchovy")
mtext("Time", side =1, outer = TRUE, line = -1)
mtext("Log Anuual Average \nForage Biomass", side =2, outer = TRUE, line = 1)


par(mfrow = c(3,2))
par(mar = c(3,3,1,1))

performance.measures <- c("LTmeancatch","LTnonzeromeancatch","SDcatch",
                          "n.5yrclose","n.10yrclose","nyrs0catch","meanbiomass",
                          "good4preds","SDbiomass","very.bad4preds","meanDepl",
                          "overallMaxCollapseLength","overallMaxBonanzaLength",
                          "BonanzaLength","CollapseLength","Prob.Collapse",
                          "Collapse.Severity","CV.Catch","Sustained.collapse")

# NEED TO REMOVE FIRST 100 yrs or so???? 
library(matrixStats)
results <- as.data.frame(matrix(0, nrow = 6, ncol = 13))
colnames(results) = c("hcr", "bird1", "bird2", "bird3", "bird4", "bird5","bird6",
                      "meancatch", "sdcatch", "meanbiomass", "sdbiomass" ,"Probcollapse", "years0catch")
HCR = c("nofish","constF","C1","constFdelay", "C1delay", "constFhi") # Took out trend because it was unrealistic-- but using trend in CPUE as adjustment (data-poor method) might be a good idea!
results$hcr = HCR

results$bird1 = c(1, 0.54, 0.8, 0.28, 0.52, 0)
results$bird2 = c(1, 0.95, 0.98, 0.93, 0.96, 0.77)
results$bird3 = c(1, 0.81, 0.92, 0.69, 0.81, 0.03)
results$bird4 = c(1, 0.85, 0.94, 0.77, 0.86, 0.19)
results$bird5 = c(1, 0.9, 0.96, 0.85, 0.91, 0.49)
results$bird6 = c(1, 0.9, 0.95, 0.85, 0.90, 0.48)
results$meancatch = c(mean(rowMeans(nofish$total.catch[,51:1001])), mean(rowMeans(constF$total.catch[,51:1001])),
                mean(rowMeans(C1$total.catch[,51:1001])), mean(rowMeans(constFdelay$total.catch[,51:1001])),
                mean(rowMeans(C1delay$total.catch[,51:1001])), mean(rowMeans(constFhi$total.catch[,51:1001])))
results$sdcatch = c(mean(rowSds(nofish$total.catch[,51:1001])), mean(rowSds(constF$total.catch[,51:1001])),
                mean(rowSds(C1$total.catch[,51:1001])), mean(rowSds(constFdelay$total.catch[,51:1001])),
                mean(rowSds(C1delay$total.catch[,51:1001])), mean(rowSds(constFhi$total.catch[,51:1001])))
results$meanbiomass = c(mean(rowMeans(nofish$biomass.total.true[,51:1001])), mean(rowMeans(constF$biomass.total.true[,51:1001])),
                mean(rowMeans(C1$biomass.total.true[,51:1001])), mean(rowMeans(constFdelay$biomass.total.true[,51:1001])),
                mean(rowMeans(C1delay$biomass.total.true[,51:1001])), mean(rowMeans(constFhi$biomass.total.true[,51:1001])))
results$sdbiomass = c(mean(rowSds(nofish$biomass.total.true[,51:1001])), mean(rowSds(constF$biomass.total.true[,51:1001])),
                      mean(rowSds(C1$biomass.total.true[,51:1001])), mean(rowSds(constFdelay$biomass.total.true[,51:1001])),
                      mean(rowSds(C1delay$biomass.total.true[,51:1001])), mean(rowSds(constFhi$biomass.total.true[,51:1001])))
#thres = 0.2*mean(rowMeans(nofish$biomass.total.true[,51:1001]))
thres = 0.2*equilib$B0
for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(nofish$biomass.total.true[i,51:1001] < thres)) 
}
results[1,12] = mean(temp)/1000

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(constF$biomass.total.true[i,51:1001] < thres)) 
}
results[2,12] = mean(temp)/1000

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(C1$biomass.total.true[i,51:1001] < thres)) 
}
results[3,12] = mean(temp)/1000

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(constFdelay$biomass.total.true[i,51:1001] < thres)) 
}
results[4,12] = mean(temp)/1000

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(C1delay$biomass.total.true[i,51:1001] < thres)) 
}
results[5,12] = mean(temp)/1000

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(constFhi$biomass.total.true[i,51:1001] < thres)) 
}
results[6,12] = mean(temp)/1000

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(C1$total.catch[i,51:1001] == 0)) 
}
results[3,13] = mean(temp)

for(i in 1:nsims) {
  temp = rep(0, length = nsims)
  temp[i] = length(which(C1delay$total.catch[i,51:1001] == 0)) 
}
results[5,13] = mean(temp)


resultslow = results[-1,]
resultslow = resultslow[-5,]

resultslow$bird1 = rescale(resultslow$bird1, to = c(0,3))
resultslow$bird2 = rescale(resultslow$bird2, to = c(0,3))
resultslow$bird3 = rescale(resultslow$bird3, to = c(0,3))
resultslow$bird4 = rescale(resultslow$bird4, to = c(0,3))
resultslow$bird5 = rescale(resultslow$bird5, to = c(0,3))
resultslow$bird6 = rescale(resultslow$bird6, to = c(0,3))
resultslow$meancatch = rescale(resultslow$meancatch, to = c(0,3))
resultslow$sdcatch = rescale(resultslow$sdcatch, to = c(3,0))
resultslow$meanbiomass = rescale(resultslow$meanbiomass, to = c(0,3))
resultslow$sdbiomass = rescale(resultslow$sdbiomass, to = c(3,0))
resultslow$Probcollapse = rescale(resultslow$Probcollapse, to = c(3,0))
resultslow$years0catch = rescale(resultslow$years0catch, to = c(3,0))

library(fmsb)
resultslow = rbind(rep(3, 13), rep(0,13), resultslow)
nofishtest = as.data.frame(results[1:3,2:12])
radarchart(nofishtest)

hcrnew = c("Constant F low", "Hockey stick low", "Constant F low, delayed detection",
           "Hockey stick low, delayed detection")
par(mfrow = c(2,2))
par(mar = c(1,1,1,1))
for(i in 1:4) {
testing = as.data.frame(rbind(resultslow[1:2,2:13], resultslow[i+2,2:13]))

radarchart( testing  , axistype=1 , 
            #custom polygon
            pcol=rgb(0.5,0.5,0.5,0.9) , pfcol=rgb(0.5,0.5,0.5,0.5) , plwd=4 , 
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,3,1), cglwd=0.8,
            vlabels = c("Bird 1", "Bird 2", "Bird 3", "Bird 4", "Bird 5", "Bird 6",
                            "Mean catch", "Minimize SD catch", "Mean prey bio", "Min. SD prey bio", "Min. p(collapse)", "Min. years 0 catch"),
            
            title = hcrnew[i],
            #custom labels
            vlcex=0.7, seg = 3 
)
}
ggradar(results[,7:11],font.radar = "Helvetica",
        grid.label.size=3,axis.label.size=7, 
        legend.text.size = 4,
        axis.labels = axis.labels,
        plot.legend=legend.presence,palette.vec = hcr.colors,
        manual.levels = levels(final.tab$group),
        axis.label.offset=1.1)

X = c(mean(rowSds(nofish$biomass.total.true[,201:1001])), mean(rowSds(constF$biomass.total.true[,201:1001])),
                       mean(rowSds(constFhi$biomass.total.true[,201:1001])))
