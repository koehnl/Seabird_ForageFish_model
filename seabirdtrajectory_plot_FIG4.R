# TO MAKE FIGURE 4 FROM KOEHN ET AL. 2021
# Need to run seabird model with anchovy prey (Run_seabird_3types.R); 
#Run this code with anchovy prey, save seabird-anchovy runs
# Then run base seabird model (Run_seabird_3types.R) with sardine prey; 
# run this with sardine, save seabird-sardine runs

#fish = constF$penguinfood
fish = constF$avgbiomass
# constant low, C1, constant Hi, C2, C3


# 1
# FIND STABLE DISTRIBUTION
maxage_use = maxage[1]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[1], 
                                          K = K_init, 
                                          maxage = maxage[1], egg_s = eggsur, chick_s = chicksur, clutch = clutch[1])

total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,450000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total

sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdoutB = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                   recruits = matrix(nrow = nsims, ncol = Nyear),
                   totalpop = matrix(nrow = nsims, ncol = Nyear),
                   fledglings = matrix(nrow = nsims, ncol = Nyear),
                   eggs = matrix(nrow = nsims, ncol = Nyear),
                   testing = matrix(nrow = nsims, ncol = Nyear),
                   testing2 = matrix(nrow = nsims, ncol = Nyear),
                   testing3 = matrix(nrow = nsims, ncol = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0initB = B0init
#B02initB = B02init
#P0noninitB = P0noninit

set.seed(216)
for(ss in 1:nsims) {
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  #   firsthalfdelay = E[c(TRUE, FALSE)]
  #   secondhalfdelay = E[c(FALSE, TRUE)]
  #   
  #   fishinghalf1 = ((firsthalf*fishing2[ss])-(firsthalfdelay*fishing[ss]))
  #   fishinghalf2 = ((secondhalf*fishing2[ss])-(secondhalfdelay*fishing[ss]))
  #make fishing2 a vector of 1 if no true biomass fishing and only use fishing vector (delayed detection)
  #make fishing vector a vector of 0 if no delayed detection 
  #   
  #   ffbiomass[1,,ss] = fishinghalf1*props[1,,ss]
  #   ffbiomass[2,,ss] = fishinghalf2*props[2,,ss]
  #   ffbiomass_non[1,,ss] = fishinghalf1*nonbreedprops[1,,ss]
  #   ffbiomass_non[2,,ss] = fishinghalf2*nonbreedprops[2,,ss]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  
  #set.seed(821)
  # specialist
  seabirdB = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                             K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  # seabirdB = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
  #                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
  #                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
  #                                            clutch = clutch[1],P0 = B0initB[ss], B0 = B0initB[ss], P02 = B02initB[ss], P0non = P0noninitB[ss],
  #                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
  #                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
  #                                            param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
  #                                            param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
  #                                            param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutB[["breeders"]][ss,] = colSums(seabirdB$pop[2,5:31,])
  seabirdoutB[["recruits"]][ss,] = seabirdB$recruits
  seabirdoutB[["totalpop"]][ss,] = colSums(seabirdB$pop[2,1:31,])
  seabirdoutB[["fledglings"]][ss,] = seabirdB$pop[2,1,]
  seabirdoutB[["eggs"]][ss,] = seabirdB$pop[1,1,]
  seabirdoutB[["testing"]][ss,] = seabirdB$testing
  seabirdoutB[["testing2"]][ss,] = seabirdB$testing2
  seabirdoutB[["testing3"]][ss,] = seabirdB$testing3
  
}


########################## 2 ###################################
maxage_use = maxage[2]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])
stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, 
                                          agebreeding = agebreeding[2], 
                                          K = K_init, 
                                          maxage = maxage[2], egg_s = eggsur, chick_s = chicksur, clutch = clutch[2])
total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,350000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdoutB2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = matrix(nrow = nsims, ncol = Nyear),
                    testing2 = matrix(nrow = nsims, ncol = Nyear),
                    testing3 = matrix(nrow = nsims, ncol = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0initB2 = B0init2
#B02initB2 = B02init2
#P0noninitB2 = P0noninit2

set.seed(216)
for(ss in 1:nsims) {
  
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props_low[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props_low[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops_low[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops_low[2,,ss]
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  # generalist
  seabirdB2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean,  
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutB2[["breeders"]][ss,] = colSums(seabirdB2$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdoutB2[["recruits"]][ss,] = seabirdB2$recruits
  seabirdoutB2[["totalpop"]][ss,] = colSums(seabirdB2$pop[2,1:(maxage[2]+1),])
  seabirdoutB2[["fledglings"]][ss,] = seabirdB2$pop[2,1,]
  seabirdoutB2[["eggs"]][ss,] = seabirdB2$pop[1,1,]
  seabirdoutB2[["testing"]][ss,] = seabirdB2$testing
  seabirdoutB2[["testing2"]][ss,] = seabirdB2$testing2
  seabirdoutB2[["testing3"]][ss,] = seabirdB2$testing3
}

#fish = constFhi$penguinfood #turn off scenario 1 polygon
fish = constFhi$avgbiomass

# 1
# FIND STABLE DISTRIBUTION
maxage_use = maxage[1]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[1], 
                                          K = K_init, 
                                          maxage = maxage[1], egg_s = eggsur, chick_s = chicksur, clutch = clutch[1])

total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,450000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total

sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdoutD = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                   recruits = matrix(nrow = nsims, ncol = Nyear),
                   totalpop = matrix(nrow = nsims, ncol = Nyear),
                   fledglings = matrix(nrow = nsims, ncol = Nyear),
                   eggs = matrix(nrow = nsims, ncol = Nyear),
                   testing = matrix(nrow = nsims, ncol = Nyear),
                   testing2 = matrix(nrow = nsims, ncol = Nyear),
                   testing3 = matrix(nrow = nsims, ncol = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0initB = B0init
#B02initB = B02init
#P0noninitB = P0noninit

set.seed(216)
for(ss in 1:nsims) {
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  #   firsthalfdelay = E[c(TRUE, FALSE)]
  #   secondhalfdelay = E[c(FALSE, TRUE)]
  #   
  #   fishinghalf1 = ((firsthalf*fishing2[ss])-(firsthalfdelay*fishing[ss]))
  #   fishinghalf2 = ((secondhalf*fishing2[ss])-(secondhalfdelay*fishing[ss]))
  #make fishing2 a vector of 1 if no true biomass fishing and only use fishing vector (delayed detection)
  #make fishing vector a vector of 0 if no delayed detection 
  #   
  #   ffbiomass[1,,ss] = fishinghalf1*props[1,,ss]
  #   ffbiomass[2,,ss] = fishinghalf2*props[2,,ss]
  #   ffbiomass_non[1,,ss] = fishinghalf1*nonbreedprops[1,,ss]
  #   ffbiomass_non[2,,ss] = fishinghalf2*nonbreedprops[2,,ss]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  
  #set.seed(821)
  # specialist
  seabirdD = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                             K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  # seabirdB = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
  #                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
  #                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
  #                                            clutch = clutch[1],P0 = B0initB[ss], B0 = B0initB[ss], P02 = B02initB[ss], P0non = P0noninitB[ss],
  #                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
  #                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
  #                                            param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
  #                                            param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
  #                                            param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutD[["breeders"]][ss,] = colSums(seabirdD$pop[2,5:31,])
  seabirdoutD[["recruits"]][ss,] = seabirdD$recruits
  seabirdoutD[["totalpop"]][ss,] = colSums(seabirdD$pop[2,1:31,])
  seabirdoutD[["fledglings"]][ss,] = seabirdD$pop[2,1,]
  seabirdoutD[["eggs"]][ss,] = seabirdD$pop[1,1,]
  seabirdoutD[["testing"]][ss,] = seabirdD$testing
  seabirdoutD[["testing2"]][ss,] = seabirdD$testing2
  seabirdoutD[["testing3"]][ss,] = seabirdD$testing3
  
}


########################## 2 ###################################
maxage_use = maxage[2]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])
stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, 
                                          agebreeding = agebreeding[2], 
                                          K = K_init, 
                                          maxage = maxage[2], egg_s = eggsur, chick_s = chicksur, clutch = clutch[2])
total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,350000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdoutD2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = matrix(nrow = nsims, ncol = Nyear),
                    testing2 = matrix(nrow = nsims, ncol = Nyear),
                    testing3 = matrix(nrow = nsims, ncol = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0initB2 = B0init2
#B02initB2 = B02init2
#P0noninitB2 = P0noninit2

set.seed(216)
for(ss in 1:nsims) {
  
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  ffbiomass[1,,ss] = firsthalf*props_low[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props_low[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops_low[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops_low[2,,ss]
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  # generalist
  seabirdD2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean,  
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutD2[["breeders"]][ss,] = colSums(seabirdD2$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdoutD2[["recruits"]][ss,] = seabirdD2$recruits
  seabirdoutD2[["totalpop"]][ss,] = colSums(seabirdD2$pop[2,1:(maxage[2]+1),])
  seabirdoutD2[["fledglings"]][ss,] = seabirdD2$pop[2,1,]
  seabirdoutD2[["eggs"]][ss,] = seabirdD2$pop[1,1,]
  seabirdoutD2[["testing"]][ss,] = seabirdD2$testing
  seabirdoutD2[["testing2"]][ss,] = seabirdD2$testing2
  seabirdoutD2[["testing3"]][ss,] = seabirdD2$testing3
}

seabirdout$totalpop[is.na(seabirdout$totalpop)] = 0
seabirdout2$totalpop[is.na(seabirdout2$totalpop)] = 0
seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutB2$totalpop[is.na(seabirdoutB2$totalpop)] = 0
seabirdoutD$totalpop[is.na(seabirdoutD$totalpop)] = 0
seabirdoutD2$totalpop[is.na(seabirdoutD2$totalpop)] = 0

sardine1 = seabirdout
sardine2 = seabirdout2
sardineB1 = seabirdoutB
sardineB2 = seabirdoutB2
sardineD1 = seabirdoutD
sardineD2 = seabirdoutD2

anchovy1 = seabirdout
anchovy2 = seabirdout2
anchovyB1 = seabirdoutB
anchovyB2 = seabirdoutB2
anchovyD1 = seabirdoutD
anchovyD2 = seabirdoutD2

# if wanted to run a menhaden run
# men1 = seabirdout
# men2 = seabirdout2
# menB1 = seabirdoutB
# menB2 = seabirdoutB2
# menD1 = seabirdoutD
# menD2 = seabirdoutD2

library(matrixStats)

#### new 4.2020 Tim graph ####
# anchovy1 = seabirdout
# anchovy2 = seabirdout2
# anchovyB1 = seabirdoutB
# anchovyB2 = seabirdoutB2
# anchovyD1 = seabirdoutD
# anchovyD2 = seabirdoutD2
meanB1A = rowMeans(anchovyB1$totalpop[,201:1000]/anchovy1$totalpop[,201:1000])

meanB2A = rowMeans(anchovyB2$totalpop[,201:1000]/anchovy2$totalpop[,201:1000])

meanD1A = rowMeans(anchovyD1$totalpop[,201:1000]/anchovy1$totalpop[,201:1000])

meanD2A = rowMeans(anchovyD2$totalpop[,201:1000]/anchovy2$totalpop[,201:1000])
anchovy = cbind(meanB2A,meanD2A,meanB1A,  meanD1A)

meanB1S = rowMeans(sardineB1$totalpop[,201:1000]/sardine1$totalpop[,201:1000])

meanB2S = rowMeans(sardineB2$totalpop[,201:1000]/sardine2$totalpop[,201:1000])

meanD1S = rowMeans(sardineD1$totalpop[,201:1000]/sardine1$totalpop[,201:1000])

meanD2S = rowMeans(sardineD2$totalpop[,201:1000]/sardine2$totalpop[,201:1000])

#probability of extinction of <10% un-fished
birdthres1 = median(rowMeans(anchovy1$totalpop[,201:1000]))*0.1
birdthres2 = median(rowMeans(anchovy2$totalpop[,201:1000]))*0.1 # note flexible actually has smaller population
# size even without fishing. makes sense because shorter life span. 
prob1= median((apply(anchovyB1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1))}))/800)
prob2 = median((apply(anchovyB2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2))}))/800)
prob3= median((apply(anchovyD1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1))}))/800)
prob4 = median((apply(anchovyD2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2))}))/800)

quan1 = ((apply(anchovyB1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1))}))/800)
quan2 = ((apply(anchovyB2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2))}))/800)
quan3 = ((apply(anchovyD1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1))}))/800)
quan4 = ((apply(anchovyD2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2))}))/800)


birdthres1S = median(rowMeans(sardine1$totalpop[,201:1000]))*0.1
birdthres2S= median(rowMeans(sardine2$totalpop[,201:1000]))*0.1 # note flexible actually has smaller population
# size even without fishing. makes sense because shorter life span. 
prob1S= median((apply(sardineB1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1S))}))/800)
prob2S = median((apply(sardineB2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2S))}))/800)
prob3S= median((apply(sardineD1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1S))}))/800)
prob4S = median((apply(sardineD2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2S))}))/800)


quan1S = ((apply(sardineB1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1S))}))/800)
quan2S = ((apply(sardineB2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2S))}))/800)
quan3S = ((apply(sardineD1$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres1S))}))/800)
quan4S = ((apply(sardineD2$totalpop[,201:1000], 1, function(x){ length(which(x< birdthres2S))}))/800)

par(mfrow =c (3,2))
par(xpd = NA)
par(mar = c(1,2,2,2)) #
par(oma = c(7,9,2,10)) # bottom left top right
x = c(0.98,1,1.1,1.12)
plot(x, c( NA, median(meanB2A),median(meanD2A), NA), ylim = c(0,1), xlab = "", ylab = "", main  = "Anchovy prey",
     axes = FALSE, cex = 2)
points(x, c(NA, median(meanB1A), median(meanD1A), NA), pch = 16, cex = 2)
segments(x[2], median(meanB2A), x[3], median(meanD2A) )
segments(x[2], median(meanB1A), x[3], median(meanD1A) )
text(x = 0.88, y = 0.5, "Median relative \n abundance", cex = 1.25,xpd = NA)
#Axis(side=1, at = c(1,1.1), labels=FALSE)
Axis(side=2, labels=TRUE, las = 1)


plot(x, c( NA, median(meanB2S),median(meanD2S), NA), ylim = c(0,1), xlab = "", ylab = "", main  = "Sardine prey",
     axes = FALSE, cex = 2)
points(x, c(NA, median(meanB1S), median(meanD1S), NA), pch = 16, cex = 2)
segments(x[2], median(meanB2S), x[3], median(meanD2S) )
segments(x[2], median(meanB1S), x[3], median(meanD1S) )
#text(x = 0.88, y = 0.5, "Median relative \n abundance", cex = 1.25,xpd = NA)
#Axis(side=1, at = c(1,1.1), labels=FALSE)
Axis(side=2, labels=TRUE, las = 1)
#legend(0.85,2, legend = c("Flexible seabird", "Restricted Seabird"), col = c("black", "black"), pch = c(1,16),
     #  ncol = 2, cex = 1.2, pt.cex = 2,xpd = NA)


 plot(x, c(NA, prob2,prob4, NA), ylim = c(0,1), xlab = "", ylab = "",
      axes = FALSE, cex = 2)
 points(x,c(NA, prob1, prob3,NA), pch = 16, cex = 2)
 segments(x[2], prob2, x[3], prob4 )
 segments(x[2], prob1, x[3], prob3 )
 #Axis(side=1, at = c(1,1.1), labels=FALSE)
 Axis(side=2, labels=TRUE, las = 1)
 text(x = 0.88, y = 0.5, "Prob(Extinct)", cex = 1.25,xpd = NA)
 
 legend(1.34,2.3, legend = c("Flexible seabird", "Restricted seabird"), col = c("black", "black"), pch = c(1,16),
       cex = 1.25, pt.cex = 2,xpd = NA)
 
 plot(x, c(NA, prob2S,prob4S, NA), ylim = c(0,1), xlab = "", ylab = "",
      axes = FALSE, cex = 2)
 points(x,c(NA, prob1S, prob3S,NA), pch = 16, cex = 2)
 segments(x[2], prob2S, x[3], prob4S )
 segments(x[2], prob1S, x[3], prob3S )
 #Axis(side=1, at = c(1,1.1), labels=FALSE)
 Axis(side=2, labels=TRUE, las = 1)
# text(x = 0.88, y = 0.5, "Prob(Extinct)", cex = 1.25,xpd = NA)

 plot(x, c(NA, var(meanB2A),var(meanD2A),NA), ylim = c(0,0.055), xlab = "", ylab = "",
      axes = FALSE, cex = 2)
 points(x, c(NA, var(meanB1A), var(meanD1A), NA), pch = 16, cex = 2)
 segments(x[2], var(meanB2A), x[3], var(meanD2A) )
 segments(x[2], var(meanB1A), x[3], var(meanD1A) )
 #segments(x[1], 0, x[2], 0)
 #Axis(side=1, labels=FALSE, tick = FALSE, line = TRUE)
 Axis(side=2, labels=TRUE, las = 1)
 text(x = 0.88, y = 0.025, "Variance", cex = 1.25,xpd = NA)
 
 axis(side = 1, labels = c(expression("0.25 F"["msy"]), expression("0.5 F"["msy"])),
      at = c(1.0,1.1), line = 1.5, cex.axis = 1.5, tick = FALSE)
 #mtext("Constant Fishing Rate", side =1, line = 4)

 
 plot(x, c(NA, var(meanB2S),var(meanD2S),NA), ylim = c(0,0.055), xlab = "", ylab = "",
      axes = FALSE, cex = 2)
 points(x, c(NA, var(meanB1S), var(meanD1S), NA), pch = 16, cex = 2)
 segments(x[2], var(meanB2S), x[3], var(meanD2S) )
 segments(x[2], var(meanB1S), x[3], var(meanD1S) )
 #segments(x[1], 0, x[2], 0)
 #Axis(side=1, labels=FALSE, tick = FALSE, line = TRUE)
 Axis(side=2, labels=TRUE, las = 1)
 #text(x = 0.88, y = 0.01, "Variance", cex = 1.25,xpd = NA)
 
 axis(side = 1, labels = c(expression("0.25 F"["msy"]), expression("0.5 F"["msy"])),
      at = c(1.0,1.1), line = 1.5, cex.axis = 1.5, tick = FALSE)
 text(0.94, -0.05, "Constant Fishing Rate", xpd = NA, cex = 1.5)

 
 ## OLD TEST CODE NOT USED ##### 
 #### confidence intervals #####
 
 par(mfrow =c (2,2))
 par(xpd = NA)
 par(mar = c(1,2,2,2)) #
 par(oma = c(7,9,2,10)) # bottom left top right
 x = c(0.98,1,1.1,1.12)
 plot(x, c( NA, median(meanB2A),median(meanD2A), NA), ylim = c(0,1), xlab = "", ylab = "", main  = "Anchovy prey",
      axes = FALSE, cex = 2)
 points(x, c(NA, median(meanB1A), median(meanD1A), NA), pch = 16, cex = 2)

 arrows(x, c(NA, quantile(meanB1A, 0.75), quantile(meanD1A, 0.75),NA),
        x, c(NA, quantile(meanB1A, 0.25), quantile(meanD1A, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 arrows(x, c(NA, quantile(meanB2A, 0.75), quantile(meanD2A, 0.75),NA),
        x, c(NA, quantile(meanB2A, 0.25), quantile(meanD2A, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 segments(x[2], median(meanB2A), x[3], median(meanD2A) )
 segments(x[2], median(meanB1A), x[3], median(meanD1A) )
 text(x = 0.88, y = 0.5, "Median relative \n abundance", cex = 1.25,xpd = NA)
 #Axis(side=1, at = c(1,1.1), labels=FALSE)
 Axis(side=2, labels=TRUE, las = 1)
 
 
 plot(x, c( NA, median(meanB2S),median(meanD2S), NA), ylim = c(0,1), xlab = "", ylab = "", main  = "Sardine prey",
      axes = FALSE, cex = 2)
 points(x, c(NA, median(meanB1S), median(meanD1S), NA), pch = 16, cex = 2)
 
 arrows(x, c(NA, quantile(meanB1S, 0.75), quantile(meanD1S, 0.75),NA),
       x, c(NA, quantile(meanB1S, 0.25), quantile(meanD1S, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 arrows(x, c(NA, quantile(meanB2S, 0.75), quantile(meanD2S, 0.75),NA),
        x, c(NA, quantile(meanB2S, 0.25), quantile(meanD2S, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 
 segments(x[2], median(meanB2S), x[3], median(meanD2S) )
 segments(x[2], median(meanB1S), x[3], median(meanD1S) )
 #text(x = 0.88, y = 0.5, "Median relative \n abundance", cex = 1.25,xpd = NA)
 #Axis(side=1, at = c(1,1.1), labels=FALSE)
 Axis(side=2, labels=TRUE, las = 1)
 #legend(0.85,2, legend = c("Flexible seabird", "Restricted Seabird"), col = c("black", "black"), pch = c(1,16),
 #  ncol = 2, cex = 1.2, pt.cex = 2,xpd = NA)
 
 
 plot(x, c(NA, prob2,prob4, NA), ylim = c(0,1), xlab = "", ylab = "",
      axes = FALSE, cex = 2)
 points(x,c(NA, prob1, prob3,NA), pch = 16, cex = 2)
 arrows(x, c(NA, quantile(quan2, 0.75), quantile(quan4, 0.75),NA),
        x, c(NA, quantile(quan2, 0.25), quantile(quan4, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 arrows(x, c(NA, quantile(quan1, 0.75), quantile(quan3, 0.75),NA),
        x, c(NA, quantile(quan1, 0.25), quantile(quan3, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 segments(x[2], prob2, x[3], prob4 )
 segments(x[2], prob1, x[3], prob3 )
 #Axis(side=1, at = c(1,1.1), labels=FALSE)
 Axis(side=2, labels=TRUE, las = 1)
 text(x = 0.88, y = 0.5, "Prob(Extinct)", cex = 1.25,xpd = NA)
 
 legend(1.34,2.3, legend = c("Flexible seabird", "Restricted seabird"), col = c("black", "black"), pch = c(1,16),
        cex = 1.25, pt.cex = 2,xpd = NA)
 
 plot(x, c(NA, prob2S,prob4S, NA), ylim = c(0,1), xlab = "", ylab = "",
      axes = FALSE, cex = 2)
 points(x,c(NA, prob1S, prob3S,NA), pch = 16, cex = 2)
 arrows(x, c(NA, quantile(quan2S, 0.75), quantile(quan4S, 0.75),NA),
        x, c(NA, quantile(quan2S, 0.25), quantile(quan4S, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 arrows(x, c(NA, quantile(quan1S, 0.75), quantile(quan3S, 0.75),NA),
        x, c(NA, quantile(quan1S, 0.25), quantile(quan3S, 0.25),NA),length=0.1, angle=90, code=3, lwd = 1)
 segments(x[2], prob2S, x[3], prob4S )
 segments(x[2], prob1S, x[3], prob3S )
 #Axis(side=1, at = c(1,1.1), labels=FALSE)
 Axis(side=2, labels=TRUE, las = 1)
 # text(x = 0.88, y = 0.5, "Prob(Extinct)", cex = 1.25,xpd = NA)
 
 axis(side = 1, labels = c(expression("0.25 F"["msy"]), expression("0.5 F"["msy"])),
      at = c(1.0,1.1), line = 1.5, cex.axis = 1.5, tick = FALSE)
 text(0.94, -0.2, "Constant Fishing Rate", xpd = NA, cex = 1.5)
 
 
##### metrics for last 100 years 
#mean1S = rowMeans(sardine1$totalpop[,900:1000])

#mean2S = rowMeans(sardine2$totalpop[,900:1000])

meanB1S = rowMeans(sardineB1$totalpop[,900:1000]/sardine1$totalpop[,900:1000])

meanB2S = rowMeans(sardineB2$totalpop[,900:1000]/sardine2$totalpop[,900:1000])

meanD1S = rowMeans(sardineD1$totalpop[,900:1000]/sardine1$totalpop[,900:1000])

meanD2S = rowMeans(sardineD2$totalpop[,900:1000]/sardine2$totalpop[,900:1000])

#mean1A = rowMeans(anchovy1$totalpop[,900:1000])
#mean2A = rowMeans(anchovy2$totalpop[,900:1000])

meanB1A = rowMeans(anchovyB1$totalpop[,900:1000]/anchovy1$totalpop[,900:1000])

meanB2A = rowMeans(anchovyB2$totalpop[,900:1000]/anchovy2$totalpop[,900:1000])

meanD1A = rowMeans(anchovyD1$totalpop[,900:1000]/anchovy1$totalpop[,900:1000])

meanD2A = rowMeans(anchovyD2$totalpop[,900:1000]/anchovy2$totalpop[,900:1000])

#mean1M = rowMeans(men1$totalpop[,900:1000])
#mean2M = rowMeans(men2$totalpop[,900:1000])

meanB1M = rowMeans(menB1$totalpop[,900:1000]/men1$totalpop[,900:1000])

meanB2M = rowMeans(menB2$totalpop[,900:1000]/men2$totalpop[,900:1000])

meanD1M = rowMeans(menD1$totalpop[,900:1000]/men1$totalpop[,900:1000])

meanD2M = rowMeans(menD2$totalpop[,900:1000]/men2$totalpop[,900:1000])

sardine = cbind(meanB2S,meanD2S,meanB1S,  meanD1S)
anchovy = cbind(meanB2A,meanD2A,meanB1A,  meanD1A)
menhaden = cbind(meanB2M,meanD2M,meanB1M,  meanD1M)

sardine = as.data.frame(sardine)
names(sardine) = c("Low F", "High F", "Low F", "High F")
anchovy = as.data.frame(anchovy)
names(anchovy) = c("Low F", "High F", "Low F", "High F")
menhaden = as.data.frame(menhaden)
names(menhaden) = c("Low F", "High F", "Low F", "High F")
boxplot(cbind(sardine, anchovy, menhaden))
par(mfrow = c(2,1))
#par(oma = c(2,4,7,1))
par(mar = c(2,2,2,2))
par(oma = c(3,3,4,1))
par(xpd = FALSE)
colors = gray.colors(n = 5)
boxplot(sardine, frame = F, axes = F, col = colors[4:1], border = colors[4:1])#c("darkgrey", "black","darkgrey", "black"),
       #xlab = "Sardine prey")
axis(side = 2, line = -1)
mtext("Prey = Sardine", side = 3, line = 1)
abline(v = 2.5, lty = 2)
boxplot(anchovy, frame = F, axes = F, col = colors[4:1], border = colors[4:1])#c("darkgrey", "black","darkgrey", "black"),
        #xlab= "Anchovy prey")
axis(side = 2, line = -1)
abline(v = 2.5, lty = 2, xpd = FALSE)
mtext(side = 2, "Mean seabird abundance relative to unfished", xpd = NA, line = -0.5, outer = TRUE)
mtext("Prey = Anchovy", side = 3, line = 1)
mtext(side = 1, line = 3, "Constant fishing rate")
# boxplot(menhaden, xlab = "", axes = F, col = colors[4:1], border = colors[4:1],#c("darkgrey", "black","darkgrey", "black"),
#        main = "", frame = F, cex.axis = 1.5)
# mtext("Prey = Menhaden", side = 3, line =1)
# mtext("Constant Fishing Rate", side = 1, line =5)

# legend(1,-1, legend = c("Generalist, Low Const. F", "Generalist, High Const. F",
#                             "Specialist, Low Const. F", "Specialist, High Const. F"),
#        lty = 1, col = colors[4:1], ncol = 2,
#        xpd = NA, cex = 1)
abline(v = 2.5, lty = 2)
axis(side = 2, line = -1)
axis(side = 1, labels = c(expression("0.25 F"["msy"]), expression("0.5 F"["msy"]),
                          expression("0.25 F"["msy"]), expression("0.5 F"["msy"])), 
     at = c(1,2,3,4), line = 1, cex.axis = 1.5)
#axis(side = 1, tck = 0, labels = FALSE)
library(pBrackets)
brackets(0.5,3.3,2.4,3.3, h = 0.3,  ticks = 0.5, curvature = 0.1, type = 2
         ,col = 1, lwd = 1, lty = 1, xpd = NA)
brackets(2.6,3.3,4.5,3.3, h = 0.3,  ticks = 0.5, curvature = 0.1, type = 1,col = 1, lwd = 1, lty = 1, xpd = NA)
text(x = 1.5,y= 3.75,label = "Flexible Seabird", xpd = NA, cex = 1.5)
text(x = 3.5,y= 3.75,label = "Restricted Seabird", xpd = NA, cex = 1.5)

test = boxplot(sardine)
test2 = boxplot(anchovy)
test3 = boxplot(menhaden)
sardine2 = as.data.frame(sardine)
names(sardine2) = c("Generalist, low const. F", "Generalist, high const. F",
                    "Specialist, low const. F", "Specialist, high const. F")
####################3
mean1 = colMedians(sardine1$totalpop[,201:1000])
CI1 = apply(sardine1$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean2 = colMedians(sardine2$totalpop[,201:1000])
CI2 = apply(sardine2$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanB1 = colMedians(sardineB1$totalpop[,201:1000])
CIB1 = apply(sardineB1$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanB2 = colMedians(sardineB2$totalpop[,201:1000])
CIB2 = apply(sardineB2$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanD1 = colMedians(sardineD1$totalpop[,201:1000])
CID1 = apply(sardineD1$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanD2 = colMedians(sardineD2$totalpop[,201:1000])
CID2 = apply(sardineD2$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean1A = colMedians(sardine1$totalpop[,201:1000])
CI1A= apply(sardine1$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean2A = colMedians(sardine2$totalpop[,201:1000])
CI2A = apply(sardine2$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanB1A = colMedians(sardineB1$totalpop[,201:1000])
CIB1A = apply(sardineB1$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanB2A = colMedians(sardineB2$totalpop[,201:1000])
CIB2A = apply(sardineB2$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanD1A = colMedians(sardineD1$totalpop[,201:1000])
CID1A = apply(sardineD1$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

meanD2A = colMedians(sardineD2$totalpop[,201:1000])
CID2A = apply(sardineD2$totalpop[,201:1000], 2, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

par(mar = c(6,4,1,8))
x.seq = 1:800
plot(x.seq, meanB1/mean1, type = 'n', ylim = c(0,1), bty = "n", xlab = "", ylab = "")
polygon(c(x.seq, rev(x.seq)), c(CIB1[1,]/CI1[1,], rev(CIB1[2,]/CI1[2,])), col = rgb(0,0,0,0.3), border = NA)
lines(x.seq, meanB1/mean1, lty = 2)

polygon(c(x.seq, rev(x.seq)), c(CIB2[1,]/CI2[1,], rev(CIB2[2,]/CI2[2,])), col = rgb(0,0,0,0.3), border = NA)
lines(x.seq, meanB2/mean2, lty = 1)

polygon(c(x.seq, rev(x.seq)), c(CID1[1,]/CI1[1,], rev(CID1[2,]/CI1[2,])), col = rgb(0,0,0,0.3), border = NA)
lines(x.seq, meanD1/mean1, lty = 2)

polygon(c(x.seq, rev(x.seq)), c(CID2[1,]/CI2[1,], rev(CID2[2,]/CI2[2,])), col = rgb(0,0,0,0.3), border = NA)
lines(x.seq, meanD2/mean2, lty = 1)

#### 2nd fish
rednew = rgb(100,0,0,alpha = 30, maxColorValue = 100)
polygon(c(x.seq, rev(x.seq)), c(CIB1A[1,]/CI1A[1,], rev(CIB1A[2,]/CI1A[2,])), col = rednew, border = NA)
lines(x.seq, meanB1A/mean1A, lty = 2, col = "red")

polygon(c(x.seq, rev(x.seq)), c(CIB2A[1,]/CI2A[1,], rev(CIB2A[2,]/CI2A[2,])), col = rednew, border = NA)
lines(x.seq, meanB2A/mean2A, lty = 1, col = "red")

polygon(c(x.seq, rev(x.seq)), c(CID1A[1,]/CI1A[1,], rev(CID1A[2,]/CI1A[2,])), col = rednew, border = NA)
lines(x.seq, meanD1A/mean1A, lty = 2, col = "red")

polygon(c(x.seq, rev(x.seq)), c(CID2A[1,]/CI2A[1,], rev(CID2A[2,]/CI2A[2,])), col = rednew, border = NA)
lines(x.seq, meanD2A/mean2A, lty = 1, col = "red")
mtext("Population", side = 2, line = 2.5)
mtext("Time", side = 1, line = 2.5)

text(x = 960, y = 0.95, "Generalist, 0.25 Fmsy", xpd = NA, cex = 0.8)
text(x = 950, y = 0.87, "Generalist, 0.5 Fmsy", xpd = NA, cex =0.8)

text(x = 950, y = 0.4, "Specialist, 0.25 Fmsy", xpd = NA, cex = 0.8)
text(x = 950, y = 0.0, "Specialist, 0.5 Fmsy", xpd = NA, cex =0.8)
legend(x = 600, y = -0.2, legend = c("Anchovy prey", "Sardine prey"), col = c("Black", "red"), lty = 1, xpd = NA, bty = 'n')
