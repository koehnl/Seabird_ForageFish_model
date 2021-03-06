par(mfrow = c(2,2))
par(oma =c(2,3,1,2))
par(mar = c(3,3,1,2))

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

################ cut-off rule 1 ##################
#fish = C1$penguinfood
fish = C1$avgbiomass

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

seabirdoutC = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdC = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
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
  seabirdoutC[["breeders"]][ss,] = colSums(seabirdC$pop[2,5:31,])
  seabirdoutC[["recruits"]][ss,] = seabirdC$recruits
  seabirdoutC[["totalpop"]][ss,] = colSums(seabirdC$pop[2,1:31,])
  seabirdoutC[["fledglings"]][ss,] = seabirdC$pop[2,1,]
  seabirdoutC[["eggs"]][ss,] = seabirdC$pop[1,1,]
  seabirdoutC[["testing"]][ss,] = seabirdC$testing
  seabirdoutC[["testing2"]][ss,] = seabirdC$testing2
  seabirdoutC[["testing3"]][ss,] = seabirdC$testing3
  
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

seabirdoutC2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdC2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean,  
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutC2[["breeders"]][ss,] = colSums(seabirdC2$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdoutC2[["recruits"]][ss,] = seabirdC2$recruits
  seabirdoutC2[["totalpop"]][ss,] = colSums(seabirdC2$pop[2,1:(maxage[2]+1),])
  seabirdoutC2[["fledglings"]][ss,] = seabirdC2$pop[2,1,]
  seabirdoutC2[["eggs"]][ss,] = seabirdC2$pop[1,1,]
  seabirdoutC2[["testing"]][ss,] = seabirdC2$testing
  seabirdoutC2[["testing2"]][ss,] = seabirdC2$testing2
  seabirdoutC2[["testing3"]][ss,] = seabirdC2$testing3
}

####################### high constant fishing ###################
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
seabirdoutD$totalpop[is.na(seabirdoutD$totalpop)] = 0


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

seabirdoutD2$totalpop[is.na(seabirdoutD2$totalpop)] = 0

median(rowMeans(seabirdoutD2$totalpop[,900:1000])/rowMeans(seabirdout2$totalpop[,900:1000]))
plot(rowMeans(seabirdoutD2$totalpop[,900:1000])/rowMeans(seabirdout2$totalpop[,900:1000]))
####################### cut-off rule 2 - low Blim ###########
#fish = C2$penguinfood
fish = C2$avgbiomass

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

seabirdoutE = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdE = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
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
  seabirdoutE[["breeders"]][ss,] = colSums(seabirdE$pop[2,5:31,])
  seabirdoutE[["recruits"]][ss,] = seabirdE$recruits
  seabirdoutE[["totalpop"]][ss,] = colSums(seabirdE$pop[2,1:31,])
  seabirdoutE[["fledglings"]][ss,] = seabirdE$pop[2,1,]
  seabirdoutE[["eggs"]][ss,] = seabirdE$pop[1,1,]
  seabirdoutE[["testing"]][ss,] = seabirdE$testing
  seabirdoutE[["testing2"]][ss,] = seabirdE$testing2
  seabirdoutE[["testing3"]][ss,] = seabirdE$testing3
  
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

seabirdoutE2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdE2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean,  
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutE2[["breeders"]][ss,] = colSums(seabirdE2$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdoutE2[["recruits"]][ss,] = seabirdE2$recruits
  seabirdoutE2[["totalpop"]][ss,] = colSums(seabirdE2$pop[2,1:(maxage[2]+1),])
  seabirdoutE2[["fledglings"]][ss,] = seabirdE2$pop[2,1,]
  seabirdoutE2[["eggs"]][ss,] = seabirdE2$pop[1,1,]
  seabirdoutE2[["testing"]][ss,] = seabirdE2$testing
  seabirdoutE2[["testing2"]][ss,] = seabirdE2$testing2
  seabirdoutE2[["testing3"]][ss,] = seabirdE2$testing3
}


##################### cut-off rule 3 - high F max #############
#fish = C3$penguinfood

fish = C3$avgbiomass
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

seabirdoutF = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdF = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
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
  seabirdoutF[["breeders"]][ss,] = colSums(seabirdF$pop[2,5:31,])
  seabirdoutF[["recruits"]][ss,] = seabirdF$recruits
  seabirdoutF[["totalpop"]][ss,] = colSums(seabirdF$pop[2,1:31,])
  seabirdoutF[["fledglings"]][ss,] = seabirdF$pop[2,1,]
  seabirdoutF[["eggs"]][ss,] = seabirdF$pop[1,1,]
  seabirdoutF[["testing"]][ss,] = seabirdF$testing
  seabirdoutF[["testing2"]][ss,] = seabirdF$testing2
  seabirdoutF[["testing3"]][ss,] = seabirdF$testing3
  
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

seabirdoutF2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdF2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean,  
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutF2[["breeders"]][ss,] = colSums(seabirdF2$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdoutF2[["recruits"]][ss,] = seabirdF2$recruits
  seabirdoutF2[["totalpop"]][ss,] = colSums(seabirdF2$pop[2,1:(maxage[2]+1),])
  seabirdoutF2[["fledglings"]][ss,] = seabirdF2$pop[2,1,]
  seabirdoutF2[["eggs"]][ss,] = seabirdF2$pop[1,1,]
  seabirdoutF2[["testing"]][ss,] = seabirdF2$testing
  seabirdoutF2[["testing2"]][ss,] = seabirdF2$testing2
  seabirdoutF2[["testing3"]][ss,] = seabirdF2$testing3
}

##### MSY for generalist #######
fish = msycatch$avgbiomass

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

seabirdoutmsy = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
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
  seabirdmsy = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean,  
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdoutmsy[["breeders"]][ss,] = colSums(seabirdmsy$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdoutmsy[["recruits"]][ss,] = seabirdmsy$recruits
  seabirdoutmsy[["totalpop"]][ss,] = colSums(seabirdmsy$pop[2,1:(maxage[2]+1),])
  seabirdoutmsy[["fledglings"]][ss,] = seabirdmsy$pop[2,1,]
  seabirdoutmsy[["eggs"]][ss,] = seabirdmsy$pop[1,1,]
  seabirdoutmsy[["testing"]][ss,] = seabirdmsy$testing
  seabirdoutmsy[["testing2"]][ss,] = seabirdmsy$testing2
  seabirdoutmsy[["testing3"]][ss,] = seabirdmsy$testing3
}


seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutB2$totalpop[is.na(seabirdoutB2$totalpop)] = 0
seabirdoutC$totalpop[is.na(seabirdoutC$totalpop)] = 0
seabirdoutC2$totalpop[is.na(seabirdoutC2$totalpop)] = 0
seabirdoutD$totalpop[is.na(seabirdoutD$totalpop)] = 0
seabirdoutD2$totalpop[is.na(seabirdoutD2$totalpop)] = 0
seabirdoutE$totalpop[is.na(seabirdoutE$totalpop)] = 0
seabirdoutE2$totalpop[is.na(seabirdoutE2$totalpop)] = 0
seabirdoutF$totalpop[is.na(seabirdoutF$totalpop)] = 0
seabirdoutF2$totalpop[is.na(seabirdoutF2$totalpop)] = 0
seabirdoutmsy$totalpop[is.na(seabirdoutmsy$totalpop)] = 0


library(matrixStats)
mean1 = rowMeans(seabirdout$totalpop[,200:1000])
CI1 = apply(seabirdout$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean2 = rowMeans(seabirdout2$totalpop[,200:1000])
CI2 = apply(seabirdout2$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})


mean1B = rowMeans(seabirdoutB$totalpop[,200:1000])
CI1B = apply(seabirdoutB$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
test = rowMedians(seabirdoutB$totalpop[,200:1000])

mean2B = rowMeans(seabirdoutB2$totalpop[,200:1000])
CI2B = apply(seabirdoutB2$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean1C = rowMeans(seabirdoutC$totalpop[,200:1000])
CI1C = apply(seabirdoutC$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
mean2C = rowMeans(seabirdoutC2$totalpop[,200:1000])
CI2C = apply(seabirdoutC2$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean1D = rowMeans(seabirdoutD$totalpop[,200:1000])
CI1D = apply(seabirdoutD$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
test = rowMedians(seabirdoutD$totalpop[,200:1000])
mean2D = rowMeans(seabirdoutD2$totalpop[,200:1000])
CI2D = apply(seabirdoutD2$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean1E = rowMeans(seabirdoutE$totalpop[,200:1000])
CI1E = apply(seabirdoutE$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
mean2E = rowMeans(seabirdoutE2$totalpop[,200:1000])
CI2E = apply(seabirdoutE2$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean1F = rowMeans(seabirdoutF$totalpop[,200:1000])
CI1F = apply(seabirdoutF$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
mean2F = rowMeans(seabirdoutF2$totalpop[,200:1000])
CI2F = apply(seabirdoutF2$totalpop[,200:1000], 1, function(x){median(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

# constant low = B, C1 = C, constant high = D, C2 = E, C3 = F
# A for anchovy, S for sardine runs
# middleA = c(median(mean1C/mean1), median(mean2C/mean2),
#             median(mean1E/mean1), median(mean2E/mean2),
#             median(mean1F/mean1), median(mean2F/mean2),
#             median(mean1B/mean1), median(mean2B/mean2),
#             median(mean1D/mean1), median(mean2D/mean2))
# # lowerA = c(median(CI1C[1,]/CI1[1,]), median(CI2C[1,]/CI2[1,]),
# #            median(CI1E[1,]/CI1[1,]), median(CI2E[1,]/CI2[1,]),
# #            median(CI1F[1,]/CI1[1,]), median(CI2F[1,]/CI2[1,]),
# #            median(CI1B[1,]/CI1[1,]), median(CI2B[1,]/CI2[1,]),
# #            median(CI1D[1,]/CI1[1,]), median(CI2D[1,]/CI2[1,]))
# # upperA = c(median(CI1C[2,]/CI1[2,]), median(CI2C[2,]/CI2[2,]),
# #            median(CI1E[2,]/CI1[2,]), median(CI2E[2,]/CI2[2,]),
# #            median(CI1F[2,]/CI1[2,]), median(CI2F[2,]/CI2[2,]),
# #            median(CI1B[2,]/CI1[2,]), median(CI2B[2,]/CI2[2,]),
# #            median(CI1D[2,]/CI1[2,]), median(CI2D[2,]/CI2[2,]))
# lowerA =  c(quantile(mean1C/mean1, 0.025), quantile(mean2C/mean2, 0.025),
#             quantile(mean1E/mean1,0.025), quantile(mean2E/mean2, 0.025),
#             quantile(mean1F/mean1, 0.025), quantile(mean2F/mean2, 0.025),
#             quantile(mean1B/mean1, 0.025), quantile(mean2B/mean2, 0.025),
#             quantile(mean1D/mean1, 0.025), quantile(mean2D/mean2, 0.025))
# upperA =  c(quantile(mean1C/mean1, 0.975), quantile(mean2C/mean2, 0.975),
#             quantile(mean1E/mean1,0.975), quantile(mean2E/mean2, 0.975),
#             quantile(mean1F/mean1, 0.975), quantile(mean2F/mean2, 0.975),
#             quantile(mean1B/mean1, 0.975), quantile(mean2B/mean2, 0.975),
#             quantile(mean1D/mean1, 0.975), quantile(mean2D/mean2, 0.975))

middleA = c(median(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])))
lowerA = c(quantile(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025))
upperA = c(quantile(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975))

middleS = c(median(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
            median(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
            median(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])))
lowerS = c(quantile(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
           quantile(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025))
upperS = c(quantile(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
           quantile(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975))

setwd(here::here())
figwd <- file.path(getwd(),"Plots/")

library(RColorBrewer)
library(scales)

# Demo figure for control rules
source(file.path("ForageFishModel/Control Rules/smith_oceana.R"))
source(file.path("ForageFishModel/Control Rules/cfp.R"))
source(file.path("ForageFishModel/Control Rules/hockey-stick.R"))
source(file.path("ForageFishModel/Control Rules/trend-based-rule.R"))



B <- 1:10000
B0 <- 7000
Bmsy <- 3000
Fmsy = 0.7
m = 0.8

C1.vec <- C2.vec <- C3.vec <- cfp.vec <- vector()
for(i in 1:length(B)){
  C1.vec[i] <- calc.F.oceana(Bt = B[i],Blim = 0.4*B0,Btarget = 0.8*B0, M = m)
  C2.vec[i] <- calc.F.oceana(Bt = B[i], Blim = 0.1*B0, Btarget = 0.8*B0, M = m)
  C3.vec[i] <- calc.F.stick(Bt = B[i], Blim = 0.4*B0, Btarget = 0.8*B0, Fmax = Fmsy)
  cfp.vec[i] <- calc.F.stick(Bt = B[i],Blim = 0.5*Bmsy, Btarget = 0.4*B0,Fmax = Fmsy)
}

constf.vec <- rep(0.5*Fmsy,times=length(B)) # was Fmsy before
constf.vec_lo <- rep(0.25*Fmsy,times=length(B)) # was 0.5 Fmsy before 


palette <- brewer.pal(6,"Spectral")
hcr.colors <- palette[c(6,5,4,3,1,2)]



layout(matrix(c(1,2,3,1,2,3,1,2,0), nrow = 3, byrow = TRUE))
#par(mfrow = c(1,3))
par(oma = c(4,6,2,2))
par(mar = c(2,2,2,1))
x = c(1,1,2,2,3,3,4,4,5,5)
names = c("Hockey-stick A","Hockey-stick B", "Hockey-stick C", "Constant Low", "Constant \nmoderate")
plot(middleA, x, frame = F, xlab = "", ylab = "", 
     xlim = c(0,1), axes = F, pch = 19, col = c("Black", "Red"), main = "Anchovy Fishery")
arrows(lowerA, x, 
       upperA, x, length=0, angle=90, code=3, lwd = 1,col = c("Black", "Red"))
#text(middleA, x, labels =round(middleA,2) )

axis(2, at = c(1,2,3,4,5), labels = names, las = 2, pos = -0.1, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = 0.8)
legend(-0.5,0.4,xpd=NA, col = c("Black", "Red"), pch = 19, legend = c("Restricted", "Flexible"))
#text(-0.4,5.4, label = "(A)", xpd = NA)

#par(mar = c(6,8.5,2,2))
#x = c(1,1,2,2,3,3,4,4,5,5)
#names = c("Hockey-stick -\n Fmax 0.5 Fmsy","Hockey-stick - \nLow Blim", "Hockey-stick - \nFmax Fmsy", "Constant Low - \n0.25 Fmsy", "Constant Mod. - \n0.5 Fmsy")
plot(middleS, x, frame = F, xlab = "", ylab = "", 
     xlim = c(0,1), axes = F, pch = 19, col = c("Black", "Red"), main = "Sardine Fishery")
arrows(lowerS, x, 
       upperS, x, length=0, angle=90, code=3, lwd = 1,col = c("Black", "Red"))
#axis(2, at = c(1,2,3,4,5), labels = names, las = 2, pos = -0.1, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = 0.8)
#text(middleS, x, labels =round(middleS,2) )
#legend("bottomleft", col = c("Black", "Red"), pch = 19, legend = c("Specialist", "Generalist"))
#text(-0.4,5.4, label = "(B)", xpd = NA)
text(-0.1,0.2, "Average relative seabird abundance", xpd = NA, cex = 1.5)



lwdp = 3
par(mar=c(6,4,2,2)) # Rotate axis labels 
# margins: c(5, 4, 4, 2) + 0.1
adj = 0.05/4
plot(C1.vec+adj,type='l',col=hcr.colors[4],lwd=lwdp, ylim=c(0,0.9),
     axes = FALSE, xlab="",ylab="",xaxs="i",yaxs="i")
xticks = c(min(B),max(B))
yticks = c(0,Fmsy,0.8)
axis(side = 1, at = xticks,labels = c("Low","High"), las = 1)
axis(side = 2, at = yticks, labels = c("0","Fmsy"," "), las = 2)
mtext( "Biomass", side = 1, line = 2, cex = 0.75)
mtext( "Fishing \nrate", side = 2.5,  line = 2, cex = 0.75)

lines(C2.vec+adj+0.006, col=hcr.colors[5],lwd=lwdp)         # green
lines(C3.vec+adj+0.008, col=hcr.colors[6],lwd=lwdp)         # pale green
#lines(cfp.vec+adj+0.01,col=hcr.colors[6],lwd=lwdp)          # orange
#lines(constf.vec,col=add.alpha(hcr.colors[4],alpha = 0.6),lwd=lwdp) # red line, constF.high
#lines(constf.vec_lo,col=add.alpha(hcr.colors[5],alpha = 0.6),lwd=lwdp)  
legend(2,1.1,legend = c("A) Moderate","B) Low cut-off","C) High Max"),
       bty = "n",lwd=rep(lwdp,times=4),col=hcr.colors[c(4,5,6)],lty=rep(1, times=3), xpd = NA)


####### what does convential fishing levels do to the generalist?
sardine1 = seabirdout2
sardineMSY = seabirdoutmsy

sar1 =  rowMeans(sardine1$totalpop[,201:1000])
sarmsy = rowMeans(sardineMSY$totalpop[,201:1000]/sardine1$totalpop[,201:1000])

anchovy1 = seabirdout2
anchovyMSY = seabirdoutmsy
#anchovyMSY$totalpop[is.na(anchovyMSY$totalpop)] = 0


anc1 =  rowMeans(anchovy1$totalpop[,201:1000])
ancmsy = rowMeans(anchovyMSY$totalpop[,201:1000]/anchovy1$totalpop[,201:1000])

# menhaden1 = seabirdout2
# menhadenMSY = seabirdoutD2
# men1 =  rowMeans(menhaden1$totalpop[,200:1000])
# menmsy = rowMeans(menhadenMSY$totalpop[,200:1000]/menhaden1$totalpop[,200:1000])

middle = c(median(sarmsy), median(ancmsy))#, median(menmsy))
lower = c(quantile(sarmsy, 0.025), quantile(ancmsy, 0.025))#, quantile(menmsy, 0.025))
upper = c(quantile(sarmsy, 0.975), quantile(ancmsy, 0.975))#, quantile(menmsy, 0.975))

par(mfrow = c(1,1))
par(oma = c(2,2,1,1))
par(xpd = FALSE)
colors = gray.colors(n = 5)
allff = as.data.frame(cbind(sarmsy, ancmsy)) #, menmsy))
names(allff) = c("Sardine", "Anchovy")#, "Menhaden")
boxplot(allff, ylim = c(0,1), frame = F, xlab = "Constant Fishing at MSY")#c("darkgrey", "black","darkgrey", "black"),
#xlab = "Sardine prey")
#axis(side = 2, line = -1)
#mtext("Prey = Sardine", side = 3, line = 1)
mtext(side = 2, "Mean Relative Flexible Seabird Abundance", xpd = NA, line = 2.6)
mtext(side = 1, "Fishing Rate Fmsy", xpd = NA, line = 2)

#MSY menahden restricted
median(rowMeans(seabirdoutD$totalpop[,200:1000])/rowMeans(seabirdout$totalpop[,200:1000]))

par(mfrow = c(1,1))
par(oma = c(1,1,1,1))
par(xpd = FALSE)
colors = gray.colors(n = 5)
allff = as.data.frame(cbind(sarmsy, ancmsy))
names(allff) = c("Sardine", "Anchovy")
boxplot(allff, ylim = c(0,1), frame = F, xlab = "Constant Fishing at MSY")#c("darkgrey", "black","darkgrey", "black"),
#xlab = "Sardine prey")
axis(side = 2, line = -1)
#mtext("Prey = Sardine", side = 3, line = 1)
mtext(side = 2, "Mean Relative Seabird Abundance", xpd = NA, line = 2.6)
output = boxplot(allff)
