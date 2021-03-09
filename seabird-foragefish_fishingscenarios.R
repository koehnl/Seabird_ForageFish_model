# When sourced by other code, runs simulations of the seabird population with
# a forage fish prey fished under different harvest control rules
# runs for both a restricted seabird life history and a flexible seabird life history

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