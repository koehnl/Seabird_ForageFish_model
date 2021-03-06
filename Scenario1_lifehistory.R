# NO FISHING

#fish = nofish$penguinfood
fish = nofish$avgbiomass
B0init_low = rep(NA, length = nsims)
B02init_low = rep(NA, length = nsims)
P0noninit_low = rep(NA, length = nsims)

B0init_h = rep(NA, length = nsims)
B02init_h = rep(NA, length = nsims)
P0noninit_h = rep(NA, length = nsims)

#harsher 3-egg fledging response - 
#param_fledge_special = c(-0.3,30,0.5) 

# for(j in 1:nsims) {
#   B0init_low[j] = mean(fish[1,201:1001,j]*props_low[1,201:1001,j])
#   B02init_low[j] = mean(fish[2,201:1001,j]*props_low[2,201:1001,j])
#   P0noninit_low[j] = mean(fish[1,201:1001,j]*nonbreedprops_low[1,201:1001,j])
#   B0init_h[j] = mean(fish[1,201:1001,j]*props[1,201:1001,j])
#   B02init_h[j] = mean(fish[2,201:1001,j]*props[2,201:1001,j])
#   P0noninit_h[j] = mean(fish[1,201:1001,j]*nonbreedprops[1,201:1001,j])
# }
for(j in 1:nsims) {
  B0init_low[j] = mean(fish[j,201:1001]*props_low[1,201:1001,j])
  B02init_low[j] = mean(fish[j,201:1001]*props_low[2,201:1001,j])
  P0noninit_low[j] = mean(fish[j,201:1001]*nonbreedprops_low[1,201:1001,j])
  B0init_h[j] = mean(fish[j,201:1001]*props[1,201:1001,j])
  B02init_h[j] = mean(fish[j,201:1001]*props[2,201:1001,j])
  P0noninit_h[j] = mean(fish[j,201:1001]*nonbreedprops[1,201:1001,j])
}

B0init_low_mean = mean(B0init_low)
B02init_low_mean = mean(B02init_low)
P0noninit_low_mean = mean(P0noninit_low)

B0init_h_mean = mean(B0init_h)
B02init_h_mean = mean(B02init_h)
P0noninit_h_mean = mean(P0noninit_h)

# BASE
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

seabirdout = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                  recruits = matrix(nrow = nsims, ncol = Nyear),
                  totalpop = matrix(nrow = nsims, ncol = Nyear),
                  fledglings = matrix(nrow = nsims, ncol = Nyear),
                  eggs = matrix(nrow = nsims, ncol = Nyear),
                  testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init = rep(NA, length = nsims)
#B02init = rep(NA, length = nsims)
#P0noninit = rep(NA, length = nsims)

set.seed(216)
for(ss in 1:nsims) {
  #totalff = fish[,ss]
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  #B0init[ss] = mean(ffbiomass[1,,ss])
  #B02init[ss] = mean(ffbiomass[2,,ss])
  #P0noninit[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #specialist
  seabird = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                            clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                            param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                            param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                            param_juv = param_juv_special, param_fledge3 = param_fledge3_special)

  seabirdout[["breeders"]][ss,] = colSums(seabird$pop[2,5:31,])
  seabirdout[["recruits"]][ss,] = seabird$recruits
  seabirdout[["totalpop"]][ss,] = colSums(seabird$pop[2,1:31,])
  seabirdout[["fledglings"]][ss,] = seabird$pop[2,1,]
  seabirdout[["eggs"]][ss,] = seabird$pop[1,1,]
  seabirdout[["testing"]] = seabird$testing
}

# low variance
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

seabirdout2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                  recruits = matrix(nrow = nsims, ncol = Nyear),
                  totalpop = matrix(nrow = nsims, ncol = Nyear),
                  fledglings = matrix(nrow = nsims, ncol = Nyear),
                  eggs = matrix(nrow = nsims, ncol = Nyear),
                  testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init = rep(NA, length = nsims)
#B02init = rep(NA, length = nsims)
#P0noninit = rep(NA, length = nsims)

set.seed(216)
for(ss in 1:nsims) {
  #totalff = fish[,ss]
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props_low[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props_low[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops_low[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops_low[2,,ss]
  
  #B0init[ss] = mean(ffbiomass[1,,ss])
  #B02init[ss] = mean(ffbiomass[2,,ss])
  #P0noninit[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  #specialist
  seabird2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                            clutch = clutch[1],P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean, 
                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                            param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                            param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                            param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
 
  seabirdout2[["breeders"]][ss,] = colSums(seabird2$pop[2,5:31,])
  seabirdout2[["recruits"]][ss,] = seabird2$recruits
  seabirdout2[["totalpop"]][ss,] = colSums(seabird2$pop[2,1:31,])
  seabirdout2[["fledglings"]][ss,] = seabird2$pop[2,1,]
  seabirdout2[["eggs"]][ss,] = seabird2$pop[1,1,]
  seabirdout2[["testing"]] = seabird2$testing
}

# different clutch
maxage_use = maxage[1]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[1], 
                                          K = K_init, 
                                          maxage = maxage[1], egg_s = eggsur, chick_s = chicksur, clutch = clutch[2])

total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,450000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total

sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdout3 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                  recruits = matrix(nrow = nsims, ncol = Nyear),
                  totalpop = matrix(nrow = nsims, ncol = Nyear),
                  fledglings = matrix(nrow = nsims, ncol = Nyear),
                  eggs = matrix(nrow = nsims, ncol = Nyear),
                  testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init = rep(NA, length = nsims)
#B02init = rep(NA, length = nsims)
#P0noninit = rep(NA, length = nsims)

set.seed(216)
for(ss in 1:nsims) {
  #totalff = fish[,ss]
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  #B0init[ss] = mean(ffbiomass[1,,ss])
  #B02init[ss] = mean(ffbiomass[2,,ss])
  #P0noninit[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  #specialist
  seabird3 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                            clutch = clutch[2],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                            param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                            param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                            param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  # seabird = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
  #                                           juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
  #                                           K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
  #                                           clutch = clutch[1],P0 = B0init[ss], B0 = B0init[ss], P02 = B02init[ss], P0non = P0noninit[ss], 
  #                                           ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
  #                                           ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
  #                                           param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
  #                                           param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
  #                                           param_juv = param_juv_general, param_fledge3 = param_fledge3_general)  
  
  seabirdout3[["breeders"]][ss,] = colSums(seabird3$pop[2,5:31,])
  seabirdout3[["recruits"]][ss,] = seabird3$recruits
  seabirdout3[["totalpop"]][ss,] = colSums(seabird3$pop[2,1:31,])
  seabirdout3[["fledglings"]][ss,] = seabird3$pop[2,1,]
  seabirdout3[["eggs"]][ss,] = seabird3$pop[1,1,]
  seabirdout3[["testing"]] = seabird3$testing
}

# generalist instead of specialist
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

seabirdout4 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                  recruits = matrix(nrow = nsims, ncol = Nyear),
                  totalpop = matrix(nrow = nsims, ncol = Nyear),
                  fledglings = matrix(nrow = nsims, ncol = Nyear),
                  eggs = matrix(nrow = nsims, ncol = Nyear),
                  testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init = rep(NA, length = nsims)
#B02init = rep(NA, length = nsims)
#P0noninit = rep(NA, length = nsims)

set.seed(216)
for(ss in 1:nsims) {
  #totalff = fish[,ss]
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  #B0init[ss] = mean(ffbiomass[1,,ss])
  #B02init[ss] = mean(ffbiomass[2,,ss])
  #P0noninit[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  #specialist
  seabird4 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                            clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                            param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                            param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                            param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  # seabird = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
  #                                           juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
  #                                           K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
  #                                           clutch = clutch[1],P0 = B0init[ss], B0 = B0init[ss], P02 = B02init[ss], P0non = P0noninit[ss], 
  #                                           ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
  #                                           ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
  #                                           param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
  #                                           param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
  #                                           param_juv = param_juv_general, param_fledge3 = param_fledge3_general)  
  
  seabirdout4[["breeders"]][ss,] = colSums(seabird4$pop[2,5:31,])
  seabirdout4[["recruits"]][ss,] = seabird4$recruits
  seabirdout4[["totalpop"]][ss,] = colSums(seabird4$pop[2,1:31,])
  seabirdout4[["fledglings"]][ss,] = seabird4$pop[2,1,]
  seabirdout4[["eggs"]][ss,] = seabird4$pop[1,1,]
  seabirdout4[["testing"]] = seabird4$testing
}

# age at first breeding
maxage_use = maxage[1]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[2], 
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

seabirdout5 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                   recruits = matrix(nrow = nsims, ncol = Nyear),
                   totalpop = matrix(nrow = nsims, ncol = Nyear),
                   fledglings = matrix(nrow = nsims, ncol = Nyear),
                   eggs = matrix(nrow = nsims, ncol = Nyear),
                   testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init = rep(NA, length = nsims)
#B02init = rep(NA, length = nsims)
#P0noninit = rep(NA, length = nsims)

set.seed(216)
for(ss in 1:nsims) {
  #totalff = fish[,ss]
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  #B0init[ss] = mean(ffbiomass[1,,ss])
  #B02init[ss] = mean(ffbiomass[2,,ss])
  #P0noninit[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  #specialist
  seabird5 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                             K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdout5[["breeders"]][ss,] = colSums(seabird5$pop[2,(agebreeding[2]+1):(maxage[1]+1),])
  seabirdout5[["recruits"]][ss,] = seabird5$recruits
  seabirdout5[["totalpop"]][ss,] = colSums(seabird5$pop[2,1:31,])
  seabirdout5[["fledglings"]][ss,] = seabird5$pop[2,1,]
  seabirdout5[["eggs"]][ss,] = seabird5$pop[1,1,]
  seabirdout5[["testing"]] = seabird5$testing
}

# lifespan
maxage_use = maxage[2]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[1], 
                                          K = K_init, 
                                          maxage = maxage[2], egg_s = eggsur, chick_s = chicksur, clutch = clutch[1])

total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,450000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total

sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdout6 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                   recruits = matrix(nrow = nsims, ncol = Nyear),
                   totalpop = matrix(nrow = nsims, ncol = Nyear),
                   fledglings = matrix(nrow = nsims, ncol = Nyear),
                   eggs = matrix(nrow = nsims, ncol = Nyear),
                   testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init = rep(NA, length = nsims)
#B02init = rep(NA, length = nsims)
#P0noninit = rep(NA, length = nsims)

set.seed(216)
for(ss in 1:nsims) {
  #totalff = fish[,ss]
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  ffbiomass[1,,ss] = firsthalf*props[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops[2,,ss]
  
  #B0init[ss] = mean(ffbiomass[1,,ss])
  #B02init[ss] = mean(ffbiomass[2,,ss])
  #P0noninit[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  #specialist
  seabird6 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                             K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdout6[["breeders"]][ss,] = colSums(seabird6$pop[2,(agebreeding[1]+1):(maxage[2]+1),])
  seabirdout6[["recruits"]][ss,] = seabird6$recruits
  seabirdout6[["totalpop"]][ss,] = colSums(seabird6$pop[2,1:(maxage[2]+1),])
  seabirdout6[["fledglings"]][ss,] = seabird6$pop[2,1,]
  seabirdout6[["eggs"]][ss,] = seabird6$pop[1,1,]
  seabirdout6[["testing"]] = seabird6$testing
}

seabirdout$totalpop[is.na(seabirdout$totalpop)] = 0
seabirdout2$totalpop[is.na(seabirdout2$totalpop)] = 0
seabirdout3$totalpop[is.na(seabirdout3$totalpop)] = 0
seabirdout4$totalpop[is.na(seabirdout4$totalpop)] = 0

mean1 = colMeans(seabirdout$totalpop[,200:1000])
CI1 = apply(seabirdout$totalpop[,200:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean2 = colMeans(seabirdout2$totalpop[,200:1000])
CI2 = apply(seabirdout2$totalpop[,200:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean3 = colMeans(seabirdout3$totalpop[,200:1000])
CI3 = apply(seabirdout3$totalpop[,200:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean4 = colMeans(seabirdout4$totalpop[,200:1000])
CI4 = apply(seabirdout4$totalpop[,200:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

library(matrixStats)
seabird1total = seabirdout$totalpop[,200:1000]
quant1 = colQuantiles(seabird1total, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quant2 = colQuantiles(seabirdout2$totalpop[,200:1000], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quant3 = colQuantiles(seabirdout3$totalpop[,200:1000], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
quant4 = colQuantiles(seabirdout4$totalpop[,200:1000], probs = c(0.025, 0.25, 0.5, 0.75, 0.975))


#fish = constF$penguinfood
fish = constF$avgbiomass
# BASE
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
                   testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))


set.seed(216)
for(ss in 1:nsims) {
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] 
  secondhalf = fish[ss,]
  
  
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
  
  seabirdoutB[["breeders"]][ss,] = colSums(seabirdB$pop[2,5:31,])
  seabirdoutB[["recruits"]][ss,] = seabirdB$recruits
  seabirdoutB[["totalpop"]][ss,] = colSums(seabirdB$pop[2,1:31,])
  seabirdoutB[["fledglings"]][ss,] = seabirdB$pop[2,1,]
  seabirdoutB[["eggs"]][ss,] = seabirdB$pop[1,1,]
  seabirdoutB[["testing"]] = seabirdB$testing
}

### 2 - low variance instead of high
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

seabirdoutB2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = rep(NA, length = Nyear))

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
  
  ffbiomass[1,,ss] = firsthalf*props_low[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props_low[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops_low[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops_low[2,,ss]
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  
  #set.seed(821)
  # specialist
  seabirdB2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdoutB2[["breeders"]][ss,] = colSums(seabirdB2$pop[2,5:31,])
  seabirdoutB2[["recruits"]][ss,] = seabirdB2$recruits
  seabirdoutB2[["totalpop"]][ss,] = colSums(seabirdB2$pop[2,1:31,])
  seabirdoutB2[["fledglings"]][ss,] = seabirdB2$pop[2,1,]
  seabirdoutB2[["eggs"]][ss,] = seabirdB2$pop[1,1,]
  seabirdoutB2[["testing"]] = seabirdB2$testing
}

### 3 - different productivity - more eggs
maxage_use = maxage[1]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[2], 
                                          K = K_init, 
                                          maxage = maxage[1], egg_s = eggsur, chick_s = chicksur, clutch = clutch[2])

total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,450000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total

sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdoutB3 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = rep(NA, length = Nyear))

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
  seabirdB3 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[2],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdoutB3[["breeders"]][ss,] = colSums(seabirdB3$pop[2,5:31,])
  seabirdoutB3[["recruits"]][ss,] = seabirdB3$recruits
  seabirdoutB3[["totalpop"]][ss,] = colSums(seabirdB3$pop[2,1:31,])
  seabirdoutB3[["fledglings"]][ss,] = seabirdB3$pop[2,1,]
  seabirdoutB3[["eggs"]][ss,] = seabirdB3$pop[1,1,]
  seabirdoutB3[["testing"]] = seabirdB3$testing
}

### 4 - general functional instead of specialist
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

seabirdoutB4 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = rep(NA, length = Nyear))

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
  seabirdB4 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  
  seabirdoutB4[["breeders"]][ss,] = colSums(seabirdB4$pop[2,5:31,])
  seabirdoutB4[["recruits"]][ss,] = seabirdB4$recruits
  seabirdoutB4[["totalpop"]][ss,] = colSums(seabirdB4$pop[2,1:31,])
  seabirdoutB4[["fledglings"]][ss,] = seabirdB4$pop[2,1,]
  seabirdoutB4[["eggs"]][ss,] = seabirdB4$pop[1,1,]
  seabirdoutB4[["testing"]] = seabirdB4$testing
}

### 5 - age at first breeding
maxage_use = maxage[1]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[2], 
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

seabirdoutB5 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = rep(NA, length = Nyear))

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
  seabirdB5 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdoutB5[["breeders"]][ss,] = colSums(seabirdB5$pop[2,(agebreeding[2]+1):(maxage[1]+1),])
  seabirdoutB5[["recruits"]][ss,] = seabirdB5$recruits
  seabirdoutB5[["totalpop"]][ss,] = colSums(seabirdB5$pop[2,1:31,])
  seabirdoutB5[["fledglings"]][ss,] = seabirdB5$pop[2,1,]
  seabirdoutB5[["eggs"]][ss,] = seabirdB5$pop[1,1,]
  seabirdoutB5[["testing"]] = seabirdB5$testing
}

### 6 maxage
maxage_use = maxage[2]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])


stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, agebreeding = agebreeding[1], 
                                          K = K_init, 
                                          maxage = maxage[2], egg_s = eggsur, chick_s = chicksur, clutch = clutch[1])

total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]), ylim = c(0,450000))
stable = stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total

sum(stable[2:(maxage_use+1)])
testinits = stable*(1000000) # hypothetical starting population
#K = sum(testinits[2:(maxage_use+1)])
K = 1000000

seabirdoutB6 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                    recruits = matrix(nrow = nsims, ncol = Nyear),
                    totalpop = matrix(nrow = nsims, ncol = Nyear),
                    fledglings = matrix(nrow = nsims, ncol = Nyear),
                    eggs = matrix(nrow = nsims, ncol = Nyear),
                    testing = rep(NA, length = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0initB = B0init
#B02initB = B02init
#P0noninitB = P0noninit

set.seed(216)
for(ss in 1:nsims) {
 # firsthalf = fish[1,,ss] 
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
  seabirdB6 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdoutB6[["breeders"]][ss,] = colSums(seabirdB6$pop[2,(agebreeding[1]+1):(maxage[2]+1),])
  seabirdoutB6[["recruits"]][ss,] = seabirdB6$recruits
  seabirdoutB6[["totalpop"]][ss,] = colSums(seabirdB6$pop[2,1:(maxage[2]+1),])
  seabirdoutB6[["fledglings"]][ss,] = seabirdB6$pop[2,1,]
  seabirdoutB6[["eggs"]][ss,] = seabirdB6$pop[1,1,]
  seabirdoutB6[["testing"]] = seabirdB6$testing
}
seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutB2$totalpop[is.na(seabirdoutB2$totalpop)] = 0
seabirdoutB3$totalpop[is.na(seabirdoutB3$totalpop)] = 0
seabirdoutB4$totalpop[is.na(seabirdoutB4$totalpop)] = 0
seabirdoutB5$totalpop[is.na(seabirdoutB5$totalpop)] = 0
seabirdoutB6$totalpop[is.na(seabirdoutB6$totalpop)] = 0

par(mfrow = c(1,2))
par(mar= c(2,6,1,1))
par(oma = c(3,6,1,1))
prob1 = rep(0, length = nsims)
prob2 = rep(0, length = nsims)
prob3 = rep(0, length = nsims)
prob4 = rep(0, length = nsims)
prob5 = rep(0, length = nsims)
prob6 = rep(0, length = nsims)
thres = 0.5
for(i in 1:nsims) {
  prob1[i] = (length(which(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000] < thres)))/800
  prob2[i] = (length(which(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000] < thres)))/800
  prob3[i] = (length(which(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000] < thres)))/800
  prob4[i] = (length(which(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000] < thres)))/800 
  prob5[i] = (length(which(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000] < thres)))/800 
  prob6[i] = (length(which(seabirdoutB6$totalpop[i,201:1000]/seabirdout6$totalpop[i,201:1000] < thres)))/800 
  
}

probs = rbind(prob1,prob2, prob3, prob4)
probmed = c(quantile(prob1, 0.5), quantile(prob2, 0.5), quantile(prob3, 0.5), quantile(prob4, 0.5), quantile(prob5, 0.5),quantile(prob6, 0.5))
problow = c(quantile(prob1, 0.025), quantile(prob2, 0.025), quantile(prob3, 0.025), quantile(prob4, 0.025), quantile(prob5, 0.025),quantile(prob6, 0.025))
probhigh = c(quantile(prob1, 0.975), quantile(prob2, 0.975), quantile(prob3, 0.975), quantile(prob4, 0.975), quantile(prob5, 0.975),quantile(prob6, 0.975))
problow2 = c(quantile(prob1, 0.25), quantile(prob2, 0.25), quantile(prob3, 0.25), quantile(prob4, 0.25), quantile(prob5, 0.25),quantile(prob6, 0.25))
probhigh2 = c(quantile(prob1, 0.75), quantile(prob2, 0.75), quantile(prob3, 0.75), quantile(prob4, 0.75), quantile(prob5, 0.75),quantile(prob6, 0.75))


x = c(0,1,2,3,4)
names = c("Base Scenario 1", "Low Variance", "Low 1st Breeding Age", "Larger Clutch",
          "Generalist")
plot(1-probmed[c(4,3,5,2,1)], x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19)
arrows(1-problow[c(4,3,5,2,1)], x, 1-probhigh[c(4,3,5,2,1)], x, length=0, angle=90, code=3, lwd = 1)
arrows(1-problow2[c(4,3,5,2,1)], x, 1-probhigh2[c(4,3,5,2,1)], x, length=0, angle=90, code=3, lwd = 4)
axis(2, at = x, labels = c(names[5:1]), las = 1, pos = -0.2, cex.axis = 0.8)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
mtext("% Years with >50% of base \nmodel population",side = 1, line = 4, cex = 1)
mtext("(A)", side = 3, at = c(-0.1,4))

meanA = rep(0, length = nsims)
meanB0 = rep(0, length = nsims)
meanC = rep(0, length = nsims)
meanD = rep(0, length = nsims)
meanE = rep(0, length = nsims)
meanF = rep(0, length = nsims)
for(i in 1:nsims) {
  meanA[i] = mean(seabirdoutB$totalpop[i,200:1000]/seabirdout$totalpop[i,200:1000])
  meanB0[i] = mean(seabirdoutB2$totalpop[i,200:1000]/seabirdout2$totalpop[i,200:1000])
  meanC[i] = mean(seabirdoutB3$totalpop[i,200:1000]/seabirdout3$totalpop[i,200:1000])
  meanD[i] = mean(seabirdoutB4$totalpop[i,200:1000]/seabirdout4$totalpop[i,200:1000])
  meanE[i] = mean(seabirdoutB5$totalpop[i,200:1000]/seabirdout5$totalpop[i,200:1000])
  meanF[i] = mean(seabirdoutB6$totalpop[i,200:1000]/seabirdout6$totalpop[i,200:1000])
}
meanA[is.na(meanA)] = 0
meanB0[is.na(meanB0)] = 0
meanC[is.na(meanC)] = 0
meanD[is.na(meanD)] = 0
meanE[is.na(meanE)] = 0
meanF[is.na(meanF)] = 0

meanmed = c(quantile(meanA, 0.5), quantile(meanB0, 0.5), quantile(meanC, 0.5), quantile(meanD, 0.5),quantile(meanE, 0.5),quantile(meanF, 0.5))
meanlow = c(quantile(meanA, 0.025), quantile(meanB0, 0.025), quantile(meanC, 0.025), quantile(meanD, 0.025),quantile(meanE, 0.025),quantile(meanF, 0.025))
meanhigh = c(quantile(meanA, 0.975), quantile(meanB0, 0.975), quantile(meanC, 0.975), quantile(meanD, 0.975), quantile(meanE, 0.975),quantile(meanF, 0.975))
meanlow2 = c(quantile(meanA, 0.25), quantile(meanB0, 0.25), quantile(meanC, 0.25), quantile(meanD, 0.25), quantile(meanE, 0.25),quantile(meanF, 0.25))
meanhigh2 = c(quantile(meanA, 0.75), quantile(meanB0, 0.75), quantile(meanC, 0.75), quantile(meanD, 0.75), quantile(meanE, 0.75), quantile(meanF, 0.75))

x = c(0,1,2,3,4)
par(mar = c(6,10,4,2))
plot(meanmed[c(4,3,5,2,1)], x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19)
arrows(meanlow[c(4,3,5,2,1)], x, meanhigh[c(4,3,5,2,1)], x, length=0, angle=90, code=3, lwd = 1)
arrows(meanlow2[c(4,3,5,2,1)], x, meanhigh2[c(4,3,5,2,1)], x, length=0, angle=90, code=3, lwd = 4)
axis(2, at = x, labels = c(names[5:1]), las = 2, pos = -0.05, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
mtext("Average relative seabird \nabundance",side = 1, line = 3.5, cex = 1)
mtext("(B) Sardine Prey", side = 3, at = c(-0.3,10))

mean(colMeans(seabirdoutB$fledglings)/colMeans(seabirdout$fledglings))
mean(colMeans(seabirdoutB3$fledglings)/colMeans(seabirdout3$fledglings))

### FOR FIGURE 5 #########
####### NEW sardine and anchovy ######

seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutB2$totalpop[is.na(seabirdoutB2$totalpop)] = 0
seabirdoutB3$totalpop[is.na(seabirdoutB3$totalpop)] = 0
seabirdoutB4$totalpop[is.na(seabirdoutB4$totalpop)] = 0
seabirdoutB5$totalpop[is.na(seabirdoutB5$totalpop)] = 0
seabirdoutB6$totalpop[is.na(seabirdoutB6$totalpop)] = 0

#runs with one forage fish
meanA = rep(0, length = nsims)
meanB0 = rep(0, length = nsims)
meanC = rep(0, length = nsims)
meanD = rep(0, length = nsims)
meanE = rep(0, length = nsims)
meanF = rep(0, length = nsims)


for(i in 1:nsims) {
  meanA[i] = mean(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000])
  meanB0[i] = mean(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000])
  meanC[i] = mean(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000])
  meanD[i] = mean(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000])
  meanE[i] = mean(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000])
  meanF[i] = mean(seabirdoutB6$totalpop[i,201:1000]/seabirdout6$totalpop[i,201:1000])
}
meanA[is.na(meanA)] = 0
meanB0[is.na(meanB0)] = 0
meanC[is.na(meanC)] = 0
meanD[is.na(meanD)] = 0
meanE[is.na(meanE)] = 0
meanF[is.na(meanF)] = 0

# runs with second forage fish
meanA2 = rep(0, length = nsims)
meanB02 = rep(0, length = nsims)
meanC2 = rep(0, length = nsims)
meanD2 = rep(0, length = nsims)
meanE2 = rep(0, length = nsims)
meanF2 = rep(0, length = nsims)


for(i in 1:nsims) {
  meanA2[i] = mean(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000])
  meanB02[i] = mean(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000])
  meanC2[i] = mean(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000])
  meanD2[i] = mean(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000])
  meanE2[i] = mean(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000])
  meanF2[i] = mean(seabirdoutB6$totalpop[i,201:1000]/seabirdout6$totalpop[i,201:1000])
}
meanA2[is.na(meanA2)] = 0
meanB02[is.na(meanB02)] = 0
meanC2[is.na(meanC2)] = 0
meanD2[is.na(meanD2)] = 0
meanE2[is.na(meanE2)] = 0
meanF2[is.na(meanF2)] = 0

meanmed = c(quantile(meanA, 0.5), quantile(meanA2, 0.5), quantile(meanB0, 0.5), quantile(meanB02, 0.5),
            quantile(meanC, 0.5), quantile(meanC2, 0.5), quantile(meanD, 0.5), quantile(meanD2, 0.5),
            quantile(meanE, 0.5), quantile(meanE2, 0.5), quantile(meanF, 0.5), quantile(meanF2, 0.5))
meanlow = c(quantile(meanA, 0.025), quantile(meanA2, 0.025), quantile(meanB0, 0.025), quantile(meanB02, 0.025),
            quantile(meanC, 0.025), quantile(meanC2, 0.025), quantile(meanD, 0.025), quantile(meanD2, 0.025),
            quantile(meanE, 0.025), quantile(meanE2, 0.025),quantile(meanF, 0.025), quantile(meanF2, 0.025))
meanhigh = c(quantile(meanA, 0.975), quantile(meanA2, 0.975), quantile(meanB0, 0.975), quantile(meanB02, 0.975),
             quantile(meanC, 0.975), quantile(meanC2, 0.975), quantile(meanD, 0.975), quantile(meanD2, 0.975),
             quantile(meanE, 0.975), quantile(meanE2, 0.975), quantile(meanF, 0.975), quantile(meanF2, 0.975))
meanlow2 = c(quantile(meanA, 0.25), quantile(meanA2, 0.25), quantile(meanB0, 0.25), quantile(meanB02, 0.25),
             quantile(meanC, 0.25), quantile(meanC2, 0.25), quantile(meanD, 0.25), quantile(meanD2, 0.25),
             quantile(meanE, 0.25), quantile(meanE2, 0.25), quantile(meanF, 0.25), quantile(meanF2, 0.25))
meanhigh2 = c(quantile(meanA, 0.75), quantile(meanA2, 0.75), quantile(meanB0, 0.75), quantile(meanB02, 0.75),
              quantile(meanC, 0.75), quantile(meanC2, 0.75), quantile(meanD, 0.75), quantile(meanD2, 0.75),
              quantile(meanE, 0.75), quantile(meanE2, 0.75), quantile(meanF, 0.75), quantile(meanF2, 0.75))

##### BOTH LIFE HISTORY AND FUNCTIONAL RESPONSE SENSITIVITY ######
# which order plotted, anchovy vs. sardine, would need to match with legend
#par(mfrow = c(1,2))
names = c("Lower Max Age", "Base scenario", "Lower variance", "Lower age at\n first breeding", "Larger clutch",
          "All functional \nresponse generalist", "Reproductive success \ngeneralist",
          "Juvenile survival \ngeneralist","Breeder attendance\n generalist",
          "Adult survival \ngeneralist")
x = c(0,0.2,1,1.2,2,2.2,3,3.2,4,4.2, 5,5.2,6,6.2,7,7.2,8,8.2,9,9.2)
#par(mfrow = c(1,1))
par(mar = c(6,11,1,2))
plot(c(meanmed2[c(3,4,9,10,5,6,7,8)],meanmed[c(7,8,5,6,9,10,3,4,1,2,11,12)]), x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19, 
     col = c("black", "grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))

#text(meanmed[c(7,8,5,6,9,10,3,4,1,2)], x, labels = meanmed[c(7,8,5,6,9,10,3,4,1,2)] )

arrows(c(meanlow2B[c(3,4,9,10,5,6,7,8)], meanlow[c(7,8,5,6,9,10,3,4,1,2,11,12)]), x, 
       c(meanhigh2B[c(3,4,9,10,5,6,7,8)],meanhigh[c(7,8,5,6,9,10,3,4,1,2,11,12)]), x, length=0, angle=90, code=3, lwd = 1,
       col = c("black", "grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
arrows(c(meanlow22[c(3,4,9,10,5,6,7,8)],meanlow2[c(7,8,5,6,9,10,3,4,1,2,11,12)]), x,
       c(meanhigh22[c(3,4,9,10,5,6,7,8)],meanhigh2[c(7,8,5,6,9,10,3,4,1,2,11,12)]), x, length=0, angle=90, code=3, lwd = 4, 
       col = c("black", "grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
axis(2, at = c(0.1,1.1,2.1,3.1,4.1, 5.1,6.1,7.1,8.1,9.1), labels = c(names[10:1]), las = 2, pos = -0.05, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
mtext("Average relative seabird \nabundance",side = 1, line = 3.5, cex = 1)
#mtext("(B) Sardine Prey", side = 3, at = c(-0.3,10))
legend(-0.5,-1.4,xpd=NA, legend = c("Anchovy prey","Sardine prey"), col = c("grey50", "black"), lty = 1, pch = 1, lwd = 4)
#text(-0.6,4.5, label = "(A)", xpd = NA, cex = 1.5)
abline(h = 3.5, lty = 2)

meanmed2 = c(quantile(meanAX, 0.5), quantile(meanA2X, 0.5), quantile(meanB0X, 0.5), quantile(meanB02X, 0.5),
             quantile(meanCX, 0.5), quantile(meanC2X, 0.5), quantile(meanDX, 0.5), quantile(meanD2X, 0.5),
             quantile(meanEX, 0.5), quantile(meanE2X, 0.5))
meanlow2B = c(quantile(meanAX, 0.025), quantile(meanA2X, 0.025), quantile(meanB0X, 0.025), quantile(meanB02X, 0.025),
              quantile(meanCX, 0.025), quantile(meanC2X, 0.025), quantile(meanDX, 0.025), quantile(meanD2X, 0.025),
              quantile(meanEX, 0.025), quantile(meanE2X, 0.025))
meanhigh2B = c(quantile(meanAX, 0.975), quantile(meanA2X, 0.975), quantile(meanB0X, 0.975), quantile(meanB02X, 0.975),
               quantile(meanCX, 0.975), quantile(meanC2X, 0.975), quantile(meanDX, 0.975), quantile(meanD2X, 0.975),
               quantile(meanEX, 0.975), quantile(meanE2X, 0.975))
meanlow22 = c(quantile(meanAX, 0.25), quantile(meanA2X, 0.25), quantile(meanB0X, 0.25), quantile(meanB02X, 0.25),
              quantile(meanCX, 0.25), quantile(meanC2X, 0.25), quantile(meanDX, 0.25), quantile(meanD2X, 0.25),
              quantile(meanEX, 0.25), quantile(meanE2X, 0.25))
meanhigh22 = c(quantile(meanAX, 0.75), quantile(meanA2X, 0.75), quantile(meanB0X, 0.75), quantile(meanB02X, 0.75),
               quantile(meanCX, 0.75), quantile(meanC2X, 0.75), quantile(meanDX, 0.75), quantile(meanD2X, 0.75),
               quantile(meanEX, 0.75), quantile(meanE2X, 0.75))

## JUST LIFE HISTORY SENSITIVITY ####
x = c(0,0.2,1,1.2,2,2.2,3,3.2,4,4.2)
#par(mar = c(6,8,2,2))
plot(meanmed2[c(3,4,9,10,5,6,7,8, 1,2)], x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19, 
     col = c("black" ,"grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))

#text(meanmed2[c(4,3,10,9,6,5,8, 7,2, 1)], x, labels = meanmed2[c(4,3,10,9,6,5,8, 7,2, 1)])
arrows(meanlow2B[c(3,4,9,10,5,6,7,8, 1,2)], x, meanhigh2B[c(3,4,9,10,5,6,7,8, 1,2)], x, length=0, angle=90, code=3, lwd = 1,
       col = c("black","grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
arrows(meanlow22[c(3,4,9,10,5,6,7,8, 1,2)], x, meanhigh22[c(3,4,9,10,5,6,7,8, 1,2)], x, length=0, angle=90, code=3, lwd = 4, 
       col = c("black","grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
axis(2, at = c(0.1,1.1,2.1,3.1,4.1), labels = c(names[5:1]), las = 2, pos = -0.05, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
#mtext("Average relative seabird \nabundance",side = 1, line = 3.5, cex = 1)
#mtext("(B) Sardine Prey", side = 3, at = c(-0.3,10))
#legend(-0.55,-0.5,xpd=NA, legend = c("Anchovy prey", "Sardine prey"), col = c("Black", "grey50"), lty = 1, pch = 1, lwd = 4)
text(-0.8,4.5, label = "(B)", xpd = NA, cex = 1.5)


mtext('Average relative seabird abundance', side = 1, outer = TRUE, line = -2.5)
