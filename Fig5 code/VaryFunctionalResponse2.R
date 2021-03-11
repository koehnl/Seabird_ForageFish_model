# NO FISHING

#fish = nofish$penguinfood
fish = nofish$avgbiomass
B0init_low = rep(NA, length = nsims)
B02init_low = rep(NA, length = nsims)
P0noninit_low = rep(NA, length = nsims)

B0init_h = rep(NA, length = nsims)
B02init_h = rep(NA, length = nsims)
P0noninit_h = rep(NA, length = nsims)

# testing
#param_breeders_special = c(0.1, 20, 0.1)
#param_breeders_general = c(0.6, 20, 0.1)


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

set.seed(216, sample.kind = "Rejection")
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

# adult generalist
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

set.seed(216, sample.kind = "Rejection")
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
  seabird2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                             K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_general, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdout2[["breeders"]][ss,] = colSums(seabird2$pop[2,5:31,])
  seabirdout2[["recruits"]][ss,] = seabird2$recruits
  seabirdout2[["totalpop"]][ss,] = colSums(seabird2$pop[2,1:31,])
  seabirdout2[["fledglings"]][ss,] = seabird2$pop[2,1,]
  seabirdout2[["eggs"]][ss,] = seabird2$pop[1,1,]
  seabirdout2[["testing"]] = seabird2$testing
}

# juvenile general 
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

set.seed(216,sample.kind = "Rejection")
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
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                             param_juv = param_juv_general, param_fledge3 = param_fledge3_special)
 
  seabirdout3[["breeders"]][ss,] = colSums(seabird3$pop[2,5:31,])
  seabirdout3[["recruits"]][ss,] = seabird3$recruits
  seabirdout3[["totalpop"]][ss,] = colSums(seabird3$pop[2,1:31,])
  seabirdout3[["fledglings"]][ss,] = seabird3$pop[2,1,]
  seabirdout3[["eggs"]][ss,] = seabird3$pop[1,1,]
  seabirdout3[["testing"]] = seabird3$testing
}

# fledge general
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
                                             param_breeders = param_breeders_special, param_fledge = param_fledge_general, 
                                             param_fledge2 = param_fledge2_general, param_adult = param_adult_special, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_general)
  
  seabirdout4[["breeders"]][ss,] = colSums(seabird4$pop[2,5:31,])
  seabirdout4[["recruits"]][ss,] = seabird4$recruits
  seabirdout4[["totalpop"]][ss,] = colSums(seabird4$pop[2,1:31,])
  seabirdout4[["fledglings"]][ss,] = seabird4$pop[2,1,]
  seabirdout4[["eggs"]][ss,] = seabird4$pop[1,1,]
  seabirdout4[["testing"]] = seabird4$testing
}

# breeder attendance general
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
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                             K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_general, param_fledge = param_fledge_special, 
                                             param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                             param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdout5[["breeders"]][ss,] = colSums(seabird5$pop[2,5:31,])
  seabirdout5[["recruits"]][ss,] = seabird5$recruits
  seabirdout5[["totalpop"]][ss,] = colSums(seabird5$pop[2,1:31,])
  seabirdout5[["fledglings"]][ss,] = seabird5$pop[2,1,]
  seabirdout5[["eggs"]][ss,] = seabird5$pop[1,1,]
  seabirdout5[["testing"]] = seabird5$testing
}


mean1 = colMeans(seabirdout$totalpop[,201:1000])
CI1 = apply(seabirdout$totalpop[,201:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean2 = colMeans(seabirdout2$totalpop[,201:1000])
CI2 = apply(seabirdout2$totalpop[,201:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean3 = colMeans(seabirdout3$totalpop[,201:1000])
CI3 = apply(seabirdout3$totalpop[,201:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})

mean4 = colMeans(seabirdout4$totalpop[,201:1000])
CI4 = apply(seabirdout4$totalpop[,201:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
mean5 = colMeans(seabirdout5$totalpop[,201:1000])
CI5 = apply(seabirdout5$totalpop[,201:1000], 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})



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


set.seed(216, sample.kind = "Rejection")
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
  print(length(which(seabirdB$pop < 0)))
  seabirdoutB[["breeders"]][ss,] = colSums(seabirdB$pop[2,5:31,])
  seabirdoutB[["recruits"]][ss,] = seabirdB$recruits
  seabirdoutB[["totalpop"]][ss,] = colSums(seabirdB$pop[2,1:31,])
  seabirdoutB[["fledglings"]][ss,] = seabirdB$pop[2,1,]
  seabirdoutB[["eggs"]][ss,] = seabirdB$pop[1,1,]
  seabirdoutB[["testing"]] = seabirdB$testing
}

### adult general
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

set.seed(216, sample.kind = "Rejection")
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
  seabirdB2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_general, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdoutB2[["breeders"]][ss,] = colSums(seabirdB2$pop[2,5:31,])
  seabirdoutB2[["recruits"]][ss,] = seabirdB2$recruits
  seabirdoutB2[["totalpop"]][ss,] = colSums(seabirdB2$pop[2,1:31,])
  seabirdoutB2[["fledglings"]][ss,] = seabirdB2$pop[2,1,]
  seabirdoutB2[["eggs"]][ss,] = seabirdB2$pop[1,1,]
  seabirdoutB2[["testing"]] = seabirdB2$testing
}

### juvenile general
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

set.seed(216, sample.kind = "Rejection")
for(ss in 1:nsims) {
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
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
  seabirdB3 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                              param_juv = param_juv_general, param_fledge3 = param_fledge3_special)
  
  seabirdoutB3[["breeders"]][ss,] = colSums(seabirdB3$pop[2,5:31,])
  seabirdoutB3[["recruits"]][ss,] = seabirdB3$recruits
  seabirdoutB3[["totalpop"]][ss,] = colSums(seabirdB3$pop[2,1:31,])
  seabirdoutB3[["fledglings"]][ss,] = seabirdB3$pop[2,1,]
  seabirdoutB3[["eggs"]][ss,] = seabirdB3$pop[1,1,]
  seabirdoutB3[["testing"]] = seabirdB3$testing
}

### fledging general
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
 # firsthalf = fish[1,,ss] 
#  secondhalf = fish[2,,ss]
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
  seabirdB4 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_special, param_fledge = param_fledge_general, 
                                              param_fledge2 = param_fledge2_general, param_adult = param_adult_special, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_general)
  
  seabirdoutB4[["breeders"]][ss,] = colSums(seabirdB4$pop[2,5:31,])
  seabirdoutB4[["recruits"]][ss,] = seabirdB4$recruits
  seabirdoutB4[["totalpop"]][ss,] = colSums(seabirdB4$pop[2,1:31,])
  seabirdoutB4[["fledglings"]][ss,] = seabirdB4$pop[2,1,]
  seabirdoutB4[["eggs"]][ss,] = seabirdB4$pop[1,1,]
  seabirdoutB4[["testing"]] = seabirdB4$testing
}

### breeding attendance general 
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
  seabirdB5 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                              juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                              K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                              clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                              ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                              ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                              param_breeders = param_breeders_general, param_fledge = param_fledge_special, 
                                              param_fledge2 = param_fledge2_special, param_adult = param_adult_special, 
                                              param_juv = param_juv_special, param_fledge3 = param_fledge3_special)
  
  seabirdoutB5[["breeders"]][ss,] = colSums(seabirdB5$pop[2,5:31,])
  seabirdoutB5[["recruits"]][ss,] = seabirdB5$recruits
  seabirdoutB5[["totalpop"]][ss,] = colSums(seabirdB5$pop[2,1:31,])
  seabirdoutB5[["fledglings"]][ss,] = seabirdB5$pop[2,1,]
  seabirdoutB5[["eggs"]][ss,] = seabirdB5$pop[1,1,]
  seabirdoutB5[["testing"]] = seabirdB5$testing
}

seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutB2$totalpop[is.na(seabirdoutB2$totalpop)] = 0
seabirdoutB3$totalpop[is.na(seabirdoutB3$totalpop)] = 0
seabirdoutB4$totalpop[is.na(seabirdoutB4$totalpop)] = 0
seabirdoutB5$totalpop[is.na(seabirdoutB5$totalpop)] = 0
#seabirdoutB6$totalpop[is.na(seabirdoutB6$totalpop)] = 0

# OLD CODE #########
# par(mfrow = c(1,2))
# par(mar= c(2,6,2,1))
# par(oma = c(3,6,1,1))
# # prob1 = rep(0, length = nsims)
# # prob2 = rep(0, length = nsims)
# # prob3 = rep(0, length = nsims)
# # prob4 = rep(0, length = nsims)
# # prob5 = rep(0, length = nsims)
# # thres = 0.5
# # 
# # seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
# # seabirdoutB2$totalpop[is.na(seabirdoutB2$totalpop)] = 0
# # seabirdoutB3$totalpop[is.na(seabirdoutB3$totalpop)] = 0
# # seabirdoutB4$totalpop[is.na(seabirdoutB4$totalpop)] = 0
# # seabirdoutB5$totalpop[is.na(seabirdoutB5$totalpop)] = 0
# # 
# # for(i in 1:nsims) {
# #   prob1[i] = (length(which(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000] < thres)))/800
# #   prob2[i] = (length(which(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000] < thres)))/800
# #   prob3[i] = (length(which(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000] < thres)))/800
# #   prob4[i] = (length(which(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000] < thres)))/800 
# #   prob5[i] = (length(which(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000] < thres)))/800 
# # }
# # 
# # probs = rbind(prob1,prob2, prob3, prob4, prob5)
# # probmed = c(quantile(prob1, 0.5), quantile(prob2, 0.5), quantile(prob3, 0.5), quantile(prob4, 0.5),quantile(prob5, 0.5) )
# # problow = c(quantile(prob1, 0.025), quantile(prob2, 0.025), quantile(prob3, 0.025), quantile(prob4, 0.025), quantile(prob5, 0.025))
# # probhigh = c(quantile(prob1, 0.975), quantile(prob2, 0.975), quantile(prob3, 0.975), quantile(prob4, 0.975), quantile(prob5, 0.975))
# # problow2 = c(quantile(prob1, 0.25), quantile(prob2, 0.25), quantile(prob3, 0.25), quantile(prob4, 0.25),quantile(prob5, 0.25))
# # probhigh2 = c(quantile(prob1, 0.75), quantile(prob2, 0.75), quantile(prob3, 0.75), quantile(prob4, 0.75), quantile(prob5, 0.75))
# # 
# # 
# # x = c(0,1,2,3, 4)
# # names = c("Base Scenario 1", "Adult S. general", "Juv. S. general",
# #           "Fledging S. general", "Attendance general")
# # plot(c(1-probmed[5], 1-probmed[2:4], 1-probmed[1]), x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19)
# # arrows(c(1-problow[5], 1-problow[2:4], 1-problow[1]), x, 
# #        c(1-probhigh[5], 1-probhigh[2:4], 1-probhigh[1]), x, length=0, angle=90, code=3, lwd = 1)
# # arrows(c(1-problow2[5], 1-problow2[2:4], 1-problow2[1]), x,
# #        c(1-probhigh2[5], 1-probhigh2[2:4], 1-probhigh2[1]), x, length=0, angle=90, code=3, lwd = 4)
# # axis(2, at = x, labels = c(names[5], names[2:4], names[1]), las = 2, pos = -0.2, cex.axis = 1)
# # axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
# # mtext("% Years with >50% of base \nmodel population",side = 1, line = 3.5, cex = 1)
# # mtext("(A)", side = 3, at = c(-0.1,4))
# # 
# # 
# # meanA = rep(0, length = nsims)
# # meanB0 = rep(0, length = nsims)
# # meanC = rep(0, length = nsims)
# # meanD = rep(0, length = nsims)
# # meanE = rep(0, length = nsims)
# # for(i in 1:nsims) {
# #   meanA[i] = mean(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000])
# #   meanB0[i] = mean(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000])
# #   meanC[i] = mean(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000])
# #   meanD[i] = mean(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000])
# #   meanE[i] = mean(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000])
# # }
# # 
# # meanmed = c(quantile(meanA, 0.5), quantile(meanB0, 0.5), quantile(meanC, 0.5), quantile(meanD, 0.5),quantile(meanE, 0.5))
# # meanlow = c(quantile(meanA, 0.025), quantile(meanB0, 0.025), quantile(meanC, 0.025), quantile(meanD, 0.025), quantile(meanE, 0.025))
# # meanhigh = c(quantile(meanA, 0.975), quantile(meanB0, 0.975), quantile(meanC, 0.975), quantile(meanD, 0.975), quantile(meanE, 0.975))
# # meanlow2 = c(quantile(meanA, 0.25), quantile(meanB0, 0.25), quantile(meanC, 0.25), quantile(meanD, 0.25), quantile(meanE, 0.25))
# # meanhigh2 = c(quantile(meanA, 0.75), quantile(meanB0, 0.75), quantile(meanC, 0.75), quantile(meanD, 0.75), quantile(meanE, 0.75))
# # 
# # 
# # plot(c(meanmed[5], meanmed[2:4], meanmed[1]), x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19)
# # arrows(c(meanlow[5], meanlow[2:4], meanlow[1]), x, 
# #        c(meanhigh[5], meanhigh[2:4], meanhigh[1]), x, length=0, angle=90, code=3, lwd = 1)
# # arrows(c(meanlow2[5], meanlow2[2:4], meanlow2[1]), x, 
# #        c(meanhigh2[5], meanhigh2[2:4], meanhigh2[1]), x, length=0, angle=90, code=3, lwd = 4)
# # axis(2, at = x, labels = c(names2[5], names2[2:4], names2[1]), las = 2, pos = -0.2, cex.axis = 1)
# # axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
# # mtext("Average relative seabird \nabundance",side = 1, line = 3.5, cex = 1)
# # mtext("(B)", side = 3, at = c(-0.1,4))
# # 
# # names2 = c("Base Scenario 1", "Adult Sur. general", "Juv. Sur. general",
# #           "Fledge Sur. general", "Attendance general")
# 
# ############
# ####### New - both sardine and anchovy #########
# ## FUNCTIONAL RESPONSE SENSITIVITY #####
# # run for one forage fish prey
# 
# 
# meanAX = rep(0, length = nsims)
# meanB0X = rep(0, length = nsims)
# meanCX = rep(0, length = nsims)
# meanDX = rep(0, length = nsims)
# meanEX = rep(0, length = nsims)
# #meanFX = rep(0, length = nsims)
# 
# 
# for(i in 1:nsims) {
#   meanAX[i] = mean(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000])
#   meanB0X[i] = mean(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000])
#   meanCX[i] = mean(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000])
#   meanDX[i] = mean(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000])
#   meanEX[i] = mean(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000])
#   #meanFX[i] = mean(seabirdoutB6$totalpop[i,201:1000]/seabirdout6$totalpop[i,201:1000])
# }
# meanAX[is.na(meanAX)] = 0
# meanB0X[is.na(meanB0X)] = 0
# meanCX[is.na(meanCX)] = 0
# meanDX[is.na(meanDX)] = 0
# meanEX[is.na(meanEX)] = 0
# #meanFX[is.na(meanFX)] = 0
# 
# # second forage fish 
# meanA2X = rep(0, length = nsims)
# meanB02X = rep(0, length = nsims)
# meanC2X = rep(0, length = nsims)
# meanD2X = rep(0, length = nsims)
# meanE2X = rep(0, length = nsims)
# #meanF2X = rep(0, length = nsims)
# 
# 
# for(i in 1:nsims) {
#   meanA2X[i] = mean(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000])
#   meanB02X[i] = mean(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000])
#   meanC2X[i] = mean(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000])
#   meanD2X[i] = mean(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000])
#   meanE2X[i] = mean(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000])
#  # meanF2X[i] = mean(seabirdoutB6$totalpop[i,201:1000]/seabirdout6$totalpop[i,201:1000])
# }
# meanA2X[is.na(meanA2X)] = 0
# meanB02X[is.na(meanB02X)] = 0
# meanC2X[is.na(meanC2X)] = 0
# meanD2X[is.na(meanD2X)] = 0
# meanE2X[is.na(meanE2X)] = 0
# #meanF2X[is.na(meanF2X)] = 0
# # base, adult, juvenile, fledge, breeding 
# meanmed2 = c(quantile(meanAX, 0.5), quantile(meanA2X, 0.5), quantile(meanB0X, 0.5), quantile(meanB02X, 0.5),
#             quantile(meanCX, 0.5), quantile(meanC2X, 0.5), quantile(meanDX, 0.5), quantile(meanD2X, 0.5),
#             quantile(meanEX, 0.5), quantile(meanE2X, 0.5))
# meanlow2B = c(quantile(meanAX, 0.025), quantile(meanA2X, 0.025), quantile(meanB0X, 0.025), quantile(meanB02X, 0.025),
#             quantile(meanCX, 0.025), quantile(meanC2X, 0.025), quantile(meanDX, 0.025), quantile(meanD2X, 0.025),
#             quantile(meanEX, 0.025), quantile(meanE2X, 0.025))
# meanhigh2B = c(quantile(meanAX, 0.975), quantile(meanA2X, 0.975), quantile(meanB0X, 0.975), quantile(meanB02X, 0.975),
#              quantile(meanCX, 0.975), quantile(meanC2X, 0.975), quantile(meanDX, 0.975), quantile(meanD2X, 0.975),
#              quantile(meanEX, 0.975), quantile(meanE2X, 0.975))
# meanlow22 = c(quantile(meanAX, 0.25), quantile(meanA2X, 0.25), quantile(meanB0X, 0.25), quantile(meanB02X, 0.25),
#              quantile(meanCX, 0.25), quantile(meanC2X, 0.25), quantile(meanDX, 0.25), quantile(meanD2X, 0.25),
#              quantile(meanEX, 0.25), quantile(meanE2X, 0.25))
# meanhigh22 = c(quantile(meanAX, 0.75), quantile(meanA2X, 0.75), quantile(meanB0X, 0.75), quantile(meanB02X, 0.75),
#               quantile(meanCX, 0.75), quantile(meanC2X, 0.75), quantile(meanDX, 0.75), quantile(meanD2X, 0.75),
#               quantile(meanEX, 0.75), quantile(meanE2X, 0.75))
# 
# names = c("Base Scenario", "Fledge Surv. general",
#           "Juv. Surv. general","Breeder attendance\n general",
#           "Adult Surv. general")
# x = c(0,0.2,1,1.2,2,2.2,3,3.2,4,4.2)
# par(mar = c(6,8,2,2))
# plot(meanmed2[c(3,4,9,10,5,6,7,8, 1,2)], x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19, 
#      col = c("black" ,"grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
# 
# #text(meanmed2[c(4,3,10,9,6,5,8, 7,2, 1)], x, labels = meanmed2[c(4,3,10,9,6,5,8, 7,2, 1)])
# arrows(meanlow2B[c(3,4,9,10,5,6,7,8, 1,2)], x, meanhigh2B[c(3,4,9,10,5,6,7,8, 1,2)], x, length=0, angle=90, code=3, lwd = 1,
#        col = c("black","grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
# arrows(meanlow22[c(3,4,9,10,5,6,7,8, 1,2)], x, meanhigh22[c(3,4,9,10,5,6,7,8, 1,2)], x, length=0, angle=90, code=3, lwd = 4, 
#        col = c("black","grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
# axis(2, at = c(0.1,1.1,2.1,3.1,4.1), labels = c(names[5:1]), las = 2, pos = -0.05, cex.axis = 1)
# axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
# #mtext("Average relative seabird \nabundance",side = 1, line = 3.5, cex = 1)
# #mtext("(B) Sardine Prey", side = 3, at = c(-0.3,10))
# #legend(-0.55,-0.5,xpd=NA, legend = c("Anchovy prey", "Sardine prey"), col = c("Black", "grey50"), lty = 1, pch = 1, lwd = 4)
# text(-0.8,4.5, label = "(B)", xpd = NA, cex = 1.5)
# 
