# Run the basic seabird model with a non-fished forage fish prey
library(here)
setwd(here::here())

#directory = "laurakoehn"
#basedir_s <- paste("/Users/", directory, "/Dropbox/ff-mse2/Seabirds/", sep ="")
# Main universal parameters
Nyear = 1000
nsims = 100 # probably need 1000, trying 100 to test

# Set up seabird parameters and scenarios - NEED TO ADD SPECIALIST AND MIDDLE SCENARIOS
# FUNCTIONAL RESPONSES
# Parameters that vary
#param_breeders_special = c(0.1,20,0.3) # colony/breeder attendance (not survival)
#param_breeders_general = c(0.6,20,0.3)

# TRY these 2/12/2020 ####
param_breeders_special = c(0.2,30,0.3) # colony/breeder attendance (not survival)
param_breeders_general = c(0.6,30,0.3)
######

param_breeders_mid = c(0.3,20,0.3)
# OK? based on Piatt et al. 2007/Cairns should be shifted even further right.
# right shape? (see Harding et al. 2007 common murre)

param_fledge_special = c(-0.3,30,0.2)
# enough like 1/3 for the birds - i think so, and looks like Andre's
param_fledge_general = c(0.51,15,0.1) # so at some lower survival, switch and survival goes back to 1?
param_fledge_mid = c(0.27,20,0.15)

param_fledge2_special = c(-0.3,30,0.15)
param_fledge2_general = c(0.41,15,0.05) # *will 2nd and 3rd chick still apply for generalist who just switched anyway?
param_fledge2_mid = c(0.21,20,0.1)

param_fledge3_special = c(-0.3,30,0.1) 
param_fledge3_general = c(0.2,15,0) 
param_fledge3_mid = c(0.05,20,0.05)

fledge_harsh = c(-0.3, 30, 0.5)

#param_adult_special = c(0.1,20,0.15)  # pretty good based on Piatt et al. 2007
#param_adult_general = c(0.6,20,0.15)
param_adult_mid = c(0.3,20,0.15)

# TRY 2/12/2020###3
param_adult_special = c(0.2,30,0.15)  # should be 0.2,30,0.15 pretty good based on Piatt et al. 2007
param_adult_general = c(0.6,30,0.15)
#####

#param_juv_special = c(0.1,10,0.3) # need lots of prey or else only the strong survive
#param_juv_general = c(0.6,10,0.3)

# TRY 2/12/2020####
#param_juv_special = c(0.1,20,0.3) #the 0.3 is maybe making them very sensitive...
#param_juv_general = c(0.6,20,0.3)
# 3/3/2020
param_juv_special = c(0.2,30,0.2) #the 0.3 is maybe making them very sensitive...
param_juv_general = c(0.6,30,0.2)
# good right? cause probably don't want the same drop off (0.3) as the breeding
# attendance functional response which according to Piatt and Cairns would be
# shifted to the right from survivals (at least adult survival)
#####3

param_juv_mid = c(0.3,10,0.3)

func_flat = c(1, 0,0)

# life history parameters
agebreeding = c(5,3,5) 
maxage = c(30,15,30) 
clutch = c(1,3,1) 

#for testing
#agebreeding = c(3,3,3,3,3,3) #temp
#maxage = c(30,30,30,30,30,30) #temp 
#clutch = c(3,3,3,3,3,3)

#depth_threshold = c(25, 100, 100, 15, 100, 15) # double check - 100 for all depths, 25 for medium depths, 15 for shallow depths
depth_threshold = c(25,100,100,15,100,15)
dist_threshold = c(0.75, 0.5, 0.75, 0.5, 0.75, 0.3) #0.75 for far, 0.5 for medium, 0.3 for near shore?
# maybe actually only 2 distance thresholds....

# scenarios setting all prey accessibility the same
#depth_threshold = c(100,100, 100,100,100,100)
#dist_threshold = c(0.75, 0.75, 0.75,0.75, 0.75,0.75)

# scenarios = matrix(NA, nrow = 6, ncol = 6) #should be 6 scenarios so 6 parameter combinations 
# colnames(scenarios) = c("depth", "distance", "diet_type", "breeding_age", "max_age", "clutch")
# scenarios[,1] = c(100, 30, 100, 5, 100, 5)
# scenarios[,2] = c(NA, 0.5, 0.75, 0.5, 0.75, 0.25)
# scenarios[,3] = c(1, 2,3,1,3,NA) #1 = special, 2 = general, 3 = intermediate (mid)
# scenarios[,4] =  c(5, 3, 5, 2, 3, NA)
# scenarios[,5] = c(30,15, 30, 10, 10, NA)
# scenarios[,6] = c(1,3, 2,2,1, NA)
# 
# nscenarios = nrow(scenarios)

# Constant parameters
#diet_pref = 0.5
#Other = 1-diet_pref
# keep diet prop the same across scenarios, just change functional response? Yes for now   


estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

########################################################################################

#######################################################################
# Run each life history scenario
source(file.path("SeabirdModel_use/SeabirdStableAgeDistribution_GENERAL.R"))
source(file.path("SeabirdModel_use/Seabirdmodel_stochastic_general_2019.R"))
nsims = 100

# fishpoulation function in "fishingsimulations2019.R"
# fish = array(NA, c(2,Nyear+1,nsims))
# set.seed(10)
# #start = runif(nsims, 0, 2*pi)
# #start = rep(1, length = nsims)
# start = seq(0, 2*pi, length.out = 10)
# for(l in 1:nsims) {
#   unfish = fishpopulation(fishing = 0, harvest = 1, Blim = 0, Btarget = 0, Fmax = 0, delay = 1, years = 1000, sim = start[l], seed.index = 5)
#   fish[,,l] = unfish$biomass
# }

# Survival rates that are constant across scenarios
year1sur = 0.7
year1var = 0.005 #base is 0.005
adultsur = 0.92 #use 0.92, original is 0.9
adultvar = 0.001 # original 0.005 
chicksur = 0.76 #  0.75 is around the average from seabirds in my scenarios from Nelson 1980 Seabird book
# values from Nelson 1980 for chick survival but ignoring anomolous years since that will be captures in stochasticity
chickvar = 0.01 #0.05 too much variability for stochastic runs - leads to too low values
# this variance puts lowest possible at around 0.3, which is good. goes lower but due to prey not other factors

eggsur = 0.7 # seems low but is an average from numbers taken from Nelson 1980 for birds in my scenarios # make stochastic 40% to 95% ish 
eggvar = 0.01

adult_base_s = array(NA, c(2,Nyear+1,nsims))
juv_base_s = array(NA, c(2,Nyear+1,nsims))
egg_base_s = matrix(NA, nrow = (Nyear+1), ncol = nsims)
chick_base_s = matrix(NA, nrow = (Nyear+1), ncol = nsims)
set.seed(2) #originalsetseed 2
for(q in 1:(nsims)) {
  adult_s1 = estBetaParams(mu = adultsur, var = adultvar)
  adult_base_s[1,,q] = rbeta((Nyear+1), adult_s1$alpha, adult_s1$beta)
  adult_s2 = estBetaParams(mu = adultsur, var = adultvar)
  adult_base_s[2,,q] = rbeta((Nyear+1), adult_s2$alpha, adult_s2$beta)
  juv_s1 = estBetaParams(mu = year1sur, var = year1var)
  juv_base_s[1,,q] = rbeta((Nyear+1), juv_s1$alpha, juv_s1$beta)
  juv_s2 = estBetaParams(mu = year1sur, var = year1var)
  juv_base_s[2,,q] = rbeta((Nyear+1), juv_s2$alpha, juv_s2$beta)
  egg_s1 = estBetaParams(mu = eggsur, var = eggvar)
  egg_base_s[,q] = rbeta((Nyear+1), egg_s1$alpha, egg_s1$beta)
  chick_s1 = estBetaParams(mu = chicksur, var = chickvar)
  chick_base_s[,q] = rbeta((Nyear+1), chick_s1$alpha, chick_s1$beta)
}

adult_base_s = adult_base_s^(1/2)
juv_base_s = juv_base_s^(1/2)

#prey availability
# more constraint - high variance
set.seed(821) #set. seed was 821
props = array(NA, c(2,Nyear+1,nsims))
nonbreedprops =  array(NA, c(2,Nyear+1,nsims))
for(k in 1:nsims) {
  props[1,,k] = rlnorm(Nyear+1, meanlog = log(1) - ((0.2^2)/2), sdlog = 0.2)
  nonbreedprops[1,,k] = rlnorm(Nyear+1, log(1) - ((0.05^2)/2),  0.05)
  props[2,,k] = rlnorm(Nyear+1, log(1) - ((0.05^2)/2), 0.05)
  nonbreedprops[2,,k] = props[2,,k]
}

#less constraint - less variance
set.seed(821) #set. seed was 821
props_low = array(NA, c(2,Nyear+1,nsims))
nonbreedprops_low =  array(NA, c(2,Nyear+1,nsims))
for(k in 1:nsims) {
  props_low[1,,k] = rlnorm(Nyear+1, log(1) - ((0.1^2)/2), 0.1)
  nonbreedprops_low[1,,k] = rlnorm(Nyear+1, log(1) - ((0.025^2)/2), 0.025)
  props_low[2,,k] = rlnorm(Nyear+1, log(1) - ((0.025^2)/2), 0.025)
  nonbreedprops_low[2,,k] = props_low[2,,k]
}

#fish = nofish$penguinfood
fish = nofish$avgbiomass
#fish = nofish$biomass.oneplus.true
B0init_low = rep(NA, length = nsims)
B02init_low = rep(NA, length = nsims)
P0noninit_low = rep(NA, length = nsims)

B0init_h = rep(NA, length = nsims)
B02init_h = rep(NA, length = nsims)
P0noninit_h = rep(NA, length = nsims)
# 
# for(j in 1:nsims) {
#   B0init_low[j] = mean(fish[1,201:1001,j]*props_low[1,201:1001,j])
#   B02init_low[j] = mean(fish[2,201:1001,j]*props_low[2,201:1001,j])
#   P0noninit_low[j] = mean(fish[1,201:1001,j]*nonbreedprops_low[1,201:1001,j])
#   B0init_h[j] = mean(fish[1,201:1001,j]*props[1,201:1001,j])
#   B02init_h[j] = mean(fish[2,201:1001,j]*props[2,201:1001,j])
#   P0noninit_h[j] = mean(fish[1,201:1001,j]*nonbreedprops[1,201:1001,j])
# }

# average biomass in each year 
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

seabirdout = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                  recruits = matrix(nrow = nsims, ncol = Nyear),
                  totalpop = matrix(nrow = nsims, ncol = Nyear),
                  fledglings = matrix(nrow = nsims, ncol = Nyear),
                  eggs = matrix(nrow = nsims, ncol = Nyear),
                  testing = matrix(nrow = nsims, ncol = Nyear),
                  testing2 = matrix(nrow = nsims, ncol = Nyear),
                  testing3 = matrix(nrow = nsims, ncol = Nyear))

# still need different biomass in first half of year vs. second half of year
# because of distance birds can travel. But base fish amount will be same
# in first half as second half 
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
  firsthalf = fish[ss,] # base first half and 2nd half should be the same now 
  secondhalf = fish[ss,]
  #print(c(firsthalf, secondhalf))
  
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
  seabird = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[1],
                                            K = K, maxage = maxage[1], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                            clutch = clutch[1],P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
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
  
  seabirdout[["breeders"]][ss,] = colSums(seabird$pop[2,5:31,])
  seabirdout[["recruits"]][ss,] = seabird$recruits
  seabirdout[["totalpop"]][ss,] = colSums(seabird$pop[2,1:31,])
  seabirdout[["fledglings"]][ss,] = seabird$pop[2,1,]
  seabirdout[["eggs"]][ss,] = seabird$pop[1,1,]
  seabirdout[["testing"]][ss,] = seabird$testing
  seabirdout[["testing2"]][ss,] = seabird$testing2
  seabirdout[["testing3"]][ss,] = seabird$testing3
}

plot(colMeans(seabirdout$totalpop))
for(i in 1:nsims) {
    if(any(is.na(seabirdout$totalpop[i,])))
      print(i)
}
minbiomass = rep(NA, nsims)
for(i in 1:nsims) {
  minbiomass[i]=(min(nofish$biomass.total.true[i,]))
}
# certainty
# x.seq = 1:1000
# plot(1:1000, seabirdout$totalpop[1,], type = 'l')
# lines(1:1000, seabirdout$totalpop[2,], col = "Red")
# mean1 = colMeans(seabirdout$totalpop)
# CI1 = apply(seabirdout$totalpop, 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
# plot(x.seq, mean1, type = 'n')
# #lines(x.seq, CI1[1,], lty = 2)
# #lines(x.seq, CI1[2,], lty = 2)
# #polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
# polygon(c(x.seq, rev(x.seq)), c(CI1[1,], rev(CI1[2,])), col = 'grey80', border = NA)
# lines(x.seq, mean1)

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

seabirdout2 = list(breeders = matrix(nrow = nsims, ncol = Nyear), # at end of year
                   recruits = matrix(nrow = nsims, ncol = Nyear),
                   totalpop = matrix(nrow = nsims, ncol = Nyear),
                   fledglings = matrix(nrow = nsims, ncol = Nyear),
                   eggs = matrix(nrow = nsims, ncol = Nyear),
                   testing = matrix(nrow = nsims, ncol = Nyear),
                   testing2 = matrix(nrow = nsims, ncol = Nyear),
                   testing3 = matrix(nrow = nsims, ncol = Nyear))

ffbiomass = array(NA, c(2,Nyear+1,nsims))
ffbiomass_non = array(NA, c(2,Nyear+1,nsims))

#B0init2 = rep(NA, length = nsims)
#B02init2 = rep(NA, length = nsims)
#P0noninit2 = rep(NA, length = nsims)

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
  
  #B0init2[ss] = mean(ffbiomass[1,,ss])
  #B02init2[ss] = mean(ffbiomass[2,,ss])
  #P0noninit2[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  # generalist 
  seabird2 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[2],
                                             K = K, maxage = maxage[2], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[2], P0 = B0init_low_mean, B0 = B0init_low_mean, P02 = B02init_low_mean, P0non = P0noninit_low_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
                                             param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
                                             param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdout2[["breeders"]][ss,] = colSums(seabird2$pop[2,(agebreeding[2]+1):(maxage[2]+1),])
  seabirdout2[["recruits"]][ss,] = seabird2$recruits
  seabirdout2[["totalpop"]][ss,] = colSums(seabird2$pop[2,1:(maxage[2]+1),])
  seabirdout2[["fledglings"]][ss,] = seabird2$pop[2,1,]
  seabirdout2[["eggs"]][ss,] = seabird2$pop[1,1,]
  seabirdout2[["testing"]][ss,] = seabird2$testing
  seabirdout2[["testing2"]][ss,] = seabird2$testing2
  seabirdout2[["testing3"]][ss,] = seabird2$testing3
}


########################## 3 ###################################
maxage_use = maxage[3]
testinits = rep(10000, length = (maxage_use + 1))
K_init = sum(testinits[2:(maxage_use+1)])
stable_dis = SeabirdPop_dd_nofish_nostoch(Nyear = Nyear, inits = testinits, juv_s = year1sur, adult_s = adultsur, 
                                          agebreeding = agebreeding[3], 
                                          K = K_init, 
                                          maxage = maxage[3], egg_s = eggsur, chick_s = chicksur, clutch = clutch[3])
total = sum(stable_dis$pop[2,1:(maxage_use+1),Nyear+1])
stable_dis$pop[2,1:(maxage_use+1),Nyear+1]/total
#plot(colSums(stable_dis$pop[2,1:(maxage_use+1),]))#, ylim = c(0,350000))
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

#B0init3 = rep(NA, length = nsims)
#B02init3 = rep(NA, length = nsims)
#P0noninit3 = rep(NA, length = nsims)
set.seed(216)
for(ss in 1:nsims) {
  #firsthalf = fish[1,,ss] 
  #secondhalf = fish[2,,ss]
  firsthalf = fish[ss,] # base first half and 2nd half should be the same now 
  secondhalf = fish[ss,]
  ffbiomass[1,,ss] = firsthalf*props_low[1,,ss]
  ffbiomass[2,,ss] = secondhalf*props_low[2,,ss]
  ffbiomass_non[1,,ss] = firsthalf*nonbreedprops_low[1,,ss]
  ffbiomass_non[2,,ss] = secondhalf*nonbreedprops_low[2,,ss]
  
  #B0init3[ss] = mean(ffbiomass[1,,ss])
  #B02init3[ss] = mean(ffbiomass[2,,ss])
  #P0noninit3[ss] = mean(ffbiomass_non[1,,ss])
  
  ffbiomass_use = ffbiomass[,,ss]
  ffbiomass_non_use = ffbiomass_non[,,ss]
  
  #set.seed(821)
  seabird3 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
                                             juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[3],
                                             K = K, maxage = maxage[3], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
                                             clutch = clutch[3], P0 = B0init_h_mean, B0 = B0init_h_mean, P02 = B02init_h_mean, P0non = P0noninit_h_mean, 
                                             ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
                                             ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
                                             param_breeders = param_breeders_mid, param_fledge = param_fledge_mid, 
                                             param_fledge2 = param_fledge2_mid, param_adult = param_adult_mid, 
                                             param_juv = param_juv_mid, param_fledge3 = param_fledge3_mid)
  
  # seabird3 = SeabirdPop_simple_stochastic_RS(Nyear = Nyear, inits = testinits, ffbiomass = ffbiomass_use,
  #                                            juv_sur = juv_base_s[,,ss], adult_sur = adult_base_s[,,ss], agebreeding = agebreeding[3],
  #                                            K = K, maxage = maxage[3], egg_sur = egg_base_s[,ss], chick_sur = chick_base_s[,ss],
  #                                            clutch = clutch[3], P0 = B0init3[ss], B0 = B0init3[ss], P02 = B02init3[ss], P0non = P0noninit3[ss], 
  #                                            ffbio_nonbreeders = ffbiomass_non_use, scenario = 1,
  #                                            ffbio_init = ffbiomass[1,,ss], ffbio_init_non = ffbiomass_non[1,,ss], fscen = 1, RS_thres = 0, 
  #                                            param_breeders = param_breeders_general, param_fledge = param_fledge_general, 
  #                                            param_fledge2 = param_fledge2_general, param_adult = param_adult_general, 
  #                                            param_juv = param_juv_general, param_fledge3 = param_fledge3_general)
  seabirdout3[["breeders"]][ss,] = colSums(seabird3$pop[2,(agebreeding[3]+1):(maxage[3]+1),])
  seabirdout3[["recruits"]][ss,] = seabird3$recruits
  seabirdout3[["totalpop"]][ss,] = colSums(seabird3$pop[2,1:(maxage[3]+1),])
  seabirdout3[["fledglings"]][ss,] = seabird3$pop[2,1,]
  seabirdout3[["eggs"]][ss,] = seabird3$pop[1,1,]
  seabirdout3[["testing"]] = seabird3$testing
}



#### OLD PLOTS, OLD CODE - NOT USED #########
plot(1:1000, seabirdout$totalpop[1,], ylim = c(0, 500000), type = 'l')
lines(1:1000, seabirdout2$totalpop[1,], col = "Red")
lines(1:1000, seabirdout3$totalpop[1,], col = "Blue")


par(mar = c(4,4,3,2))
x.seq = 1:1000
mean1 = colMeans(seabirdout$totalpop)
CI1 = apply(seabirdout$totalpop, 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
plot(x.seq, mean1, type = 'n', xlab = "Years", ylab = "Population",  ylim = c(0,5500000), bty = "n")
polygon(c(x.seq, rev(x.seq)), c(CI1[1,], rev(CI1[2,])), col = rgb(0,0,0,0.3), border = NA)
lines(x.seq, mean1)

legend("bottomleft", legend = c("Scen 1", "Scen 2", "Scen 3", "Scen 4", "Scen 5", "Scen 6"),
       col=c("black", "red", "blue", "darkgreen", "goldenrod", "purple"), lty=1, cex=0.8, ncol = 3, bty = 'n')

mean2 = colMeans(seabirdout2$totalpop)
CI2 = apply(seabirdout2$totalpop, 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
rednew = rgb(100,0,0,alpha = 30, maxColorValue = 100)
polygon(c(x.seq, rev(x.seq)), c(CI2[1,], rev(CI2[2,])), col = rednew, border = NA)
lines(x.seq, mean2, col = "Red")

mean3 = colMeans(seabirdout3$totalpop)
CI3 = apply(seabirdout3$totalpop, 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})
bluenew = rgb(0,0,100,alpha = 30, maxColorValue = 100)
polygon(c(x.seq, rev(x.seq)), c(CI3[1,], rev(CI3[2,])), col = bluenew, border = NA)
lines(x.seq, mean3, col = "Blue")

