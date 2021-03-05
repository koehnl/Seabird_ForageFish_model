##########################################################
####### Stochastic Model #################################
##########################################################
##########################################################
####### Use for actual runs ##############################
##########################################################
#####################################
##########################################################

vgbf_func <- function(E.inf, adults, clutch) {
  if(clutch == 2) {
    alpha = 2/(2*E.inf)
  } else if(clutch == 1) {
    alpha = 1/(2*E.inf)
  } else {
    alpha = 3/(2*E.inf) # I think this works for 3 egg clutches
  }
  eggs = E.inf*(1-exp(-alpha*adults))
  return(eggs)
}
adults = seq(0,10000, by = 10)
test = vgbf_func(10000,adults,3)
par(mar = c(5,5,1,1))
plot(adults, test, type = "l", ylim = c(0,10000), xlab = "Adults", ylab = "# Viable Eggs")
test2 = vgbf_func(2000, adults, 3)
lines(test2)
SeabirdPop_simple_stochastic_RS <- function(Nyear, inits, ffbiomass, juv_sur,
                                            adult_sur, agebreeding, ffbiomass_non,
                                            K, maxage, egg_sur, chick_sur, clutch,
                                            P0, B0, P02, P0non, 
                                            ffbio_nonbreeders, scenario,
                                            ffbio_init, ffbio_init_non, fscen, RS_thres,
                                            param_breeders, param_fledge, param_fledge2, param_adult, 
                                            param_juv, param_fledge3) { 
  #set.seed(216) # needs to be different between sims but same for scenarios
  # should work without this seed here once I set up the run code to loop through scenarios
  #Set up age structure
  # REMEMBER! age 0 is index 1, so even though breeding age is 4 yrs old, this is at index 5
  breeders = agebreeding + 1 # breeding group index
  BirdN = array(NA, c(2,maxage+1, Nyear)) # 2 cause half a year time steps
  Pymat = matrix(nrow = 2, ncol = Nyear) 
  recruits = rep(NA, length = Nyear); recruits[1] = inits[breeders]
  breeders_contents = rep(NA, length = Nyear)
  Nmatvector = rep(NA, length = Nyear)
  RS = rep(NA, length = Nyear)
  #ws = diet_pref[2,1]/B0
    
  for(j in 2:(maxage+1)) {
    BirdN[1,j,1] = inits[j]
  }
  Nmat = sum(BirdN[1,breeders:(maxage+1),1])
  Nmatvector[1] = Nmat
  available_lay = ffbiomass[2,1] # end of first forage fish year, use to be [6,1 when there were 6 time steps]
  #Py_lay = otherprey[2,1] + available_lay*ws 
  Py_lay = available_lay
  # alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
  alpha = param_breeders[1]; slope = param_breeders[2]; beta = param_breeders[3]
  phi_attend = max(0,min(1, alpha + (1-alpha)/(1.0+exp(-slope*((Py_lay/P02)-beta)))))
  
  Females = (Nmat/2)*phi_attend # does this work? phi_attend should be between 0 and 1?
  BirdN[1,1,1] = vgbf_func(E.inf = K, adults = Females*2, clutch = clutch)
  # THE ABOVE IS # "VIABLE" eggs - eggs laid not outside of optimal habitat
  #survival = c(survival[1], survival[2:(length(survival))]^(1/2))
  breeders_contents[1] = round(BirdN[1,1,1], digits = 0)
  # DOES THIS WORK? ASSUMES EVERYONE LAID EGGS AND THEN SOME (MOST) HAD EGGS THAT HATCHED
  testvec = rep(NA, length = Nyear)
  testvec2 = rep(NA, length = Nyear)
  testvec3 = rep(NA, length = Nyear)
  
  for(i in 1:(Nyear)) {
    adult_s = adult_sur[1,i]
    juv_s = juv_sur[1,i]
    
    #ws = (diet_pref[1,i+1])/B0
    available2 = ffbiomass[1,i+1]#sum(ffbiomass[1:3,i+1]) # all age classes?
    #Py = otherprey[1,i+1] + available2*ws #right?
    Py = available2
    Pymat[1,i] = Py
    available_nonbreeder = ffbio_nonbreeders[1,i+1] 
    #Py2 = otherprey[1, i+1] + available_nonbreeder*ws 
    Py2 = available_nonbreeder
    #responsetest = Py/P0; print(responsetest)
    # LOSING CONTENTS ALL TOGETHER
    egg_s2 = egg_sur[i] 
    chick_s2 = chick_sur[i]
    
    if(fscen > 1) {
      if(i > 1 && RS[i-1] < RS_thres) {
        Py = ffbio_init[i+1] 
        Py2 = ffbio_init_non[i+1]
      }
    }
    
    if(clutch == 1) {
      alpha2 = param_fledge3[1]; slope2 = param_fledge3[2]; beta2 = param_fledge3[3]
      # should be the most lenient functional response because only have 1 chick to feed
      phi_fledge1_w1 = max(0,min(1,alpha2 + (1-alpha2)/(1.0+exp(-slope2*((Py/P0)-beta2)))))
      temp = cbind(Py, phi_fledge1_w1)
      #print(temp)
      testvec[i] = phi_fledge1_w1
      prob_survive_tofledge = chick_s2*phi_fledge1_w1*egg_s2
      
      breeders_year = breeders_contents[i]
      pairs = round(breeders_year/2, digits = 0)
      X = rbinom(1, pairs, 1-prob_survive_tofledge) # want how many pairs lose 1 chick
      breeders_lost = X*2
      BirdN[2,1,i] = max(0,(BirdN[1,1,i] - (X))) # how many chicks left 
    
    } else if(clutch == 2) {
      alpha3 = param_fledge2[1]; slope3 = param_fledge2[2]; beta3 = param_fledge2[3]
      # 2nd most constraining - have 2 chicks to feed 
      phi_fledge1_w2 = max(0,min(1,alpha3 + (1-alpha3)/(1.0+exp(-slope3*((Py/P0)-beta3)))))
      prob_survive_tofledge_w2 = chick_s2*phi_fledge1_w2*egg_s2
      
      alpha4 = param_fledge3[1]; slope4 = param_fledge3[2]; beta4 = param_fledge3[3]
      # least constraining because will use for after 1st chick is lost, prob of fledging 2nd 
      # so only need resources for 1 chick 
      phi_fledge2nd_w2 = max(0,min(1,alpha3 + (1-alpha3)/(1.0+exp(-slope3*((Py/P0)-beta3)))))
      prob_survive_tofledge2nd_w2 = chick_s2*phi_fledge2nd_w2*egg_s2
      
      breeders_year = breeders_contents[i]
      pairs = round(breeders_year/2, digits = 0)
      X = rbinom(1, pairs, 1-prob_survive_tofledge_w2) # want how many pairs had 2, lose 1 chick
      breeders_lost = X*2
      #BirdN[2,1,i] = max(0,(BirdN[1,1,i] - (X))) # how many chicks left 
      prob_lose2 = 1-prob_survive_tofledge2nd_w2 # prob. of losing a 2nd if lost the 1st 
      # = 1 - prob of fledging 2nd if lost the 1st 
      Y = rbinom(1, X, prob_lose2) # how many pairs loose 2nd chick too
      breeders_lost = Y*2 # how many individiual breeders lose 2 chicks
      BirdN[2,1,i] = max(0, (BirdN[1,1,i] - (X+Y))) # how many chicks left 
      
    } else { #clutch == 3
      alpha3 = param_fledge[1]; slope3 = param_fledge[2]; beta3 = param_fledge[3] 
      #should be most constrained - have to feed lots of chicks
      phi_fledge1_w3 = max(0,min(1,alpha3 + (1-alpha3)/(1.0+exp(-slope3*((Py/P0)-beta3)))))
      prob_survive_tofledge1_w3 = chick_s2*phi_fledge1_w3*egg_s2
      
      alpha4 = param_fledge2[1]; slope4 = param_fledge2[2]; beta4 = param_fledge2[3] 
      #down to 2/3 chicks - so less constraint but still some
      phi_fledge2nd_w3 = max(0,min(1,alpha3 + (1-alpha3)/(1.0+exp(-slope3*((Py/P0)-beta3)))))
      prob_survive_tofledge2nd_w3 = chick_s2*phi_fledge2nd_w3*egg_s2
      # prob you fledge a 2nd chick, if you lost the 1st, and still have a 3rd 
      
      alphatemp = param_fledge3[1]; slopetemp = param_fledge3[2]; betatemp = param_fledge3[3]
      # least constrainting - 1 chick left 
      phi_fledge3rd_w3 = max(0,min(1,alphatemp + (1-alphatemp)/(1.0+exp(-slopetemp*((Py/P0)-betatemp)))))
      prob_survive_tofledge3rd_w3 = chick_s2*phi_fledge3rd_w3*egg_s2 
      # prob you fledge a 3rd chick if lost the 1st 2 
      
      breeders_year = breeders_contents[i]
      pairs = round(breeders_year/2, digits = 0)
      X = rbinom(1, pairs, 1-prob_survive_tofledge1_w3) # want how many pairs lose 1 chick of 3 
      breeders_lost = X*2
      #BirdN[2,1,i] = max(0,(BirdN[1,1,i] - (X))) # how many chicks left 
      prob_lose2 = 1-prob_survive_tofledge2nd_w3 # prob lose 2nd chick if lost 1st and have a 3rd
      # which equals = 1 - prob of fledging a 2nd when you lost the 1st and have a 3rd
      Y = rbinom(1, X, prob_lose2) # how many pairs loose 2nd chick, of 3, 1 left after 
      breeders_lost = Y*2 # how many individiual breeders lose 2 chicks
      #BirdN[2,1,i] = max(0, (BirdN[1,1,i] - (X+Y))) # how many chicks left 
      prob_lose3 = 1 - prob_survive_tofledge3rd_w3 # how many pairs lose 3rd chick after losing first 2/3
      Z = rbinom(1, Y, prob_lose3)
      breeders_lost = Z*2 # how many individiual breeders lose 3 chicks
      BirdN[2,1,i] = max(0,(BirdN[1,1,i] - (X+Y+Z))) # how many chicks left
    }
    
    # breeders that fledge at least 1 chick
    breeders_contents_year = breeders_year - breeders_lost
    RS[i] = BirdN[2,1,i]
    # all forage fish biomass is i + 1 because biomass in first year determines bird eggs to start
    
    alphaA = param_adult[1]; slopeA = param_adult[2]; betaA = param_adult[3]
    S_adult = max(0,min(1,alphaA + (1-alphaA)/(1.0+exp(-slopeA*((Py/P0)-betaA)))))
    testvec2[i] = S_adult
    S_adult = S_adult^(1/2)
    
    alphanon = param_adult[1]; slopenon = param_adult[2]; betanon = param_adult[3]
    S_adult2 = max(0,min(1,alphanon + (1-alphanon)/(1.0+exp(-slopenon*((Py2/P0non)-betanon)))))
    S_adult2 = S_adult2^(1/2)
    
    # Split up those that lost contents equally among age groups OR
    # Have younger (less experienced) loose contents - with first time breeders losing most
    # Adults that didn't breed - split up equally
    nonbreeders = Nmatvector[i] - (breeders_contents[i]) #never bred to begin with 
    
    Num_breeding_agegroups = (maxage+1) - agebreeding
    percent_eachage = vector(length = Num_breeding_agegroups)
    percent_eachage[1] = (BirdN[1,breeders,i])/(Nmatvector[i])
    percent_eachage[2:Num_breeding_agegroups] = BirdN[1,(breeders+1):(maxage+1),i]/(Nmatvector[i])
    
    num_breed = percent_eachage*breeders_contents_year
    num_non = percent_eachage*nonbreeders
    #print(percent_eachage)
    if(scenario == 1) {
      num_lost = percent_eachage*breeders_lost
      for(k in breeders:(maxage+1)) {
        # Need breeders with contents, breeders that lost contents, and adults that never bred *********
        BirdN[2,k,i] = adult_s*(num_breed[k-agebreeding])*S_adult + 
          adult_s*(num_lost[k-agebreeding]+num_non[k-agebreeding])*S_adult2
        # those that lost contents or that never bred have same survival - assumption
      # die after maxage and 1 survival for adults
      }
          
    } else {
      Num_breeding_agegroups = (maxage+1) - agebreeding
      lostpercents = vector(length = Num_breeding_agegroups)
      #       lostpercents[1] = (2.4*(100/7.2))/100; lostpercents[2] = (2*(100/7.2))/100; lostpercents[3] = (1.8*(100/7.2))/100
      #       lostpercents[4:Num_breeding_agegroups] = ((100/7.2)/100)/(Num_breeding_agegroups - 3)
      lostpercents[1] = 2*percent_eachage[1]; lostpercents[2] = 2*percent_eachage[2]; lostpercents[3] = 1*percent_eachage[3]
      left = 1-(sum(lostpercents[1:3]))
      # break up the rest based on % in each group
      percent_eachage_new = BirdN[1,(breeders+3):(maxage+1),i]/(Nmatvector[i]-((BirdN[1,breeders,i])+BirdN[1, breeders+1,i]+BirdN[1, breeders +2,i]))
      lostpercents[4:Num_breeding_agegroups] = left*percent_eachage_new
      num_lost_young = lostpercents*breeders_lost
      
      for(k in breeders:(maxage+1)) {
        # Need breeders with contents, breeders that lost contents, and adults that never bred *********
        BirdN[2,k,i] = adult_s*(num_breed[k-agebreeding])*S_adult + 
          adult_s*(num_lost_young[k-agebreeding]+num_non[k-agebreeding])*S_adult2
        
      }      
    }
    # make sure just for the first half of the year and not whole year
    
    alphajuv = param_juv[1]; slopejuv = param_juv[2]; betajuv = param_juv[3]
    S_juv = max(0,min(1,alphajuv + (1-alphajuv)/(1.0+exp(-slopejuv*((Py2/P0non)-betajuv)))))
    # should be Py2 because it should be the same amount of prey available to non-breeders
    testvec3[i] = S_juv
    S_juv = S_juv^(1/2)
    
    # older juveniles have less trouble finding prey - same as adults that lost contents 
    for(l in 2:(breeders-1)) {
      if(l == 2) {
        BirdN[2,l,i] = BirdN[1,l,i]*S_juv*juv_s 
      } else {
        BirdN[2,l,i] = adult_s*BirdN[1,l,i]*S_adult2
      }
    }
    
    ############################################################################
    if(i < (Nyear)) {
      adult_s = adult_sur[2,i]
      juv_s = juv_sur[2,i]
      
      #ws = diet_pref[2,i+1]/B0
      available3 = ffbiomass[2,i+1]#sum(ffbiomass[4:6,i+1])
      #Py = otherprey[2, i+1] + available2*ws # right?
      Py = available2
      Pymat[2,i] = Py
      alphaA = param_adult[1]; slopeA = param_adult[2]; betaA = param_adult[3]
      S_adult = max(0,min(1,alphaA + (1-alphaA)/(1.0+exp(-slopeA*((Py/P02)-betaA)))))
      S_adult = S_adult^(1/2)
      
      alphajuv = param_juv[1]; slopejuv = param_juv[2]; betajuv = param_juv[3]
      S_juv = max(0,min(1,alphajuv + (1-alphajuv)/(1.0+exp(-slopejuv*((Py/P02)-betajuv)))))
      S_juv = S_juv^(1/2)
      
      
      recruits[i+1] = adult_s*BirdN[2,breeders-1,i]*S_adult # used in top-down indicator HCR
      BirdN[1,breeders,i+1] = recruits[i+1]
      # only 50% of age 4 breed, all breeding by age 5
      for(k in breeders:(maxage)) {
        BirdN[1,k+1,i+1] = adult_s*BirdN[2,k,i]*S_adult
      } # die after maxage, so those in the maxage spot (maxage + 1), have no survival
       
      
      for(k in 1:(breeders-2)) {
        if(k == 1) {
          BirdN[1,k+1,i+1] = BirdN[2,k,i]*S_juv*juv_s # fledglings surviving to juveniles
        } else {
          BirdN[1,k+1,i+1] = adult_s*BirdN[2,k,i]*S_adult
        }
      }
      
      available_lay = ffbiomass[2,i+1] #use to be [6,i+1]
      #Py_lay = otherprey[2,i+1] + available_lay*ws
      Py_lay = available_lay
      # alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
      alpha = param_breeders[1]; slope = param_breeders[2]; beta = param_breeders[3]
      phi_attend = max(0,min(1,alpha + (1-alpha)/(1.0+exp(-slope*((Py_lay/P02)-beta)))))
      
      Nmatvector[i+1] = (BirdN[1, breeders, i+1]) + sum(BirdN[1,(breeders+1):(maxage+1),i+1]) #Nmat
      Females = ((Nmatvector[i+1])/2)*phi_attend
      BirdN[1,1,i+1] = vgbf_func(E.inf = K, adults = Females*2, clutch = clutch)
      breeders_contents[i+1] = round(BirdN[1,1,i+1], digits = 0)
    }
  }
  
  return(list("pop" = BirdN,  "recruits" = recruits, "RS" = RS,
              "Py" = Pymat, "testing" = testvec, "testing2" = testvec2, "testing3" = testvec3))
}



# OLD CODE
# egg_s2 = egg_sur[i] 
# chick_s2 = chick_sur[i] 
# alpha2 = param_fledge[1]; slope2 = param_fledge[2]; beta2 = param_fledge[3]
# phi_fledge = max(0,min(1,alpha2 + (1-alpha2)/(1.0+exp(-slope2*((Py/P0)-beta2)))))
# temp = cbind(Py, phi_fledge)
# print(temp)
# testvec[i] = phi_fledge 
# prob_survive_tofledge = chick_s2*phi_fledge*egg_s2
# 
# if(clutch == 2 || clutch == 3) {
#   alpha3 = param_fledge2[1]; slope3 = param_fledge2[2]; beta3 = param_fledge2[3]
#   phi_fledge2nd = max(0,min(1,alpha3 + (1-alpha3)/(1.0+exp(-slope3*((Py/P0)-beta3)))))
#   prob_survive_tofledge_2nd = chick_s2*phi_fledge2nd*egg_s2
# }
# 
# if(clutch == 3) {
#   alphatemp = param_fledge3[1]; slopetemp = param_fledge3[2]; betatemp = param_fledge3[3]
#   phi_fledge3rd = max(0,min(1,alphatemp + (1-alphatemp)/(1.0+exp(-slopetemp*((Py/P0)-betatemp)))))
#   prob_survive_tofledge_3rd = chick_s2*phi_fledge3rd*egg_s2 
# }
# 
# breeders_year = breeders_contents[i]
# pairs = round(breeders_year/2, digits = 0)
# X = rbinom(1, pairs, 1-prob_survive_tofledge) # want how many pairs lose 1 chick
# if (clutch == 1) {
#   breeders_lost = X*2
#   BirdN[2,1,i] = max(0,(BirdN[1,1,i] - (X))) # how many chicks left 
# } else { # clutch = 2 or 3
#   prob_lose2 = 1-prob_survive_tofledge_2nd # is this the best?
#   Y = rbinom(1, X, prob_lose2) # how many pairs loose 2 chicks
#   if(clutch == 2) {
#     breeders_lost = Y*2 # how many individiual breeders lose 2 chicks
#     BirdN[2,1,i] = max(0, (BirdN[1,1,i] - (X+Y))) # how many chicks left 
#   } else {
#     prob_lose3 = 1 - prob_survive_tofledge_3rd
#     Z = rbinom(1, Y, prob_lose3)
#     breeders_lost = Z*2 # how many individiual breeders lose 3 chicks
#     BirdN[2,1,i] = max(0,(BirdN[1,1,i] - (X+Y+Z))) # how many chicks left
#   }
# } 
# # breeders that fledge at least 1 chick
# breeders_contents_year = breeders_year - breeders_lost
