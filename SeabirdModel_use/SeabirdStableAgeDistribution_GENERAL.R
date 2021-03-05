vgbf_func <- function(E.inf, adults, clutch) {
  if(clutch == 2) {
    alpha = 1/E.inf
  } else if(clutch == 1) {
    alpha = 1/(2*E.inf)
  } else {
    alpha = 2/(3*E.inf) # I think this works for 3 egg clutches
  }
  eggs = E.inf*(1-exp(-alpha*adults))
  return(eggs)
}
#adults = seq(100,30000, by = 100)
#clutch = 3
#K = 10000
#x = vgbf_func(E.inf = K, adults = adults, clutch = clutch)

SeabirdPop_dd_nofish_nostoch <- function(Nyear, inits, juv_s, adult_s, agebreeding, K, 
                                         maxage, egg_s, chick_s, clutch) { 
  #Set up age structure
  # REMEMBER! age 0 is index 1, so even though breeding age is 4 yrs old, this is at index 5
  breeders = agebreeding + 1 # breeding group index
  BirdN = array(NA, c(2,maxage+1, Nyear+1)) # 2 cause half a year time steps 
  recruits = rep(NA, length = Nyear + 1)
  #Set up initial population sizes - N[1,age, 1]
  
  #available = foragetype$dist*foragetype$depth*ffbiomass[1,1] # prey available to predator based on 
  # some how translate depth and dist to percentags available?
  #ws = diet_pref/B0
  #Py = otherprey + available*ws #other prey is a proportion...was in my model for 558
   
  for(j in 2:(maxage+1)) {
    BirdN[1,j,1] = inits[j]
  }

  Nmat = sum(BirdN[1,breeders:(maxage+1),1])
  recruits[1] = BirdN[1,breeders,1]
  Females = Nmat/2
  BirdN[1,1,1] = vgbf_func(E.inf = K, adults = Females*2, clutch = clutch)
  #survival = c(survival[1], survival[2:(length(survival))]^(1/2))
  
  for(i in 1:(Nyear+1)) {     
    for(k in breeders:(maxage+1)) {
      BirdN[2,k,i] = (adult_s^(1/2))*BirdN[1,k, i]#*S_adult
    }
        
    youngtemp = BirdN[1,1,i]*egg_s#*phi_egg
    BirdN[2,1,i] = youngtemp*chick_s#chick_s2^(1/2)#*phi_y
    
#   S_juv = S_juv^(1/2)
    for(l in 2:(breeders-1)) { # at the least, breeders could = 3 (so breeding age at 2)
      if(l == 2) {
        BirdN[2,l,i] = (juv_s^(1/2))*BirdN[1,l,i]
      } else {
        BirdN[2,l,i] = (adult_s^(1/2))*BirdN[1,l,i]
      }
    }
    ############################################################################
    if(i < (Nyear+1)) {
      recruits[i+1] = (adult_s^(1/2))*BirdN[2,breeders-1,i]#*S_adult # used in top-down indicator HCR
      BirdN[1,breeders,i+1] = recruits[i+1]
      
      for(k in breeders:(maxage)) { # those in maxage + 1 die, so not surviving to i+1
        BirdN[1,k+1,i+1] = (adult_s^(1/2))*BirdN[2,k, i]#*S_adult
      }
      # those at max age die afterwards
      
      for(k in 1:(breeders-2)) { # I think this works for any breeding age 
        if(k == 1) {
          BirdN[1,k+1,i+1] = (juv_s^(1/2))*BirdN[2,k,i]#*S_juv
        } else { # juvenile at the beginning of the year has adult survival in 2nd half of year
          # because has had "juvenile" survival for a full year (fledged in 2nd half of first year)
          BirdN[1,k+1,i+1] = (adult_s^(1/2))*BirdN[2,k,i]#*S_adult
        }
      }
      
      Nmat2 = (BirdN[1, breeders, i+1]) + sum(BirdN[1,(breeders+1):(maxage+1),i+1])
      Females2 = Nmat2/2
      BirdN[1,1,i+1] = vgbf_func(E.inf = K, adults = Females2*2, clutch = clutch)
    }
  }
  
  return(list("pop" = BirdN,  "recruits" = recruits))
}

# stable distribution is the same with or without impacts of fish
# if fish are constant. 