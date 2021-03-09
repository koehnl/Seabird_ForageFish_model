# Creates Figure 5 from Koehn et al. 2021 - Sensitivity analysis
# Runs sensitivity analysis for seabird model parameters
# both functional response parameter senstivitiy analysis and life history
# parameter sensitivity analysis (age at breeding, eggs laid, etc.)
# Runs for restricted seabird only as this seabird scenario was most sensitive to fishing
# runs with both sardine and anchovy prey

# warnings may occur and that's due to seabird models that crash (rbinom calls in the
# seabird model with NA)
#are set to 0 in the code.

library(here)
setwd(here::here())

foragefish =c("Anchovy", "Sardine")
for(i in 1:length(foragefish)) {
  fishuse = foragefish[i]
  # load forage fish runs
  load(file=paste(fishuse,1,"nofish",".RData",sep="_"))
  load(file=paste(fishuse,2,"constF",".RData",sep="_"))
  # run/source seabird runs - both with non-fish prey 
  source(file.path("Run_seabird_2lifehistory.R"))
  
  # Functional response sensitivity
  source(file.path("Fig5 code/VaryFunctionalResponse2.R"))
  
  meanAX = rep(0, length = nsims)
  meanB0X = rep(0, length = nsims)
  meanCX = rep(0, length = nsims)
  meanDX = rep(0, length = nsims)
  meanEX = rep(0, length = nsims)
  #meanFX = rep(0, length = nsims)
  
  for(i in 1:nsims) {
    meanAX[i] = mean(seabirdoutB$totalpop[i,201:1000]/seabirdout$totalpop[i,201:1000])
    meanB0X[i] = mean(seabirdoutB2$totalpop[i,201:1000]/seabirdout2$totalpop[i,201:1000])
    meanCX[i] = mean(seabirdoutB3$totalpop[i,201:1000]/seabirdout3$totalpop[i,201:1000])
    meanDX[i] = mean(seabirdoutB4$totalpop[i,201:1000]/seabirdout4$totalpop[i,201:1000])
    meanEX[i] = mean(seabirdoutB5$totalpop[i,201:1000]/seabirdout5$totalpop[i,201:1000])
    #meanFX[i] = mean(seabirdoutB6$totalpop[i,201:1000]/seabirdout6$totalpop[i,201:1000])
  }
  meanAX[is.na(meanAX)] = 0
  meanB0X[is.na(meanB0X)] = 0
  meanCX[is.na(meanCX)] = 0
  meanDX[is.na(meanDX)] = 0
  meanEX[is.na(meanEX)] = 0
  
  assign(paste(tolower(fishuse),"A",sep=""),  meanAX)
  assign(paste(tolower(fishuse),"B",sep=""),  meanB0X)
  assign(paste(tolower(fishuse),"C",sep=""),  meanCX)
  assign(paste(tolower(fishuse),"D",sep=""),  meanDX)
  assign(paste(tolower(fishuse),"E",sep=""),  meanEX)
  
  # Other life history parameter sensitivity
  source(file.path("Fig5 code/Scenario1_lifehistory.R"))
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
  
  assign(paste(tolower(fishuse),"A2",sep=""),  meanA)
  assign(paste(tolower(fishuse),"B2",sep=""),  meanB0)
  assign(paste(tolower(fishuse),"C2",sep=""),  meanC)
  assign(paste(tolower(fishuse),"D2",sep=""),  meanD)
  assign(paste(tolower(fishuse),"E2",sep=""),  meanE)
  assign(paste(tolower(fishuse),"F2",sep=""),  meanF)
  
}

# results above include sensitivity analysis for parameter max age but not
# presented in paper (unsure if a restricted seabird life history with lower max age exists
# usually long-lived species). Not presented here because this creates the figure from the paper
meanmed = c(quantile(sardineA2, 0.5), quantile(anchovyA2, 0.5), quantile(sardineB2, 0.5), quantile(anchovyB2, 0.5),
            quantile(sardineC2, 0.5), quantile(anchovyC2, 0.5), quantile(sardineD2, 0.5), quantile(anchovyD2, 0.5),
            quantile(sardineE2, 0.5), quantile(anchovyE2, 0.5))#, quantile(sardineF2, 0.5), quantile(anchovyF2, 0.5))
meanlow = c(quantile(sardineA2, 0.025), quantile(anchovyA2, 0.025), quantile(sardineB2, 0.025), quantile(anchovyB2, 0.025),
              quantile(sardineC2, 0.025), quantile(anchovyC2, 0.025), quantile(sardineD2, 0.025), quantile(anchovyD2, 0.025),
              quantile(sardineE2, 0.025), quantile(anchovyE2, 0.025))#, quantile(sardineF2, 0.025), quantile(anchovyF2, 0.025))
meanhigh = c(quantile(sardineA2, 0.975), quantile(anchovyA2, 0.975), quantile(sardineB2, 0.975), quantile(anchovyB2, 0.975),
             quantile(sardineC2, 0.975), quantile(anchovyC2, 0.975), quantile(sardineD2, 0.975), quantile(anchovyD2, 0.975),
             quantile(sardineE2, 0.975), quantile(anchovyE2, 0.975))#, quantile(sardineF2, 0.975), quantile(anchovyF2, 0.975))
meanlow2 = c(quantile(sardineA2, 0.25), quantile(anchovyA2, 0.25), quantile(sardineB2, 0.25), quantile(anchovyB2, 0.25),
             quantile(sardineC2, 0.25), quantile(anchovyC2, 0.25), quantile(sardineD2, 0.25), quantile(anchovyD2, 0.25),
             quantile(sardineE2, 0.25), quantile(anchovyE2, 0.25))#, quantile(sardineF2, 0.25), quantile(anchovyF2, 0.25))
meanhigh2 = c(quantile(sardineA2, 0.75), quantile(anchovyA2, 0.75), quantile(sardineB2, 0.75), quantile(anchovyB2, 0.75),
              quantile(sardineC2, 0.75), quantile(anchovyC2, 0.75), quantile(sardineD2, 0.75), quantile(anchovyD2, 0.75),
              quantile(sardineE2, 0.75), quantile(anchovyE2, 0.75))#, quantile(sardineF2, 0.75), quantile(anchovyF2, 0.75))


meanmed2 = c(quantile(sardineA, 0.5), quantile(anchovyA, 0.5), quantile(sardineB, 0.5), quantile(anchovyB, 0.5),
             quantile(sardineC, 0.5), quantile(anchovyC, 0.5), quantile(sardineD, 0.5), quantile(anchovyD, 0.5),
             quantile(sardineE, 0.5), quantile(anchovyE, 0.5))
meanlow2B = c(quantile(sardineA, 0.025), quantile(anchovyA, 0.025), quantile(sardineB, 0.025), quantile(anchovyB, 0.025),
              quantile(sardineC, 0.025), quantile(anchovyC, 0.025), quantile(sardineD, 0.025), quantile(anchovyD, 0.025),
              quantile(sardineE, 0.025), quantile(anchovyE, 0.025))
meanhigh2B = c(quantile(sardineA, 0.975), quantile(anchovyA, 0.975), quantile(sardineB, 0.975), quantile(anchovyB, 0.975),
               quantile(sardineC, 0.975), quantile(anchovyC, 0.975), quantile(sardineD, 0.975), quantile(anchovyD, 0.975),
               quantile(sardineE, 0.975), quantile(anchovyE, 0.975))
meanlow22 = c(quantile(sardineA, 0.25), quantile(anchovyA, 0.25), quantile(sardineB, 0.25), quantile(anchovyB, 0.25),
              quantile(sardineC, 0.25), quantile(anchovyC, 0.25), quantile(sardineD, 0.25), quantile(anchovyD, 0.25),
              quantile(sardineE, 0.25), quantile(anchovyE, 0.25))
meanhigh22 = c(quantile(sardineA, 0.75), quantile(anchovyA, 0.75), quantile(sardineB, 0.75), quantile(anchovyB, 0.75),
               quantile(sardineC, 0.75), quantile(anchovyC, 0.75), quantile(sardineD, 0.75), quantile(anchovyD, 0.75),
               quantile(sardineE, 0.75), quantile(anchovyE, 0.75))
##### BOTH LIFE HISTORY AND FUNCTIONAL RESPONSE SENSITIVITY ######
# which order plotted, anchovy vs. sardine, would need to match with legend
#par(mfrow = c(1,2))
names = c("Lower Max Age", "Base scenario", "Lower variance", "Lower age at\n first breeding", "Larger clutch",
          "All functional \nresponse generalist", "Reproductive success \ngeneralist",
          "Juvenile survival \ngeneralist","Breeder attendance\n generalist",
          "Adult survival \ngeneralist")
x = c(0,0.2,1,1.2,2,2.2,3,3.2,4,4.2, 5,5.2,6,6.2,7,7.2,8,8.2)
#par(mfrow = c(1,1))
par(mar = c(6,11,1,2))
plot(c(meanmed2[c(3,4,9,10,5,6,7,8)],meanmed[c(7,8,5,6,9,10,3,4,1,2)]), x, frame = F, xlab = "", ylab = "", xlim = c(0,1), axes = F, pch = 19, 
     col = c("black", "grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))

#text(meanmed[c(7,8,5,6,9,10,3,4,1,2)], x, labels = meanmed[c(7,8,5,6,9,10,3,4,1,2)] )

arrows(c(meanlow2B[c(3,4,9,10,5,6,7,8)], meanlow[c(7,8,5,6,9,10,3,4,1,2)]), x, 
       c(meanhigh2B[c(3,4,9,10,5,6,7,8)],meanhigh[c(7,8,5,6,9,10,3,4,1,2)]), x, length=0, angle=90, code=3, lwd = 1,
       col = c("black", "grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
arrows(c(meanlow22[c(3,4,9,10,5,6,7,8)],meanlow2[c(7,8,5,6,9,10,3,4,1,2)]), x,
       c(meanhigh22[c(3,4,9,10,5,6,7,8)],meanhigh2[c(7,8,5,6,9,10,3,4,1,2)]), x, length=0, angle=90, code=3, lwd = 4, 
       col = c("black", "grey50", "black","grey50", "black","grey50", "black","grey50", "black","grey50"))
axis(2, at = c(0.1,1.1,2.1,3.1,4.1, 5.1,6.1,7.1,8.1,9.1), labels = c(names[10:1]), las = 2, pos = -0.05, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = -0.2)
mtext("Average relative seabird \nabundance",side = 1, line = 3.5, cex = 1)
#mtext("(B) Sardine Prey", side = 3, at = c(-0.3,10))
legend(-0.5,-1.4,xpd=NA, legend = c("Anchovy prey","Sardine prey"), col = c("grey50", "black"), lty = 1, pch = 1, lwd = 4)
#text(-0.6,4.5, label = "(A)", xpd = NA, cex = 1.5)
abline(h = 3.5, lty = 2)



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
