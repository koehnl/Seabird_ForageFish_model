# RUN forage fish model with either anchovy or sardine

library(matrixStats)
results <- as.data.frame(matrix(0, nrow = 6, ncol = 9))
colnames(results) = c("hcr", "birdmean", "birdperlow",
                      "meancatch", "sdcatch", "meanbiomass", "sdbiomass" ,"Probcollapse", "years0catch")
HCR = c("nofish","Constant F low","Constant F mod","Hockey Stick", "Hockey Low Blim", "Hockey High F max") # Took out trend because it was unrealistic-- but using trend in CPUE as adjustment (data-poor method) might be a good idea!
results$hcr = HCR


results$meancatch = c(median(rowMeans(nofish$total.catch[,201:1000])), median(rowMeans(constF$total.catch[,201:1000])),
                      median(rowMeans(constFhi$total.catch[,201:1000])), median(rowMeans(C1$total.catch[,201:1000])),
                      median(rowMeans(C2$total.catch[,201:1000])), median(rowMeans(C3$total.catch[,201:1000])))
results$sdcatch = c(median(rowSds(nofish$total.catch[,201:1000])), median(rowSds(constF$total.catch[,201:1000])),
                    median(rowSds(constFhi$total.catch[,201:1000])), median(rowSds(C1$total.catch[,201:1000])),
                    median(rowSds(C2$total.catch[,201:1000])), median(rowSds(C3$total.catch[,201:1000])))
results$meanbiomass = c(median(rowMeans(nofish$biomass.total.true[,201:1000])), median(rowMeans(constF$biomass.total.true[,201:1000])),
                        median(rowMeans(constFhi$biomass.total.true[,201:1000])), median(rowMeans(C1$biomass.total.true[,201:1000])),
                        median(rowMeans(C2$biomass.total.true[,201:1000])),median(rowMeans(C3$biomass.total.true[,201:1000])))
results$sdbiomass = c(median(rowSds(nofish$biomass.total.true[,201:1000])), median(rowSds(constF$biomass.total.true[,201:1000])),
                      median(rowSds(constFhi$biomass.total.true[,201:1000])), median(rowSds(C1$biomass.total.true[,201:1000])),
                      median(rowSds(C2$biomass.total.true[,201:1000])), median(rowSds(C3$biomass.total.true[,201:1000])))
#thres = 0.2*mean(rowMeans(nofish$biomass.total.true[,51:1000]))
thres = 0.2*equilib$B0
temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(nofish$biomass.total.true[i,201:1000] < thres)) 
}
results[1,8] = mean(temp/800)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(constF$biomass.total.true[i,201:1000] < thres)) 
  #print(c(i,min(constF$biomass.total.true[i,201:1000])))
}
results[2,8] =mean(temp/800)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(constFhi$biomass.total.true[i,201:1000] < thres)) 
}
results[3,8] = mean(temp/800)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(C1$biomass.total.true[i,201:1000] < thres)) 
}
results[4,8] = mean(temp/800)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(C2$biomass.total.true[i,201:1000] < thres)) 
}
results[5,8] = mean(temp/800)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(C3$biomass.total.true[i,201:1000] < thres)) 
}
results[6,8] = mean(temp/800)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(C1$total.catch[i,201:1000] == 0)) 
}
results[4,9] = mean(temp)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(C2$total.catch[i,201:1000] == 0)) 
}
results[5,9] = mean(temp)

temp = rep(0, length = nsims)
for(i in 1:nsims) {
  temp[i] = length(which(C3$total.catch[i,201:1000] == 0)) 
}
results[6,9] = mean(temp)

# seabird results - specialist only, put 2 after last letter to get generalist
### RUN  Run_seabird_3types.R
### RUN fishingsimulations_Fig6.R 

#specialist
seabirdoutB$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutC$totalpop[is.na(seabirdoutC$totalpop)] = 0
seabirdoutD$totalpop[is.na(seabirdoutD$totalpop)] = 0
seabirdoutE$totalpop[is.na(seabirdoutE$totalpop)] = 0
seabirdoutF$totalpop[is.na(seabirdoutF$totalpop)] = 0

#generalist
seabirdoutB2$totalpop[is.na(seabirdoutB$totalpop)] = 0
seabirdoutC2$totalpop[is.na(seabirdoutC$totalpop)] = 0
seabirdoutD2$totalpop[is.na(seabirdoutD$totalpop)] = 0
seabirdoutE2$totalpop[is.na(seabirdoutE$totalpop)] = 0
seabirdoutF2$totalpop[is.na(seabirdoutF$totalpop)] = 0

# results for specialist only
results[,2] = c(median(rowMeans(seabirdout$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])), median(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
                median(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])), median(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
                median(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])), median(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])))

birdthres = 0.5*median(rowMeans(seabirdout$totalpop[,201:1000]))

# results[,2] = c(median(rowMeans(seabirdout$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])), median(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
#                 median(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])), median(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
#                 median(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])), median(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])))
# 
# birdthres = 0.5*median(rowMeans(seabirdout$totalpop[,201:1000]))
# results[1,3] = 1
# results[2,3] = median((apply(seabirdoutB$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
# results[3,3] = median((apply(seabirdoutD$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
# results[4,3] =median((apply(seabirdoutC$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
# results[5,3] = median((apply(seabirdoutE$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
# results[6,3] = median((apply(seabirdoutF$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)

results[1,3] = 1
results[2,3] = median((apply(seabirdoutB$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
results[3,3] = median((apply(seabirdoutD$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
results[4,3] =median((apply(seabirdoutC$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
results[5,3] = median((apply(seabirdoutE$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)
results[6,3] = median((apply(seabirdoutF$totalpop[,201:1000], 1, function(x){ length(which(x> birdthres))}))/800)


resultslow = results[-1,]


library(scales)
resultslow$birdmean = rescale(resultslow$birdmean, to = c(0,3))
resultslow$birdperlow = rescale(resultslow$birdperlow, to = c(0,3))
resultslow$meancatch = rescale(resultslow$meancatch, to = c(0,3))
resultslow$sdcatch = rescale(resultslow$sdcatch, to = c(3,0))
resultslow$meanbiomass = rescale(resultslow$meanbiomass, to = c(0,3))
resultslow$sdbiomass = rescale(resultslow$sdbiomass, to = c(3,0))
resultslow$Probcollapse = rescale(resultslow$Probcollapse, to = c(3,0))
resultslow$years0catch = rescale(resultslow$years0catch, to = c(3,0))

# save either sardine or anchovy to respective variable
library(fmsb)
sardinelow = resultslow
sardinelow2 = rbind(rep(3, 9), rep(0,9), sardinelow)
anchovylow = resultslow
anchovylow2 = rbind(rep(3, 9), rep(0,9), anchovylow)


## Scatterplots
#sardineraw = results
library(viridis)
mycolors = viridis(n = 5)
par(xpd = TRUE)
par(mfrow = c(2,2))
par(mar = c(6,2,2,2))
par(oma = c(2,6,4,2))
plot(sardinelow2$Probcollapse[3:7], sardinelow2$birdmean[3:7],col = mycolors, pch = 16,
     ylab = "", xlab = "Sardine Minimize \nProb. of Collapse (scaled)", cex= 3, bty = 'n', xlim = c(-0.2,3.3), ylim =c(-0.3,3.4))
#lines(0:3,0:3, lty = 2)
plot(sardinelow2$meanbiomass[3:7], sardinelow2$birdmean[3:7],col = mycolors, pch = 16,
     ylab = "", xlab = "Sardine Mean Biomass (scaled)", cex=3,  xlim = c(-0.2,3.3), ylim =c(-0.3,3.4), bty = 'n')
#lines(0:3,0:3, lty = 2)

plot(anchovylow2$Probcollapse[3:7], anchovylow2$birdmean[3:7],col = mycolors, pch = 16,
     xlab = "Anchovy Minimize \nProb. of Collapse (scaled)", ylab = "", cex = 3, xlim = c(-0.2,3.3), ylim =c(-0.3,3.4), bty = 'n')
mtext("Scaled Relative Seabird Abundance", side = 2, outer = TRUE, line = 1)
#lines(0:3,0:3, lty = 2)
plot(anchovylow2$meanbiomass[3:7], anchovylow2$birdmean[3:7],col = mycolors, pch = 16,
     xlab = "Anchovy Mean Biomass (scaled)", ylab = "", cex = 3, xlim = c(-0.2,3.3), ylim =c(-0.3,3.4), bty = 'n')
#lines(0:3,0:3, lty = 2)
par(xpd = NA)
legend(-7,16, legend = c("Constant low", "Constant moderate", "Hockey-stick moderate","Hockey low cut-off","Hockey high max fishing"), 
       col = mycolors, pch = 16, horiz = FALSE, ncol = 3,  bty = 'n', pt.cex = 2, cex = 1.2)

