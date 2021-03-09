# Creates Figure 8 Koehn et al. 2021
# loads forage fish runs for each forage fish (sardine and anchovy)
# and generates triangle plots to compare catch, variance in catch, and 
# probability of collapse across forge fish harvest control rules
# for each index (catch, variance, prob collapse) scales ranking of each harvest
# control rule between 0 and 3 (where 3 is "best" - most catch, low variance, low collapse prob.)

library(matrixStats)

foragefish =c("Anchovy", "Sardine")
for(i in 1:length(foragefish)) {
  fishuse = foragefish[i]
  load(file=paste(fishuse,1,"nofish",".RData",sep="_"))
  load(file=paste(fishuse,2,"constF",".RData",sep="_"))
  load(file=paste(fishuse,3,"constFhi",".RData",sep="_"))
  load(file=paste(fishuse,5,"C1",".RData",sep="_"))
  load(file=paste(fishuse,6,"C2",".RData",sep="_"))
  load(file=paste(fishuse,7,"C3",".RData",sep="_"))
  nsims = nrow(nofish$fishing)
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
  load(file=paste(fishuse,"equilib",".RData",sep="_"))
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
  
  # FISH only
  resultslow = results[-1,]
  
  
  library(scales)
  #resultslow$birdmean = rescale(resultslow$birdmean, to = c(0,3))
  #resultslow$birdperlow = rescale(resultslow$birdperlow, to = c(0,3))
  resultslow$meancatch = rescale(resultslow$meancatch, to = c(0,3))
  resultslow$sdcatch = rescale(resultslow$sdcatch, to = c(3,0))
  resultslow$meanbiomass = rescale(resultslow$meanbiomass, to = c(0,3))
  resultslow$sdbiomass = rescale(resultslow$sdbiomass, to = c(3,0))
  resultslow$Probcollapse = rescale(resultslow$Probcollapse, to = c(3,0))
  resultslow$years0catch = rescale(resultslow$years0catch, to = c(3,0))
  
  assign(paste(tolower(fishuse),"fishonly",sep=""),  resultslow)
  
}
# save for each fish you run
#sardinefishonly = resultslow
#anchovyfishonly = resultslow
sardinefishonly = as.data.frame(rbind(rep(3, 3), rep(0,3),sardinefishonly[,c(4,5,9)]))
anchovyfishonly = as.data.frame(rbind(rep(3, 3), rep(0,3),anchovyfishonly[,c(4,5,9)]))

# Make triangle plot graphs
library(fmsb)
par(mfrow = c(2,5))
par(oma = c(1,8,1,1))
par(mar = c(2,2,2,2))
par(mgp = c(4,1,0))
par(xpd = TRUE)
hcrnew = c("Constant low", "Constant moderate", "Hockey moderate",
           "Hockey low cut-off", "Hockey high F max")
for(i in 1:5) {
  testing = as.data.frame(rbind(sardinefishonly[1:2,],sardinefishonly[i+2,]))
  
  radarchart( testing , axistype=1 , 
              #custom polygon
              pcol=rgb(0.4,0.4,0.4,1) , pfcol=rgb(0,0,0,0)  , plwd=4 , 
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,3,1), cglwd=0.8,
              vlabels = c("Mean \ncatch",  "min. SD catch", "Min. years \n0 catch"),
              title = hcrnew[i],
              #custom labels
              vlcex=0.9, seg = 3
  )
}

text(x = -20, y = 0, "Sardine", cex = 1.25,xpd = NA)

for(i in 1:5) {
  #testing = as.data.frame(rbind(sardinelow[1:2,c(2,3,6,8)], sardinelow[i+2,c(2,3,6,8)]))
  testing2 = as.data.frame(rbind(anchovyfishonly[1:2,], anchovyfishonly[i+2,]))
  
  radarchart( testing2 , axistype=1 , 
              #custom polygon
              pcol=rgb(0.5,0.5,0.5,0.9) , pfcol=rgb(0,0,0,0) , plwd=4 , 
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,3,1), cglwd=0.8,
              vlabels = c(
                "Mean \ncatch",  "min SD catch", "Min. years \n0 catch"),
              #title = hcrnew[i],
              #custom labels
              vlcex=0.9, seg = 3 
  )
}
text(x = -20, y = 0, "Anchovy", cex = 1.25,xpd = NA)
