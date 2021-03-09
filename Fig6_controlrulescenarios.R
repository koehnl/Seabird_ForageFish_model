# Creates figure 6 from Koehn et al. 2021 
# seabird relative abundance with fishing scenarios on both sardine and anchovy prey
# under different harvest control rules 
# is run for both the specialist, restricted seabird and for the 
# generalist, flexible seabird

# to make the figure, also sources forage fish model code from Siple et al. 2019
# to make figure that represents visually the harvest control rules

# AND code to show fishing at MSY impacts on the flexible seabird

# NOTE: warnings are OK because some seabird model runs crash and produce
# NAs in rbinom. Those runs are set to 0. 

# loads forage fish runs and sources seabird model code for fishing scenarios
library(matrixStats)

foragefish =c("Anchovy", "Sardine")
for(i in 1:length(foragefish)) {
  fishuse = foragefish[i]
  # load forage fish runs
  load(file=paste(fishuse,1,"nofish",".RData",sep="_"))
  load(file=paste(fishuse,2,"constF",".RData",sep="_"))
  load(file=paste(fishuse,3,"constFhi",".RData",sep="_"))
  load(file=paste(fishuse,4,"msycatch",".RData",sep="_"))
  load(file=paste(fishuse,5,"C1",".RData",sep="_"))
  load(file=paste(fishuse,6,"C2",".RData",sep="_"))
  load(file=paste(fishuse,7,"C3",".RData",sep="_"))
  nsims = nrow(nofish$fishing)
  load(file=paste(fishuse,"equilib",".RData",sep="_"))
  
  # run/source seabird runs - both with non-fish prey and with fishing scenarios
  source(file.path("Run_seabird_2lifehistory.R"))
  source(file.path("seabird-foragefish_fishingscenarios.R"))
  
  # median and confidence intervals of runs
  middle = c(median(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
              median(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
              median(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
              median(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
              median(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
              median(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
              median(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
              median(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])),
              median(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000])),
              median(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000])))
  lower = c(quantile(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.025),
             quantile(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.025))
  upper = c(quantile(rowMeans(seabirdoutC$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutC2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutE$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutE2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutF$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutF2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutB$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutB2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutD$totalpop[,201:1000]/seabirdout$totalpop[,201:1000]), 0.975),
             quantile(rowMeans(seabirdoutD2$totalpop[,201:1000]/seabirdout2$totalpop[,201:1000]), 0.975))
  
  assign(paste(tolower(fishuse),"middle",sep=""),  middle)
  assign(paste(tolower(fishuse),"lower",sep=""),  lower)
  assign(paste(tolower(fishuse),"upper",sep=""),  upper)
  
  # to test MSY fishing on flexible seabird
  assign(paste(tolower(fishuse),1,sep=""),  seabirdout2)
  assign(paste(tolower(fishuse),"MSY",sep=""),  seabirdoutmsy)
  
}

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
plot(anchovymiddle, x, frame = F, xlab = "", ylab = "", 
     xlim = c(0,1), axes = F, pch = 19, col = c("Black", "Red"), main = "Anchovy Fishery")
arrows(anchovylower, x, 
       anchovyupper, x, length=0, angle=90, code=3, lwd = 1,col = c("Black", "Red"))
#text(middleA, x, labels =round(middleA,2) )

axis(2, at = c(1,2,3,4,5), labels = names, las = 2, pos = -0.1, cex.axis = 1)
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), cex.axis=1, tck=0.01, pos = 0.8)
legend(-0.5,0.4,xpd=NA, col = c("Black", "Red"), pch = 19, legend = c("Restricted", "Flexible"))
#text(-0.4,5.4, label = "(A)", xpd = NA)

#par(mar = c(6,8.5,2,2))
#x = c(1,1,2,2,3,3,4,4,5,5)
#names = c("Hockey-stick -\n Fmax 0.5 Fmsy","Hockey-stick - \nLow Blim", "Hockey-stick - \nFmax Fmsy", "Constant Low - \n0.25 Fmsy", "Constant Mod. - \n0.5 Fmsy")
plot(sardinemiddle, x, frame = F, xlab = "", ylab = "", 
     xlim = c(0,1), axes = F, pch = 19, col = c("Black", "Red"), main = "Sardine Fishery")
arrows(sardinelower, x, 
       sardineupper, x, length=0, angle=90, code=3, lwd = 1,col = c("Black", "Red"))
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


# MSY on flexible seabird
####### what does convential fishing levels do to the generalist?

sar1 =  rowMeans(sardine1$totalpop[,201:1000])
sarmsy = rowMeans(sardineMSY$totalpop[,201:1000]/sardine1$totalpop[,201:1000])


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
