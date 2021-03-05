# Laura modified from Megsie's initial code on March 19th, 2020

basedir <- "/Users/laurakoehn/Dropbox/ff-mse2"
figwd <- file.path(basedir,"Plots/")

library(RColorBrewer)
library(scales)

# Demo figure for control rules
source(file.path(basedir,"Control Rules/smith_oceana.R"))
source(file.path(basedir,"Control Rules/cfp.R"))
source(file.path(basedir,"Control Rules/hockey-stick.R"))
source(file.path(basedir,"Control Rules/trend-based-rule.R"))



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

#setwd(figwd)
#pdf("ControlRules_v11.pdf",width = 7,height = 7,useDingbats = FALSE)
lwdp = 3
par(las=2,mar=c(6,6,2,2)+0.1) # Rotate axis labels 
# margins: c(5, 4, 4, 2) + 0.1
adj = 0.05/4
plot(C1.vec+adj,type='l',col=hcr.colors[1],lwd=lwdp, ylim=c(0,0.9),
     axes=FALSE,xlab="Biomass",ylab="Fishing rate (F) \n\n",xaxs="i",yaxs="i")
xticks = c(0,0.1*B0,0.5*Bmsy,0.4*B0,0.8*B0,max(B))
yticks = c(0,0.5*m,0.25*Fmsy, 0.5*Fmsy,Fmsy,0.8)
axis(side = 1, at = xticks,labels = c(" ","0.1B0","0.5Bmsy","0.4B0","0.8B0"," "))
axis(side = 2, at = yticks, labels = c("0","0.5M","0.25Fmsy", "0.5Fmsy","Fmsy"," "))

lines(C2.vec+adj+0.006, col=hcr.colors[2],lwd=lwdp)         # green
lines(C3.vec+adj+0.008, col=hcr.colors[3],lwd=lwdp)         # pale green
#lines(cfp.vec+adj+0.01,col=hcr.colors[6],lwd=lwdp)          # orange
lines(constf.vec,col=add.alpha(hcr.colors[4],alpha = 0.6),lwd=lwdp) # red line, constF.high
lines(constf.vec_lo,col=add.alpha(hcr.colors[5],alpha = 0.6),lwd=lwdp)  
legend("topleft",legend = c("Basic hockey stick","Low Blim","High Fmax","Moderate F","Low F"),bty = "n",lwd=rep(lwdp,times=4),col=hcr.colors[c(1,2,3,4,5)],lty=rep(1, times=6))
dev.off()



# functions used ----------------------------------------------------------

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# just hockey stick -------------------------------

lwdp = 3
par(las=2,mar=c(6,6,2,2)+0.1) # Rotate axis labels 
# margins: c(5, 4, 4, 2) + 0.1
adj = 0.05/4
plot(C1.vec+adj,type='l',col=hcr.colors[1],lwd=lwdp, ylim=c(0,0.9),
     axes = FALSE, xlab="Biomass",ylab="Fishing rate (F) \n\n",xaxs="i",yaxs="i")
xticks = c(min(B),max(B))
yticks = c(0,Fmsy,0.8)
axis(side = 1, at = xticks,labels = c("Low","High"), las = 1)
axis(side = 2, at = yticks, labels = c("0","Fmsy"," "))

lines(C2.vec+adj+0.006, col=hcr.colors[2],lwd=lwdp)         # green
lines(C3.vec+adj+0.008, col=hcr.colors[3],lwd=lwdp)         # pale green
#lines(cfp.vec+adj+0.01,col=hcr.colors[6],lwd=lwdp)          # orange
#lines(constf.vec,col=add.alpha(hcr.colors[4],alpha = 0.6),lwd=lwdp) # red line, constF.high
#lines(constf.vec_lo,col=add.alpha(hcr.colors[5],alpha = 0.6),lwd=lwdp)  
legend("topleft",legend = c("Moderate hockey stick","Low Biomass cut-off","High Max. Fishing"),
       bty = "n",lwd=rep(lwdp,times=4),col=hcr.colors[c(1,2,3,4,5)],lty=rep(1, times=6))