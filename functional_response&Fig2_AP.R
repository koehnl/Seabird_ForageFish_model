# CODE TO PLAY AROUND WITH FUNCTIONAL RESPONSES
# and make Figure 2 from Koehn et al. 2021
# acknowledgments to Andre Punt

alpha <- 0 # lowest possible impact on survival 
slope <-30 # how fast it drops from 1 to the lowest possible impact on suvival
beta <- 0.2 # ratio (P/P0) where predator would switch?

# P/P0
xx <- seq(from=0,to=1.5,by=0.01)

yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
plot(xx,yy,type='l')
abline(v = 0.15, lty = 3)

######################################
alpha <- 0
slope <- 40
beta <- 0.15

# P/P0
xx <- seq(from=0,to=1.5,by=0.01)

yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
plot(xx,yy,type='l')

alpha <- 0
slope <- 30
beta <- 0.2

# P/P0

yy2 <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx, yy2, col = "red")

functionalresponse <- function(alpha, slope, beta) {
  xx <- seq(from=0,to=1.1,by=0.01)
  yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
  plot(xx,yy,type='l', ylim = c(0,1), ylab = "", xlab = "", yaxt = "n", xaxt = "n")
  #abline(v = 0.3, lty = 3)
  axis(side = 1, at = c(0,0.1,0.2,0.3,0.4,0.5,1), las = 0)
  axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), las = 2)
  return(cbind(yy,xx))
}
xx <- seq(from=0,to=1.1,by=0.01)
functionalresponse(0.1, 40, 0.3) # test
zz <- seq(from=0,to=1.1,by=0.01)
lines(zz,xx)

par(mfrow = c(3,2))
par(mar=c(3,4,3,2))
par(oma=c(3,3,2,2))
# BREEDER ATTENDANCE
temp =functionalresponse(0.1,20,0.3)
functionalresponse(0.1,30,0.3) #breeders specialist
#functionalresponse(0.6,20,0.3) #generalist
#functionalresponse(0.3,20,0.3) #middle

alpha = 0.6; slope = 20; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Breeder Attendance", side=3, line = 1)
#mtext("Proportion Prey Available", side=1, line = 3)
alpha = 0.3; slope = 20; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("Specialist", "Generalist"), 
       col = c("Black", "Red"), lty = c(1,1), y.intersp=1, bty = 'n')

# fledgling
functionalresponse(-0.3,30,0.2) # have a max(0, yy) so will give 0 if yy goes negative
alpha = 0.51; slope = 15; beta = 0.1
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Fledge 1 chick - if 3 chicks", side=3, line = 1)
#mtext("Proportion Prey Available", side=1, line = 3)
alpha = 0.27; slope = 20; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("Specialist", "Generalist", "Middle"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1)

# Second chick
functionalresponse(-0.3,30,0.15)
# vs. firrst chick
alpha = -0.3; slope = 30; beta = 0.2
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")

functionalresponse(-0.3,30,0.15)
alpha = 0.41; slope = 15; beta = 0.05
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Fledge 1 chick - if 2 chicks", side=3, line = 1)
#mtext("Proportion Prey Available", side=1, line = 3)
alpha = 0.21; slope = 20; beta = 0.1
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("Specialist", "Generalist", "Middle"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1)

# 3rd chick
functionalresponse(-0.3,30,0.1)
alpha = 0.2; slope = 15; beta = 0
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Fledge 1 chick", side=3, line = 1)
#mtext("Proportion Prey Available", side=1, line = 3)
alpha = 0.05; slope = 20; beta = 0.05
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("Specialist", "Generalist", "Middle"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1)

param_breeders_special = c(0,20,0.3) # colony/breeder attendance (not survival)
param_breeders_general = c(0.5,20,)
param_breeders_mid = c(0.25,20,)
# OK? based on Piatt et al. 2007/Cairns should be shifted even further right. 
# also right shape?

param_fledge_special = c(-0.3,30,0.2)
# enough like 1/3 for the birds - i think so, and looks like Andre's
param_fledge_general = c(0.3,15,0.1) # so at some lower survival, switch and survival goes back to 1?
param_fledge_mid = c(0.1,20,0.15)

param_fledge2_special = c(-0.3,30,0.10)
paramfledge2_general = c(0.2,20,0.05) # *will 2nd and 3rd chick still apply for generalist who just switched anyway?
paramfledge2_mid = c()

paramfledge3_special = c(-0.3,30,0.05) 

# adult survival
functionalresponse(0.1,20,0.15)
alpha = 0.6; slope = 20; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Adult Survival", side=3, line = 1)
#mtext("Proportion Prey Available", side=1, line = 3)
alpha = 0.3; slope = 20; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("Specialist", "Generalist", "Middle"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1)

param_adult_special = c(-0.1,20,0.15)  # pretty good based on Piatt et al. 2007
param_adult_general = c(0.5,20,0.15)
param_adult_mid = c(0.25,20,0.15)

functionalresponse(0.1,10,0.3)
alpha = 0.6; slope = 10; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Juvenile Survival", side=3, line = 1)
#mtext("Proportion Prey Available", side=1, line = 3)
alpha = 0.3; slope = 10; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("Specialist", "Generalist", "Middle"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1)

par(xpd = NA)
mtext(expression(paste('Relative Prey Availability (P'['y,s,l'],' / P'[0],')')),
      side = 1, outer = TRUE, line = 1)
mtext(expression(paste('Demographic Rate Modifier (',delta['y,s,l'],')')), side = 2, outer = TRUE, line = 0, padj = 0)


param_juv_special = c(0,10,0.3) # need lots of prey or else only the strong survive
param_juv_general = c(0.5,10,0.3)
param_juv_mid = c(0.25,10,0.3)

# 1st vs. 2nd vs. 3rd chick
functionalresponse(-0.3,30,0.2) # have a max(0, yy) so will give 0 if yy goes negative
alpha = -0.3; slope = 30; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")
mtext("% Chick Fledge", side=2, line = 3)
mtext("Proportion Prey Available", side=1, line = 3)
alpha = -0.3; slope = 30; beta = 0.1
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")
legend("bottomright", legend = c("1st chick", "2nd chick", "3rd chick"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1)


functionalresponse(0.1,20,0.15)
alpha = 0.1; slope = 20; beta = 0.1
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "grey50")
alpha = 0.1; slope = 20; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "grey60")
alpha = 0.3; slope = 20; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "grey70")
alpha = 0.5; slope = 20; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "grey30")

######################### new plot 3/3/2020

functionalresponse <- function(alpha, slope, beta) {
  xx <- seq(from=0,to=0.5,by=0.01)
  yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
  plot(xx,yy,type='l', ylim = c(0,1), ylab = "", xlab = "", yaxt = "n", xaxt = "n")
  #abline(v = 0.3, lty = 3)
  axis(side = 1, at = c(0,0.1,0.2,0.3,0.4,0.5,1), las = 0)
  axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), las = 2)
  return(cbind(yy,xx))
}
xx <- seq(from=0,to=0.5,by=0.01)

par(xpd = FALSE)
par(mfrow = c(2,1))
par(mar = c(2,2,1,1))
par(oma = c(3,4,1,1))
functionalresponse(0.2,30,0.3) #breeders specialist
alpha = 0.6; slope = 30; beta = 0.3
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "black", lty = 2)
#mtext(expression(paste("Rate Modifier (", delta, ")" )), side=2, line = 2.5)
#mtext("Proportion Prey Available", side=1, line = 3)
#breeder = expression(paste("Breeder attendance", gamma))
legend("bottomright", 
       legend = c("Specialist", "Generalist", expression(paste("Breeder attendance ",gamma)), "Adult survival", "Juvenile survival"), 
       col = c("grey", "grey","black", "red", "blue"), lty = c(1,2,1,1,1), y.intersp=1, bty = 'n')

alpha = 0.6; slope = 30; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red", lty = 2)
alpha = 0.2; slope = 30; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red", lty = 1)

alpha = 0.6; slope = 30; beta = 0.2
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue", lty = 2)
alpha = 0.2; slope = 30; beta = 0.2
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue", lty = 1)

text(0,0.95, label = "(A)", xpd = NA)

functionalresponse(-0.3,30,0.2) # have a max(0, yy) so will give 0 if yy goes negative

alpha = 0.51; slope = 15; beta = 0.1
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "black", lty = 2)

alpha = -0.3; slope = 30; beta = 0.15
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red")

alpha = 0.41; slope = 15; beta = 0.05
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "red", lty = 2)

alpha = -0.3; slope = 30; beta = 0.1
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue")

alpha = 0.2; slope = 15; beta = 0
yy <- alpha + (1-alpha)/(1.0+exp(-slope*(xx-beta)))
lines(xx,yy, col = "blue", lty = 2)
legend("bottomright", legend = c("Fledge 1 of 3 chicks", "Fledge 1 of 2 chicks", "Fledge 1 of 1 chick"), 
       col = c("Black", "Red", "Blue"), lty = c(1,1,1), y.intersp=1, bty = 'n')
abline(v = 0.33, col = 'grey', lty = 4)
par(xpd = NA)
mtext(expression(paste('Relative Prey Availability (P'['y,s,l'],' / ',tilde(P)['l'],')')),
      side = 1, outer = TRUE, line = 1)
mtext(expression(paste('Demographic Rate Modifier (',delta['y,s,l'],' or ',gamma['y,s,l'], ')')), side = 2, outer = TRUE, line = 1, padj = 0)

text(0,0.95, label = "(B)", xpd = NA)


# param_fledge_special = c(-0.3,30,0.2)
# # enough like 1/3 for the birds - i think so, and looks like Andre's
# param_fledge_general = c(0.51,15,0.1) # so at some lower survival, switch and survival goes back to 1?
# 
# param_fledge2_special = c(-0.3,30,0.15)
# param_fledge2_general = c(0.41,15,0.05) # *will 2nd and 3rd chick still apply for generalist who just switched anyway?
# 
# param_fledge3_special = c(-0.3,30,0.1) 
# param_fledge3_general = c(0.2,15,0) 
