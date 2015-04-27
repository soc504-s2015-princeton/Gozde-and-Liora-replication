##Replication Code for "The Unholy Trinity: Immigration, Unemployment and Extremism in Europe, 1980-2002"
##Jeffrey Javed and Noam Gidron
##May 8, 2011

##Load libraries.

library(Zelig)
library(lme4)
library(MASS)
library(nlme)
library(plm)
library(xtable)

##Load Arzheimer (2009) data.

data <- as.data.frame(read.delim("nonimp.tab.tsv"))

##Center variables (done as according to the author).

salienzmean.c <- data$salienzmean - 3.84568
rvar.c <- data$rvar - 21.75423

##Create country dummies.

country <- factor(data$sortcountry, labels=c("AT", "BE","DE-E","DE-W","DK","ES","FI","FR","GR","IT","LU","NL","NO","PT","SE"))

##Run with PQL.
 
out.pql <- glmmPQL(rexvote ~ male + age1 + age2 + age4 + mye1 + mye2 + farmerown + worker + retired + unemployed + zlrs + euschlecht + zsatisdmo + disp + lfed1 + zasylumseekers + zsur + zasylumseekers:zsur + zreplacementrate + zsur:zreplacementrate + zasylumseekers:zreplacementrate + rmax + salienzmean.c + rvar.c + rvar.c:salienzmean.c + country - 1, random=~1|kontext, family=binomial(link="logit"), data=data, verbose=TRUE)

summary(out.pql)

####Simulations for the Original Model with Centered Variables##

##Simulate betas for centered data set.

set.seed(02138)
simbetas <- mvrnorm(n=10000, mu=out.pql$coef$fixed, Sigma=out.pql$var)

##Create matrix for all variables.

mat1 <- cbind(data$rexvote, data$male, data$age1, data$age2, data$age4, data$mye1, data$mye2, data$farmerown, 
data$worker, data$retired, data$unemployed, data$zlrs, data$euschlecht, data$zsatisdmo, data$disp, data$lfed1, 
data$zasylumseekers, data$zsur, data$zreplacementrate, data$rmax, salienzmean.c, rvar.c, data$at, data$be, data$deo, 
data$dew, data$dk, data$es, data$fi, data$fr, data$gr, data$it, data$nl, data$no, data$pt, data$se, data$zasylumseekers*data$zsur, 
data$zsur*data$zreplacementrate, data$zasylumseekers*data$zreplacementrate, salienzmean.c*rvar.c)

##Omit data and create x matrix.

mat2 <- na.omit(mat1)[,-1]

##Create means for all variables

means <- as.vector(apply(mat2, 2, mean))

##Simulate change in prob of ER vote due to change in immigration (asylum seekers) with all other variables held at mean.

asy <- seq(from = -.98, to=4.46, by=.01)
ests <- matrix(data=NA, ncol=length(asy) ,nrow=10000)

for(j in 1:length(asy)){
	data.asy <- means
	data.asy[16] <- asy[j] 
	data.asy[36] <- asy[j]*data.asy[17]  ##Interaction between asylumseekers and zsur (17).
	data.asy[38] <- asy[j]*data.asy[18]  ##Interaction between asylumseekers and replacementrate (18). 
	ests[,j] <- (1+exp(-data.asy%*%t(simbetas)))^-1
	}

pdf(file="asylumseekers.pdf")
plot(NA,NA, ylim=c(0,.04), xlim=c(-1,5), xlab="Asylum Seekers (centered)", ylab="Probability of Voting for the ER", main="The Effect of Asylum Seekers on Voting for the ER", cex.main=1)
lines(asy, apply(ests, 2, mean))
lines(asy, apply(ests, 2, quantile, 0.025),col="red",lty=3)
lines(asy, apply(ests, 2, quantile, 0.975),col="red",lty=3)
abline(h=0)
dev.off()


####################Net Migration Variable########################

##Load migration data.

mig.data <- read.csv("migration.csv")

##Label migration data.

colnames(mig.data) <- c("year","migration","at","be","dk","fi","fr","deo","dew","gr","it","nl","no","pt","es","se","lu","sortcountry")

##Merge migration data with Arzheimer dataset by year and country.

new.data <- merge(data, mig.data, by=c("year","at","be","dk","fi","fr","deo","dew","gr","it","nl","no","pt","es","se","lu","sortcountry"))

##Create country dummies for merged dataset.

country.mig <- factor(new.data$sortcountry, labels=c("AT", "BE","DE-E","DE-W","DK","ES","FI","FR","GR","IT","LU","NL","NO","PT","SE"))

##Run with PQL.

out.mig <- glmmPQL(rexvote ~ male + age1 + age2 + age4 + mye1 + mye2 + farmerown + worker + retired + unemployed + zlrs + euschlecht + zsatisdmo + disp + lfed1 + migration + sur + migration:sur + replacementrate + sur:replacementrate + migration:replacementrate + rmax + salienzmean + rvar + rvar:salienzmean + country.mig - 1, random=~1|kontext, family=binomial(link="logit"), data=new.data, verbose=TRUE)

summary(out.mig)


#############Simulations with Migration Variable################

##Simulate betas.

set.seed(02138)
simbetas.mig <- mvrnorm(n=10000, mu=out.mig$coef$fixed, Sigma=out.mig$var)

##Create matrix for all variables.

mat.mig <- cbind(new.data$rexvote, new.data$male, new.data$age1, new.data$age2, new.data$age4, 
new.data$mye1, new.data$mye2, new.data$farmerown, new.data$worker, new.data$retired, new.data$unemployed, 
new.data$zlrs, new.data$euschlecht, new.data$zsatisdmo, new.data$disp, new.data$lfed1, new.data$migration, 
new.data$sur, new.data$replacementrate, new.data$rmax, new.data$salienzmean, new.data$rvar, new.data$at, new.data$be, 
new.data$deo, new.data$dew, new.data$dk, new.data$es, new.data$fi, new.data$fr, new.data$gr, new.data$it, new.data$nl, 
new.data$no, new.data$pt, new.data$se, new.data$migration*new.data$sur, new.data$sur*new.data$replacementrate, new.data$migration*new.data$replacementrate, 
new.data$salienzmean*new.data$rvar)

##Label matrix.

colnames(mat.mig) <- c("rexvote","male","age1","age2","age4","mye1","mye2","farmerown","worker","retired","unemployed","zlrs","euschlecht","zsatisdmo","disp","lfed1","migration","sur","replacementrate","rmax","salienzmean","rvar","AT", "BE","DE-E","DE-W","DK","ES","FI","FR","GR","IT","NL","NO","PT","SE","migration*sur", "sur*replacementrate","migration*replacementrate", "rvar*salienzmean")

##Omit missing data and create x matrix.

mat.mig2 <- na.omit(mat.mig)
mat.mig3 <- mat.mig2[,-1]

##Create means for all variables.

means.mig <- as.vector(apply(mat.mig3, 2, mean))

##Simulate change in prob of ER vote due to change in net migration rate with all other variables held at mean.

mig <- seq(from = -4, to=16.3, by=.1)
ests.mig <- matrix(data=NA, ncol=length(mig) ,nrow=10000)

for(j in 1:length(mig)){
	data.mig <- means.mig
	data.mig[16] <- mig[j] 
	data.mig[36] <- mig[j]*data.mig[17]  ##Interaction between migration and sur (17).
	data.mig[38] <- mig[j]*data.mig[18]  ##Interaction between migration and replacementrate (18). 
	ests.mig[,j] <- (1+exp(-data.mig%*%t(simbetas.mig)))^-1
	}

pdf(file="migration.pdf")
plot(NA,NA, ylim= c(0,.06),xlim=c(-4,16.5), xlab="Net Migration Rate", ylab="Probability of Voting for the ER", main="The Effect of Net Migration on Voting for the ER", cex.main=1)
lines(mig, apply(ests.mig, 2, mean))
lines(mig, apply(ests.mig, 2, quantile, 0.025),col="red",lty=3)
lines(mig, apply(ests.mig, 2, quantile, 0.975),col="red",lty=3)
abline(h=0)
dev.off()

################Foreign Population Variable##############

##Load foreign population data.

fp.data <- read.csv("foreign population.csv")

##Merge foreign population data with Arzheimer dataset by year and country.

new.data4 <- merge(data, fp.data, by=c("year", "sortcountry"))

##Run with PQL.

out.fp <- glmmPQL(rexvote ~ male + age1 + age2 + age4 + mye1 + mye2 + farmerown + worker + retired + unemployed + zlrs + euschlecht + zsatisdmo + disp + lfed1 + foreignpop + sur + foreignpop:sur + replacementrate + sur:replacementrate + foreignpop:replacementrate + rmax + salienzmean + rvar + rvar:salienzmean + factor(sortcountry) - 1, random=~1|kontext, family=binomial(link="logit"), data=new.data4, verbose=TRUE)

summary(out.fp)

############Simulations with Foreign Population Variable#######

##Simulate betas.

set.seed(02138)
simbetas.fp <- mvrnorm(n=10000, mu=out.fp$coef$fixed, Sigma=out.fp$var)

##Create matrix for all variables.

mat.fp <- cbind(new.data4$rexvote, new.data4$male, new.data4$age1, new.data4$age2, new.data4$age4, new.data4$mye1, new.data4$mye2, new.data4$farmerown, new.data4$worker, new.data4$retired, new.data4$unemployed, new.data4$zlrs, new.data4$euschlecht, new.data4$zsatisdmo, new.data4$disp, new.data4$lfed1, new.data4$foreignpop, new.data4$sur, new.data4$replacementrate, new.data4$rmax, new.data4$salienzmean, new.data4$rvar, new.data4$at, new.data4$be, new.data4$deo, new.data4$dew, new.data4$dk, new.data4$es, new.data4$fi, new.data4$it, new.data4$nl, new.data4$no, new.data4$pt, new.data4$se, new.data4$foreignpop*new.data4$sur, new.data4$sur*new.data4$replacementrate, new.data4$foreignpop*new.data4$replacementrate, new.data4$salienzmean*new.data4$rvar)

##Label matrix.

colnames(mat.fp) <- c("rexvote","male","age1","age2","age4","mye1","mye2","farmerown","worker","retired","unemployed","zlrs","euschlecht","zsatisdmo","disp","lfed1","foreignpop","sur","replacementrate","rmax","salienzmean","rvar","AT", "BE","DE-E","DE-W","DK","ES","FI","IT","NL","NO","PT","SE","foreignpop*sur", "sur*replacementrate","foreignpop*replacementrate", "rvar*salienzmean")

##Omit missing data and create x matrix.

mat.fp2 <- na.omit(mat.fp)
mat.fp3 <- mat.fp2[,-1]

##Create means for all variables

means.fp <- as.vector(apply(mat.fp3, 2, mean))

##Simulate change in prob of ER vote due to change in foreign population with all other variables held at mean.

fp <- seq(from = 0.6, to=9.3, by=.1)
ests.fp <- matrix(data=NA, ncol=length(fp) ,nrow=10000)

for(j in 1:length(fp)){
	data.fp <- means.fp
	data.fp[16] <- fp[j] 
	data.fp[34] <- fp[j]* data.fp[17]  ##Interaction between fp and sur (17).
	data.fp[36] <- fp[j]* data.fp[18]  ##Interaction between fp and replacementrate (18). 
	ests.fp[,j] <- (1+exp(-data.fp%*%t(simbetas.fp)))^-1
	}

pdf(file="foreignpop.pdf")
plot(NA,NA, ylim= c(0,.2),xlim=c(0.5,9.3), xlab="Foreign Population (%)", ylab="Probability of Voting for the ER", main="The Effect of Foreign Population on Voting for the ER", cex.main=1)
lines(fp, apply(ests.fp, 2, mean))
lines(fp, apply(ests.fp, 2, quantile, 0.025),col="red",lty=3)
lines(fp, apply(ests.fp, 2, quantile, 0.975),col="red",lty=3)
abline(h=0)
dev.off()


#################Marginal Effects####################

##Interaction between asylum seekers (centered) and unemployment.

unemp <- seq(from = -5.9, to=12.29, by=.01)
ests.unemp <- matrix(data=NA, ncol=length(unemp) ,nrow=10000)

for(j in 1:length(unemp)){
	data.unemp <- means
	data.unemp[17] <- unemp[j] 
	ests.unemp[,j] <- data.unemp[17]*(simbetas[,36]) + simbetas[,16]
	}


asy <- seq(from = -0.98, to=4.46, by=.01)
ests.asy <- matrix(data=NA, ncol=length(asy) ,nrow=10000)

for(j in 1:length(asy)){
	data.asy <- means
	data.asy[16] <- asy[j]
	ests.asy[,j] <- data.asy[16]*(simbetas[,36]) + simbetas[,17]
	}

pdf(file="interaction_asy.pdf")
par(mfrow=c(2,1))
plot(NA, NA, xlab="Asylum Seekers", ylab="Effect: Unemployment", xlim=c(-1,4.5), ylim=c(-.2,.2), main="Marginal Effect of Unemployment")
curve(mean(simbetas[,17]) + mean(simbetas[,36])*x, add=TRUE, col="red")
lines(asy, apply(ests.asy, 2, quantile, 0.025),col="red", lty=3)
lines(asy, apply(ests.asy, 2, quantile, 0.975),col="red", lty=3)
abline(h=0, lwd=1)
rug(mat2[,16])
plot(NA, NA, xlab="Unemployment", ylab="Effect: Asylum Seekers", xlim=c(-6,12.5), ylim=c(-.5,.6), main="Marginal Effect of Asylum Seekers")
curve(mean(simbetas[,16]) + mean(simbetas[,36])*x, add=TRUE, col="blue")
lines(unemp, apply(ests.unemp, 2, quantile, 0.025),col="blue", lty=3)
lines(unemp, apply(ests.unemp, 2, quantile, 0.975),col="blue", lty=3)
abline(h=0, lwd=1)
rug(mat2[,17])
dev.off()


##Interaction effect between Net Migration and unemployment.

unemp.mig <- seq(from = 1.6, to=19.8, by=.1)
ests.unemp.mig <- matrix(data=NA, ncol=length(unemp.mig) ,nrow=10000)

for(j in 1:length(unemp.mig)){
	data.unemp.mig <- means.mig
	data.unemp.mig[17] <- unemp.mig[j]
	ests.unemp.mig[,j] <- data.unemp.mig[17]*(simbetas.mig[,36]) + simbetas.mig[,16] + simbetas.mig[,38]*data.unemp.mig[18]
	}

mig <- seq(from = -4, to=16.3, by=.1)
ests.mig <- matrix(data=NA, ncol=length(mig) ,nrow=10000)

for(j in 1:length(mig)){
	data.mig <- means.mig
	data.mig[16] <- mig[j]
	ests.mig[,j] <- data.mig[16]*(simbetas.mig[,36]) + simbetas.mig[,17] + simbetas.mig[,37]*data.mig[18]
	}

pdf(file="interaction_mig.pdf")
par(mfrow=c(2,1))
plot(NA, NA, xlab="Net Migration", ylab="Effect: Unemployment", xlim=c(-4, 16.3), ylim=c(-0.2,.4), main="Marginal Effect of Unemployment")
curve(mean(simbetas.mig[,17]) + mean(simbetas.mig[,36])*x + mean(simbetas.mig[,37])*data.mig[18], add=TRUE, col="red")
lines(mig, apply(ests.mig, 2, quantile, 0.025),col="red", lty=3)
lines(mig, apply(ests.mig, 2, quantile, 0.975),col="red", lty=3)
abline(h=0, lwd=1)
rug(mat.mig3[,16])
plot(NA, NA, xlab="Unemployment", ylab="Effect: Net Migration", xlim=c(1,20), ylim=c(-0.1,.3), main="Marginal Effect of Net Migration")
curve(mean(simbetas.mig[,16]) + mean(simbetas.mig[,36])*x + mean(simbetas.mig[,38])*data.unemp.mig[18], add=TRUE, col="blue")
lines(unemp.mig, apply(ests.unemp.mig, 2, quantile, 0.025),col="blue", lty=3)
lines(unemp.mig, apply(ests.unemp.mig, 2, quantile, 0.975),col="blue", lty=3)
abline(h=0, lwd=1)
rug(mat.mig3[,17])
dev.off()

##Interaction effect between Foreign Population and unemployment.

unemp.fp <- seq(from = 1.6, to=19.8, by=.1)
ests.unemp.fp <- matrix(data=NA, ncol=length(unemp.fp) ,nrow=10000)

for(j in 1:length(unemp.fp)){
	data.unemp.fp <- means.fp
	data.unemp.fp[17] <- unemp.fp[j]
	ests.unemp.fp[,j] <- data.unemp.fp[17]*(simbetas.fp[,34]) + simbetas.fp[,16] + simbetas.fp[,36]*data.unemp.fp[18]
	}
	
fp <- seq(from = 0.6, to=9.3, by=.1)
ests.fp <- matrix(data=NA, ncol=length(fp) ,nrow=10000)

for(j in 1:length(fp)){
	data.fp <- means.fp
	data.fp[16] <- fp[j]
	ests.fp[,j] <- data.fp[16]*(simbetas.fp[,34]) + simbetas.fp[,17] + simbetas.fp[,35]*data.fp[18]
	}
	
pdf(file="interaction_fp.pdf")
par(mfrow=c(2,1))
plot(NA, NA, xlab="Foreign Population (%)", ylab="Effect: Unemployment", xlim=c(0.6, 9.3), ylim=c(-0.3,.2), main="Marginal Effect of Unemployment")
curve(mean(simbetas.fp[,17]) + mean(simbetas.fp[,34])*x + mean(simbetas.fp[,35])*data.fp[18], add=TRUE, col="red")
lines(fp, apply(ests.fp, 2, quantile, 0.025),col="red", lty=3)
lines(fp, apply(ests.fp, 2, quantile, 0.975),col="red", lty=3)
abline(h=0, lwd=1)
rug(mat.fp3[,16])
plot(NA, NA, xlab="Unemployment", ylab="Effect: Foreign Population", xlim=c(1.6,20), ylim=c(-0.1,1.5), main="Marginal Effect of Foreign Population")
curve(mean(simbetas.fp[,16]) + mean(simbetas.fp[,34])*x + mean(simbetas.fp[,36])*data.unemp.fp[18], add=TRUE, col="blue")
lines(unemp.fp, apply(ests.unemp.fp, 2, quantile, 0.025),col="blue", lty=3)
lines(unemp.fp, apply(ests.unemp.fp, 2, quantile, 0.975),col="blue", lty=3)
abline(h=0, lwd=1)
rug(mat.fp3[,17])
dev.off()

###########Tables

##Create vectors of coefficients and standard deviations for each model.

asy.coef <- as.vector(out.pql$coef$fixed)
asy.sd <- as.vector(sqrt(diag(out.pql$varFix)))
mig.coef <- as.vector(out.mig$coef$fixed)
mig.sd <- as.vector(sqrt(diag(out.mig$varFix)))
fp.coef <- c(as.vector(out.fp$coef$fixed), 0, 0)
fp.sd <- c(as.vector(sqrt(diag(out.fp$varFix))), 0, 0)

##Create matrix of coefficients and standard deviations.

table.mat <- cbind(asy.coef, asy.sd, mig.coef, mig.sd, fp.coef, fp.sd)

##Label matrix.

rownames(table.mat) <- c("male","age1","age2","age4","mye1","mye2","farmerown","worker","retired","unemployed","zlrs","euschlecht","zsatisdmo","disp","lfed1","immigration","sur","replacementrate","rmax","salienzmean","rvar","AT", "BE","DE-E","DE-W","DK","ES","FI","FR","GR","IT","NL","NO","PT","SE","immigration*sur", "sur*replacementrate","immigration*replacementrate", "rvar*salienzmean")

##Create table.

xtable(table.mat, digits=3)
