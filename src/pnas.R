## Human demographic matrix projections
## Corey J. A. Bradshaw & Brook W. Brook
## corey.bradshaw@adelaide.edu.au
## The University of Adelaide, Adelaide, Australia
## Sep 2014
## Remove everything
rm(list = ls())
## libraries
library(boot)
## matrix operators source fille (supplied separately)
source("matrixOperators.R")
## set working directory
#setwd("~/xxx/xxx/") # user to update
###############################################
## WHO-CHOICE 0-100 per 1-yr age structure
## import 2013 data
## US Census Bureau 2013 world life table
oldwd <- setwd("../data")
dat.world13 <- read.table("world2013lifetable.csv", header=T, sep=",") # available at http://dx.doi.org/10.4227/05/53869A9434A46
# import data
dat.afrD <- read.table("afrD.csv", header=T, sep=",")
dat.afrE <- read.table("afrE.csv", header=T, sep=",")
dat.AmA <- read.table("AmA.csv", header=T, sep=",")
dat.AmB <- read.table("AmB.csv", header=T, sep=",")
dat.AmD <- read.table("AmD.csv", header=T, sep=",")
dat.eMedB <- read.table("eMedB.csv", header=T, sep=",")
dat.eMedD <- read.table("eMedD.csv", header=T, sep=",")
dat.EurA <- read.table("EurA.csv", header=T, sep=",")
dat.EurB <- read.table("EurB.csv", header=T, sep=",")
dat.EurC <- read.table("EurC.csv", header=T, sep=",")
dat.SEAB <- read.table("SEAB.csv", header=T, sep=",")
dat.SEAD <- read.table("SEAD.csv", header=T, sep=",")
dat.wPacA <- read.table("wPacA.csv", header=T, sep=",")
dat.wPacB <- read.table("wPacB.csv", header=T, sep=",")

setwd(oldwd)

## combined age structure
mal05.N <- dat.afrD$mal.N + dat.afrE$mal.N + dat.AmA$mal.N + dat.AmB$mal.N + dat.AmD$mal.N + dat.eMedB$mal.N + dat.eMedD$mal.N + dat.EurA$mal.N +
dat.EurB$mal.N + dat.EurC$mal.N + dat.SEAB$mal.N + dat.SEAD$mal.N + dat.wPacA$mal.N + dat.wPacB$mal.N
fem05.N <- dat.afrD$fem.N + dat.afrE$fem.N + dat.AmA$fem.N + dat.AmB$fem.N + dat.AmD$fem.N + dat.eMedB$fem.N + dat.eMedD$fem.N + dat.EurA$fem.N +
dat.EurB$fem.N + dat.EurC$fem.N + dat.SEAB$fem.N + dat.SEAD$fem.N + dat.wPacA$fem.N + dat.wPacB$fem.N
## 2005 pop estimate
tot.pop05 <- sum(as.numeric(mal05.N+fem05.N))
## scale to 2013 pop estimate
tot.pop13 <- sum(as.numeric(dat.world13$N))
mal13.N <- round((tot.pop13/tot.pop05) * mal05.N, 0)
fem13.N <- round((tot.pop13/tot.pop05) * fem05.N, 0)
## N-weighted mortality rates
a.ages <- length(mal13.N)

mal.mort.mat <- as.matrix(data.frame(dat.afrD$mal.M,dat.afrE$mal.M,dat.AmA$mal.M,dat.AmB$mal.M,dat.AmD$mal.M,dat.eMedB$mal.M,dat.eMedD$mal.M,dat.EurA$mal.M,dat.EurB$mal.M,dat.EurC$mal.M,dat.SEAB$mal.M,dat.SEAD$mal.M,dat.wPacA$mal.M,dat.wPacB$mal.M))
lab.vec <- c("afrD","afrE","AmA","AmB","AmD","eMedB","eMedD","EurA","EurB","EurC","SEAB","SEAD","wPacA","wPacB")
colnames(mal.mort.mat) <- lab.vec
rownames(mal.mort.mat) <- 0:100
fem.mort.mat <-
as.matrix(data.frame(dat.afrD$fem.M,dat.afrE$fem.M,dat.AmA$fem.M,dat.AmB$fem.M,dat.AmD$fem.M,dat.eMedB$fem.M,dat.eMedD$fem.M,dat.EurA$fem.M,dat.EurB$fem.M,dat.EurC$fem.M,dat.SEAB$fem.M,dat.SEAD$fem.M,dat.wPacA$fem.M,dat.wPacB$fem.M))
lab.vec <- c("afrD","afrE","AmA","AmB","AmD","eMedB","eMedD","EurA","EurB","EurC","SEAB","SEAD","wPacA","wPacB")
colnames(fem.mort.mat) <- lab.vec
rownames(fem.mort.mat) <- 0:100
mal.N.mat <-
as.matrix(data.frame(dat.afrD$mal.N,dat.afrE$mal.N,dat.AmA$mal.N,dat.AmB$mal.N,dat.AmD$mal.N,dat.eMedB$mal.N,dat.eMedD$mal.N,dat.EurA$mal.N,dat.E urB$mal.N,dat.EurC$mal.N,dat.SEAB$mal.N,dat.SEAD$mal.N,dat.wPacA$mal.N,dat.wPacB$mal.N))
colnames(mal.N.mat) <- lab.vec
rownames(mal.N.mat) <- rownames(mal.mort.mat)
fem.N.mat <-
as.matrix(data.frame(dat.afrD$fem.N,dat.afrE$fem.N,dat.AmA$fem.N,dat.AmB$fem.N,dat.AmD$fem.N,dat.eMedB$fem.N,dat.eMedD$fem.N,dat.EurA$fem.N,dat.E urB$fem.N,dat.EurC$fem.N,dat.SEAB$fem.N,dat.SEAD$fem.N,dat.wPacA$fem.N,dat.wPacB$fem.N))
colnames(fem.N.mat) <- lab.vec
rownames(fem.N.mat) <- rownames(fem.mort.mat)
mal.wM <- fem.wM <- rep(0,a.ages)
for (a in 1:a.ages) {
    mal.wM[a] <- weighted.mean(mal.mort.mat[a,],mal.N.mat[a,])
    fem.wM[a] <- weighted.mean(fem.mort.mat[a,],fem.N.mat[a,])
}
lng.age.vec <- seq(0,100,1)
## construct matrix
stagesWHO <- a.ages
popmatWHO <- matrix(0,nrow=stagesWHO,ncol=stagesWHO)
colnames(popmatWHO) <- lng.age.vec[1:stagesWHO]
rownames(popmatWHO) <- lng.age.vec[1:stagesWHO]
## populate matrix
popmatWHO[1,1:length(dat.world13$m.f)] <- dat.world13$m.f
surv.vec2 <- 1-fem.wM
diag(popmatWHO[2:stagesWHO,]) <- surv.vec2[-stagesWHO]
#popmatWHO[2,1] <- surv.vec[1] ## infant mortality from dat.world13
popmatWHO[stagesWHO,stagesWHO] <- surv.vec2[stagesWHO]
popmatWHO.orig <- popmatWHO ## save original matrix
## matrix properties
max.lambda(popmatWHO) ## 1-yr lambda
max.r(popmatWHO) # rate of population change, 1-yr
stable.stage.dist(popmatWHO) ## stable stage distribution
R.val(popmatWHO,stagesWHO) # reproductive value
G.val(popmatWHO,stagesWHO) # mean generation length
## initial population vector
initWHO.vec <- fem13.N

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2013 # update if more data available post-2010
#************************
yr.end <- 2100 # set projection end date
#************************
t <- (yr.end - yr.now)
## linear fertility trends
int.vec <- 1:t
tot.F <- sum(popmatWHO.orig[1,])
popmatWHO <- popmatWHO.orig
#****************************************************
## Choose scenarios
## Fertility change
F.scen0 <- tot.F ## no fertility change
F.scen1 <- 1 # worldwide 1-child policy
F.scen2 <- 2.00 #
#************************
end.fert <- 2100
F.scen.choose <- F.scen0
#************************
## fertilty-change vector
F.scen.ch <- ifelse(F.scen.choose == tot.F, F.scen.choose, F.scen.choose/2)
F.mult <- rep(F.scen.ch/tot.F, t)
ft <- (end.fert - yr.now)
F.mult.vec <- pmax(F.mult, (tot.F - (tot.F - F.scen.ch)*int.vec/ft)/tot.F)
#****************************************************************************
## Unwanted pregancy aversion
preg.aversion <- 0 # 1 if avoid number of unwanted pregancies; 0 if not avoid
# 2008: 208 million pregnancies; 86 million unintended
# 33 million unplanned births; 41 million abortions; 11 million miscarriages
prop.unwanted <- 33/208 # worldwide proportion of potentially avoidable births
#****************************************************************************
#****************************************************************************
## Stepped fertility reduction (as per reviewer 1's recommended scenario)
## 1 child/female by 2045
stepped <- 0 # 1 if this stepped function invoked
fts <- (2045 - yr.now); ints.vec <- 1:fts; Fs.targ <- 1.0*0.5; mults <- Fs.targ/tot.F; F.mults <- rep(mults, fts)
step.fert.mult <- pmax(F.mults, (tot.F - (tot.F - Fs.targ)*ints.vec/fts)/tot.F)
step.fert.mult.vec <- c(step.fert.mult,rep(step.fert.mult[length(step.fert.mult)], ft-fts))
#****************************************************************************

#****************************************************************************
## Stepped fertility reduction 2nd scenario
## 2.0 child/female by 2020
stepped2 <- 0 # 1 if this stepped function invoked
fts2 <- (2020 - yr.now); ints2.vec <- 1:fts2; Fs.targ2 <- 2.0*0.5; mults2 <- Fs.targ2/tot.F; F.mults2 <- rep(mults2, fts2)
step.fert.mult2 <- pmax(F.mults2, (tot.F - (tot.F - Fs.targ2)*ints2.vec/fts2)/tot.F)
step.fert.mult2.vec <- c(step.fert.mult2,rep(step.fert.mult2[length(step.fert.mult2)], ft-fts2))
##****************************************************************************

## Age at primiparity changes (alpha)
A.scen0 <- 1
A.scen1 <- 0.5 # amount of fertility redistributed from 14:24 to 25:49
#************************
end.alpha <- 2100
A.scen.choose <- A.scen0
#************************
## Non-juvenile (6:oldest) survival change (S)
D.scen0 <- 1 ## no survival change
D.scen1 <- 0.50 # 50 % reduction in stage-specific death rate
#************************
end.death <- 2100
D.scen.choose <- D.scen0
#************************
## need to create using UN life-expectancy projections (back-calculated to death rates)
## Juvenile survival change (0-5 yrs)
J.scen0 <- 1 ## no survival change
J.scen1 <- 0.50 # 50 % reduction in juvenile death rate
J.scen2 <- 1.50 # 50 % increase in juvenile death rates (e.g., famine from CC)
#************************
end.juvd <- 2100
J.scen.choose <- J.scen0
#************************
## add a catastrophic mortality event
# spread over 5 years
# implemented mid-projection
# equal likelihood of taking any age (/2 for females only)
no.toll <- 0
firstsecwars.toll <- 1.31e+8/2 #131,000,000
firstsecwars.prop.toll <- (firstsecwars.toll*2/2500000000*tot.pop13)/2 # deaths proportional to 2.5 b alive at end WWII
twobillion.toll <- 2e+9/2 #2,000,000,000 dead
sixbillion.toll <- 6e+9/2 #6,000,000,000 dead
#************************
Cat.scen <- no.toll
#************************
## change to fertility/survival if war/pandemic invoked
war.fert.mult <- 2 # fertility doubles following war/pandemic (subsequently increases linearly until 2013 values thereafter)
war.surv.mult <- 2 # mortality doubles following war (subsequently increases linearly until 2013 values thereafter)
yr.vec <- seq(yr.now,yr.end)
if (Cat.scen != sixbillion.toll) {
    eyr.cat <- yr.vec[round(t/2)+5] # year + 1 end of catastrophe
}
if (Cat.scen == sixbillion.toll) {
    eyr.cat <- yr.vec[round(t/3)+5] # year + 1 end of catastrophe
}
ws.dt <- (end.death - eyr.cat)
mid.int.vec <- 1:ws.dt
war.surv.mult.vc <- rev(pmin(war.surv.mult, 1 - (1 - war.surv.mult)*mid.int.vec/ws.dt))
## add 1s to first part of vector for no change prior to war/disease catastrophe

add.surv.mult <- rep(1,t-ws.dt)
war.surv.mult.vec <- c(add.surv.mult,war.surv.mult.vc)
war.surv.mult.vec
war.F.mult.vc <- rev(pmin(war.fert.mult, 1 - (1 - war.fert.mult)*mid.int.vec/ws.dt))
add.fert.mult <- rep(1,t-ws.dt)
war.fert.mult.vec <- c(add.fert.mult,war.F.mult.vc)
war.fert.mult.vec
## changes in female age at primparity (alpha)
## done via redistribution of fertility in 15:24 breeding classes amongst remaining (25:49)
A.mult <- rep(0, t)
at <- (end.alpha - yr.now)
A.mult.vec <- A.scen.choose^(1/t)
A.mult.vec
## changes in survival
## expressed in death rate (1 - S)
D.mult <- rep(1, t)
dt <- (end.death - yr.now)
D.mult.vec <- pmax(D.scen.choose, 1 - (1 - D.scen.choose)*int.vec/dt)
D.mult.vec
## changes in juvenile survival (expressed in mortality changes)
J.mult <- rep(1, t)
jt <- (end.juvd - yr.now)
if (J.scen.choose <= 1) {
    J.mult.vec <- pmax(J.scen.choose, 1 - (1 - J.scen.choose)*int.vec/jt)}
if (J.scen.choose > 1) {
    J.mult.vec <- 1 - (1 - J.scen.choose)*int.vec/jt
}
J.mult.vec
## set population storage matrix
n.matWHO <- matrix(0,nrow=stagesWHO,ncol=(t+1))
n.matWHO[,1] <- initWHO.vec
## fertility storage vector
fertWHO.tot <- fertWHO.15.24 <- rep(0,t)
## set up projection loop
for (i in 1:t) {
    fertWHO.tot[i] <- sum(popmatWHO[1,])
    fertWHO.15.24[i] <- sum(popmatWHO[1,16:25])
    n.matWHO[,i+1] <- popmatWHO %*% n.matWHO[,i]

    ## invoke catastrophic mortality over 5-year window
    if (Cat.scen != sixbillion.toll) {
        if (i == round(t/2,0)) {
            if (Cat.scen == firstsecwars.prop.toll) {
                Cat.scen <- firstsecwars.toll*2/2500000000*sum(n.matWHO[,i])
            }
            prop.death <- Cat.scen/5/sum(n.matWHO[,i+1])
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death) ## five-year (Cat.scen/5) catastrophic mortality event
        }
        if (i == round(t/2,0)+1) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }

        if (i == round(t/2,0)+2) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
        if (i == round(t/2,0)+3) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
        if (i == round(t/2,0)+4) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
    }
    if (Cat.scen == sixbillion.toll) {
        if (i == round(t/3,0)) {
            prop.death <- Cat.scen/5/sum(n.matWHO[,i+1])
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death) ## five-year (Cat.scen/5) catastrophic mortality event
        }
        if (i == round(t/3,0)+1) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
        if (i == round(t/3,0)+2) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
        if (i == round(t/3,0)+3) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
        if (i == round(t/3,0)+4) {
            n.matWHO[,i+1] <- n.matWHO[,i+1] - (n.matWHO[,i+1]*prop.death)
        }
    }

    ## unwanted pregnancies averted
    if (preg.aversion == 1) {
        n.matWHO[1,i+1] <- n.matWHO[1,i+1] * (1-prop.unwanted) # n * 1-proportion unwanted
    }

    ## stepped fertility decline
    if (stepped == 1) {
        F.mult.vec <- step.fert.mult.vec
    }

    ## stepped fertility decline (2nd scenario)
    if (stepped2 == 1) {
        F.mult.vec <- step.fert.mult2.vec
    }

    ## change fertility/survival following catastrophic mortality event
    if (Cat.scen > 0) {
        popmatWHO[1,] <- popmatWHO.orig[1,]*war.fert.mult.vec[i]
        diag(popmatWHO[7:stagesWHO,6:stagesWHO]) <- 1 - ((1 - diag(popmatWHO.orig[7:stagesWHO,6:stagesWHO]))*war.surv.mult.vec[i])
        popmatWHO[stagesWHO,stagesWHO] <- 1 - ((1 - popmatWHO.orig[stagesWHO,stagesWHO])*war.surv.mult.vec[i])
        diag(popmatWHO[2:6,1:5]) <- 1 - ((1 - diag(popmatWHO.orig[2:6,1:5]))*war.surv.mult.vec[i])
    }

    popmatWHO[1,] <- popmatWHO.orig[1,]*F.mult.vec[i]
    diag(popmatWHO[7:stagesWHO,6:stagesWHO]) <- 1 - ((1 - diag(popmatWHO.orig[7:stagesWHO,6:stagesWHO]))*D.mult.vec[i])
    popmatWHO[stagesWHO,stagesWHO] <- 1 - ((1 - popmatWHO.orig[stagesWHO,stagesWHO])*D.mult.vec[i])
    diag(popmatWHO[2:6,1:5]) <- 1 - ((1 - diag(popmatWHO.orig[2:6,1:5]))*J.mult.vec[i])
    A.allocWHO <- sum(popmatWHO[1,16:25])*(1-A.mult.vec)
    F1524WHO.wt <- popmatWHO[1,16:25]/sum(popmatWHO[1,16:25])
    popmatWHO[1,16:25] <- popmatWHO[1,16:25] - F1524WHO.wt*A.allocWHO
    FRWHO.wt <- popmatWHO[1,26:50]/sum(popmatWHO[1,26:50])
    popmatWHO[1,26:50] <- popmatWHO[1,26:50] + FRWHO.wt*A.allocWHO
    ## revisit the mult.vecs - might be able to do faster as for A.alloc
}

## final population size
totWHO.sex.ratio <- 1/(sum(mal13.N)/(sum(fem13.N) + sum(mal13.N)))
fin.popWHO <- totWHO.sex.ratio*(sum(n.matWHO[,(t+1)]))
fin.popWHO
## end max lambda
max.lambda(popmatWHO) # 1-yr lambda
# max population size reached over projection interval
n.maxWHO <- totWHO.sex.ratio*max(colSums(n.matWHO))
n.maxWHO

## x change (initial - final)
times.deltaWHO <- round(fin.popWHO/(totWHO.sex.ratio*sum(initWHO.vec)), 2)
times.deltaWHO
## year projection vector
yrs <- seq(yr.now,yr.end,1)
## > 65 (or 75, for sensitivity test) proportion (of total population) (choose 65 or 75 upper threshold - comment/uncomment accordingly)
#over65WHO <- colSums(n.matWHO[67:stagesWHO,])/colSums(n.matWHO)
over65WHO <- colSums(n.matWHO[77:stagesWHO,])/colSums(n.matWHO)
## < 15 proportion (of total pop)
under15WHO <- colSums(n.matWHO[1:15,])/colSums(n.matWHO)
## dependency ratio
# number of people < 15 and > 65 relative to rest (choose 65 or 75 upper threshold - comment/uncomment accordingly)
#dep.ratioWHO <- (colSums(n.matWHO[1:15,]) + colSums(n.matWHO[67:stagesWHO,])) / colSums(n.matWHO[16:66,])
dep.ratioWHO <- (colSums(n.matWHO[1:15,]) + colSums(n.matWHO[77:stagesWHO,])) / colSums(n.matWHO[16:76,])

## plots
par(mfrow=c(2,2),yaxt="s")

plot(yrs, as.vector(colSums(n.matWHO)*totWHO.sex.ratio), type="l",
     xlab="year", ylab="N", ylim=c(0,1.02*n.maxWHO), xlim=c(yr.now,yr.end))
abline(h=sum(n.matWHO[,1]*totWHO.sex.ratio),lty=2)
maxN.sub <- which(colSums(n.matWHO)==max(colSums(n.matWHO)))
abline(v=yrs[maxN.sub],lty=2)
title(main=paste("final N = ", round(fin.popWHO/1e9, 4), " b",
                 "; max N = ", round(n.maxWHO/1e9, 4), " b", sep=""),
      sub = paste("delta = ", times.deltaWHO, "x", sep=""))

plot(yrs, as.vector(over65WHO), type="l",
     xlab="year", ylab="proportion", ylim=c(0,1))
lines(yrs,as.vector(under15WHO),lty=2)
title(main=paste("prop <15 (dashed) & >65 (solid) yrs", sep=""))

plot(yrs,as.vector(dep.ratioWHO),type="l",xlab="year",ylab="dependency ratio",ylim=c(0,1))
title(main=paste("dependency ratio: initial = ", round(dep.ratioWHO[1],4), "; final = ", round(dep.ratioWHO[t+1],4), sep=""))

barplot(n.matWHO[-stagesWHO,t+1],names.arg=rownames(popmatWHO[-stagesWHO]),beside=T,xlab="age class",ylab="final N",axes=T,horiz=F,axis.lty=1)

par(mfrow=c(1,1))

