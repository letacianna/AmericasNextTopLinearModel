library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)
library(dplyr)

#load data
anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")
anole2 <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()%>%
  print()
anole.log <- anole2%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)
anole.lm <- lm(HTotal~SVL,anole2)

anole.allo <- nls(HTotal~a*SVL^b, start=list(b=1, a=1),data = anole2)

anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2,anole.log)
summary(anole.log.eco.lm)
anova(anole.log.eco.lm)
anole.log.lm  <- lm(HTotal~SVL,anole.log)
anova(anole.log.lm)
anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)
anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))
anole.tree <- read.tree("anole.tre")

#PGLS under BM, w ecomorph
pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))


anole.log
## Question 2 ## - Generate two simple linear models that assess the effect of perch diameter and height
PerchHeightModel <- lm(HTotal~SVL+PH,anole.log)
PerchDiameterModel <- lm(HTotal~SVL+ArbPD,anole.log)

##Question 3## - Plot Residuals For The Above Models
PHMRes <- PerchHeightModel$residuals
PDMRes <- PerchDiameterModel$residuals
anole.log <- anole.log%>% 
  mutate(PHMRes,PDMRes)
ggplot(anole.log, aes(PH,PHMRes)) + geom_point()
ggplot(anole.log,aes(ArbPD,PDMRes)) + geom_point()

#question 4 - phylogenetic least squares models of hindlimn-SVL relationships
#PGLS model with the hindlimb-SVL relationship + perch height
pgls.BM3 <- gls(HTotal ~SVL + PH, correlation = corBrownian(1,phy=anole.tree,form=~Species),data=anole.log, method = "ML")
#PGLS model with the hindlimb-SVL relationship + perch diameter
pgls.BM4 <- gls(HTotal ~SVL + ArbPD, correlation = corBrownian(1,phy=anole.tree,form=~Species),data=anole.log, method = "ML")
#PGSL model with the hindlimb-SVL relationship + perch height + perch diameter
pgls.BM5 <- gls(HTotal ~SVL + ArbPD + PH, correlation = corBrownian(1,phy=anole.tree,form=~Species),data=anole.log, method = "ML")

#question 5 - AIC, the most negative fit constitutes as a good fit. Therefore, pgls.BM5 is the best model for our data. Furthermore, the combined relationship of perch height and perch diameter are significant predictors of hindlimb length.
anole.phylo2.aic <- AICc(pgls.BM3,pgls.BM4,pgls.BM5)
aicw(anole.phylo2.aic$AICc)

#question 6
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM5))

HindlimbResPlot <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=phylo.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3) +xlab("Ecomorphs") + ylab("Hindlimb Model Residuals")

print(HindlimbResPlot)

#CPK: Just outstanding! Excellent work.
#CPK: Total points: 20