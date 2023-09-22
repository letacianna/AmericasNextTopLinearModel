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
anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_smooth(method="lm") ## This can be taken out
anole.lm <- lm(HTotal~SVL,anole2)
coef(anole.lm)
anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_abline(slope=coef(anole.lm)[2],intercept=coef(anole.lm)[1],col="blue")
SVL2 <- seq(min(anole2$SVL),max(anole2$SVL),0.1)

pred.lm <-tibble(
  SVL=SVL2,
  H.pred=predict(anole.lm,newdata = data.frame(SVL=SVL2))
)

anole2%>%
  ggplot(aes(SVL,HTotal))+geom_point()+geom_point(data=pred.lm,aes(SVL,H.pred),col="blue") ## this can be taken out
summary(anole.lm)
anole.allo <- nls(HTotal~a*SVL^b, start=list(b=1, a=1),data = anole2)

summary(anole.allo)
#AICc from the MuMIn package
anole.aic <- AICc(anole.lm,anole.allo)

#aicw from the geiger package
anole.aicw <- aicw(anole.aic$AICc)

print(anole.aicw)
anole.log%>%
  ggplot(aes(HTotal,SVL,col=Ecomorph2))+geom_point()+geom_smooth(method="lm")
anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2,anole.log)
summary(anole.log.eco.lm)
anova(anole.log.eco.lm)
anole.log.lm  <- lm(HTotal~SVL,anole.log)
anova(anole.log.lm)
anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)
anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm))
anole.log%>%
  ggplot(aes(Ecomorph2,res))+geom_point()
p.eco <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=res)) +geom_boxplot()
print(p.eco)
p.eco+ geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)
anole.tree <- read.tree("anole.tre")
plot(anole.tree,cex=0.4)
#PGLS under BM, no ecomorph
pgls.BM1 <- gls(HTotal ~SVL, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under BM, w ecomorph
pgls.BM2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")


#PGLS under OU, no ecomorph
pgls.OU1 <- gls(HTotal ~SVL, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#PGLS under OU, w, ecomorph
pgls.OU2 <- gls(HTotal ~SVL * Ecomorph2, correlation = corMartins(0,phy = anole.tree,form=~Species),data = anole.log, method = "ML")
anova(pgls.BM2)
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM2))

p.eco.phylo <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=phylo.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

print(p.eco.phylo)
anole.log%>%
  dplyr::select(Ecomorph2,res,phylo.res)%>%
  pivot_longer(cols=c("res","phylo.res"))%>%
  print%>%
  ggplot(aes(x=Ecomorph2,y=value)) +geom_boxplot() +
  stat_summary(fun=mean, geom="point", size=3)+
  facet_grid(name~.,scales = "free_y")+ylab("residual")
anole.log
## Question 2 ## - Generate two simple linear models that assess the effect of perch diameter and height
PerchHeightModel <- lm(HTotal~SVL*PH,anole.log)
PerchDiameterModel <- lm(HTotal~SVL*ArbPD,anole.log)

##Question 3## - Plot Residuals For The Above Models
PHMRes <- PerchHeightModel$residuals
PDMRes <- PerchDiameterModel$residuals
anole.log %>% 
  mutate(PHMRes,PDMRes)
ggplot(anole.log, aes(PH,PHMRes)) + geom_point()
ggplot(anole.log,aes(ArbPD,PDMRes)) + geom_point()

#question 4 - phylogenetic least squares models of hindlimn-SVL relationships
#PGLS model with the hindlimb-SVL relationship + perch height
pgls.BM3 <- gls(HTotal ~SVL * PHMRes, correlation = corBrownian(1,phy=anole.tree,form=~Species),data=anole.log, method = "ML")
#PGLS model with the hindlimb-SVL relationship + perch diameter
pgls.BM4 <- gls(HTotal ~SVL *PDMRes, correlation = corBrownian(1,phy=anole.tree,form=~Species),data=anole.log, method = "ML")
#PGSL model with the hindlimb-SVL relationship + perch height + perch diameter
pgls.BM5 <- gls(HTotal ~SVL *PHMRes *PDMRes, correlation = corBrownian(1,phy=anole.tree,form=~Species),data=anole.log, method = "ML")

#question 5 - AIC, the most negative fit constitutes as a good fit. Therefore, pgls.BM5 is the best model for our data. Furthermore, the combined relationship of perch height and perch diameter are significant predictors of hindlimb length.
anole.phylo2.aic <- AICc(pgls.BM3,pgls.BM4,pgls.BM5)
aicw(anole.phylo2.aic$AICc)

#question 6
anova(pgls.BM5)
anole.log <- anole.log%>%
  mutate(phylo.res=residuals(pgls.BM5))

p.eco.phylo2 <- anole.log%>%
  ggplot(aes(x=Ecomorph2,y=phylo.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)

print(p.eco.phylo2)
