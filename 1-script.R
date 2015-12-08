setwd("~/Dropbox/projects/gss-guns")

library(foreign)
library(car)
library(arm)
library(plyr)
library(lme4)
library(blme)
# library(ordinal)
library(stargazer)
library(reshape2)
library(Zelig)
library(ZeligMultilevel)

require(ggplot2)
ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE, reorder=TRUE) {
  require(ggplot2)
  f <- function(x) {
    pv   <- attr(x, "postVar")
    cols <- 1:(dim(pv)[1])
    se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
    if (reorder) {
      ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
      pDf  <- data.frame(y=unlist(x)[ord],
                         ci=1.645*se[ord],
                         nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                         ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                         ind=gl(ncol(x), nrow(x), labels=names(x)))
    } else {
      pDf  <- data.frame(y=unlist(x),
                         ci=1.645*se,
                         nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                         ID=factor(rep(rownames(x), ncol(x)), levels=rownames(x)),
                         ind=gl(ncol(x), nrow(x), labels=names(x)))
    }

    if(QQ) {  ## normal QQ-plot
      p <- ggplot(pDf, aes(nQQ, y))
      p <- p + facet_wrap(~ ind, scales="free")
      p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
    } else {  ## caterpillar dotplot
      p <- ggplot(pDf, aes(ID, y)) + coord_flip()
      if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
        p <- p + facet_wrap(~ ind)
      } else {           ## different scales for random effects
        p <- p + facet_grid(ind ~ ., scales="free_y")
      }
      p <- p + xlab("Levels of the Random Effect") + ylab("Random Effect")
    }

    p <- p + theme(legend.position="none")
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
    p <- p + geom_point(aes(size=1.2), colour="blue") 
    return(p)
  }

  lapply(re, f)
}

GSS7214 <- read.dta("~/Dropbox/data/gss/GSS7214_R4.DTA", convert.factors = FALSE)

GSS <- with(GSS7214, data.frame(region, gun, gunlaw, gunimp, guninfo, gunfirm, owngun, rowngun, gun12, gunsales, gunsdrug, semiguns, guns911, othguns, shotgun, hgunlaw, hguncrim, crimup, crimdown, gunage, gunsdrnk, pistol, rifle, age, educ, degree, race, rincome, rincom06, partyid, sex, ethnic, year))

GSS$regioncondensed <- NA
GSS$regioncondensed <- with(GSS, ifelse(region == 8 | region == 9, "West", regioncondensed))
GSS$regioncondensed <- with(GSS, ifelse(region == 3 | region == 4, "Midwest", regioncondensed))
GSS$regioncondensed <- with(GSS, ifelse(region == 5 | region == 6 | region == 7, "South", regioncondensed))
GSS$regioncondensed <- with(GSS, ifelse(region == 1 | region == 2, "Northeast", regioncondensed))

GSS$gun <- with(GSS, recode(gun, "1=1; 2=0"))
GSS$gunlaw <- with(GSS, recode(gunlaw, "1=1; 2=0"))

GSS$gunmostimp <- with(GSS, recode(gunimp, "1=1; 2:4=0"))
GSS$gunimpd <- with(GSS, recode(gunimp, "1:2=1; 3:4=0"))
GSS$gunimp <- with(GSS, recode(gunimp, "1=3; 2=2; 3=1; 4=0")) # Flip the scale; create naturally occurring zero (i.e. the not important at all).

GSS$guninfod <- with(GSS, recode(guninfo, "1:2=1; 3:4=0"))
GSS$gunfirmd <- with(GSS, recode(gunfirm, "1:2=0; 3:4=1")) # The unlikely to changes are 1

GSS$owngund <- with(GSS, recode(owngun, "1=1; 2=0; 3=NA"))

GSS$learnguninternet <- with(GSS, recode(gun12, "1=1; 2=0"))

GSS$gunsales <- with(GSS, recode(gunsales, "1=5; 2=4; 3=3; 4=2; 5=1"))
GSS$gunsalesd <- with(GSS, recode(gunsales, "4:5=1;1:3=0"))
GSS$gunsdrug <- with(GSS, recode(gunsdrug, "1=1; 2=-1; 3=0"))
GSS$gunsdrugd <- with(GSS, recode(gunsdrug, "1=1; -1=0; 0=0"))
GSS$semiguns <- with(GSS, recode(semiguns, "1=0; 2=1"))
GSS$guns911 <- with(GSS, recode(guns911, "1=1; 2=0"))
GSS$gunsdrnk <- with(GSS, recode(gunsdrnk, "1=1;2=0")) # Would support state laws making it illegal to carry a firearm while drunk.

GSS$hgunlaw <- with(GSS, recode(hgunlaw, "1=1; 2=0"))
GSS$hguncrim <- with(GSS, recode(hguncrim, "1=1; 2=-1; 3=0"))
GSS$crimup <- with(GSS, recode(crimup, "1=1; 2=0"))
GSS$crimdown <- with(GSS, recode(crimdown, "1=1; 2=0"))

GSS <- ddply(GSS, c("year"), transform, z.age = arm::rescale(age))
GSS <- ddply(GSS, c("year"), transform, z.educ = arm::rescale(educ))

GSS$collegeed <- with(GSS, recode(degree, "0:2=0; 3:4=1"))
GSS$hsed <- with(GSS, recode(degree, "0=0; 1:4=1"))
GSS$graded <- with(GSS, recode(degree, "0:3=0; 4=1"))

GSS$white <- with(GSS, recode(race, "1=1; 2:3=0"))
GSS$black <- with(GSS, recode(race, "1=0; 2=1; 3=0"))
GSS$otherrace <- with(GSS, recode(race, "1:2=0; 3=1"))

GSS <- ddply(GSS, c("year"), transform, z.rincome = arm::rescale(rincome))
GSS <- ddply(GSS, c("year"), transform, z.rincom06 = arm::rescale(rincom06))

GSS$partyid <- with(GSS, recode(partyid, "7=NA")) # Make "Other party, refuse to say" into NAs
GSS <- ddply(GSS, c("year"), transform, z.partyid = arm::rescale(partyid))

GSS$female <- with(GSS, recode(sex, "1=0; 2=1"))
GSS$hispanic <- with(GSS, ifelse(ethnic == 17 | ethnic == 22 | ethnic == 25 | ethnic == 32 | ethnic == 38, 1, 0))

YearMeans <- aggregate(cbind(gunlaw) ~ partyid + year, GSS, mean)
YearMeans$gunlaw <- YearMeans$gunlaw *100
YearMeans$gunlaw <- round(YearMeans$gunlaw, 2)

WideYearMeans <- dcast(data = YearMeans, formula = partyid ~ year)
WideYearMeans$partyid <- with(WideYearMeans, recode(partyid, "0='Strong Dem.'; 1='Not Strong Dem.'; 2='Indep., lean Dem.'; 3='Independent'; 4='Indep., lean GOP'; 5='Not Strong GOP'; 6='Strong GOP'"))
stargazer(WideYearMeans[,1:18], type="html", summary=FALSE, rownames=FALSE, font.size="scriptsize", digits=NA, style="ajps")
stargazer(WideYearMeans[,c("partyid", "1996", "1998", "2000", "2002", "2004", "2006", "2008", "2010", "2012", "2014")], type="html", summary=FALSE, rownames=FALSE, font.size="scriptsize", digits=NA, style="ajps")

OtherMeans2006 <- aggregate(cbind(gunsalesd, gunsdrugd, semiguns, gunsdrnk, guns911) ~ partyid, GSS, mean)
OtherMeans2006$gunsalesd <- OtherMeans2006$gunsalesd * 100
OtherMeans2006$gunsdrugd <- OtherMeans2006$gunsdrugd * 100
OtherMeans2006$semiguns <- OtherMeans2006$semiguns * 100
OtherMeans2006$gunsdrnk <- OtherMeans2006$gunsdrnk * 100
OtherMeans2006$guns911 <- OtherMeans2006$guns911 * 100
stargazer(OtherMeans2006, type="html", summary=FALSE, rownames=FALSE, font.size="scriptsize", digits=2, style="ajps")

## Estimate some models.
########################

summary(M1 <- glmer(gunlaw ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | year) + (1 + z.partyid | regioncondensed), data=subset(GSS), family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
summary(M1a <- glmer(gunlaw ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | year) + (1 + z.partyid | regioncondensed), data=subset(GSS, year <= 1993), family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
summary(M1b <- glmer(gunlaw ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | year) + (1 + z.partyid | regioncondensed), data=subset(GSS, year > 1993), family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))

ZeligGSS <- with(GSS, data.frame(gunlaw, z.age, collegeed, female, black, otherrace, partyid,  owngund, year, regioncondensed))
ZeligGSS <- na.omit(ZeligGSS)

Z1 <- zelig(gunlaw ~ z.age + collegeed + female + black + otherrace + partyid + owngund + tag(1 + partyid | year) + tag(1 + partyid | regioncondensed), data=ZeligGSS, model="logit.mixed")

Z1.sgopgunwomen <- setx(Z1, partyid = max(partyid), z.age=0, collegeed=median(collegeed), female=1, black=median(black), otherrace=median(otherrace), owngund=1)
Z1.sgopgunmen <- setx(Z1, partyid = max(partyid), z.age=0, collegeed=median(collegeed), female=0, black=median(black), otherrace=median(otherrace), owngund=1)
Z1.sgopwomen <- setx(Z1, partyid = max(partyid), z.age=0, collegeed=median(collegeed), female=1, black=median(black), otherrace=median(otherrace), owngund=0)
Z1.sgopmen <- setx(Z1, partyid = max(partyid), z.age=0, collegeed=median(collegeed), female=0, black=median(black), otherrace=median(otherrace), owngund=0)
Z1.sgopgunmencollege <- setx(Z1, partyid = max(partyid), z.age=0, collegeed=1, female=0, black=median(black), otherrace=median(otherrace), owngund=1)
Z1.sim.sgopgunwomen <- sim(Z1, x = Z1.sgopgunwomen)
Z1.sim.sgopgunmen <- sim(Z1, x = Z1.sgopgunmen)
Z1.sim.sgopwomen <- sim(Z1, x = Z1.sgopwomen)
Z1.sim.sgopmen <- sim(Z1, x = Z1.sgopmen)
Z1.sim.sgopgunmencollege <- sim(Z1, x = Z1.sgopgunmencollege)
summary(Z1.sim.sgopgunwomen)
summary(Z1.sim.sgopgunmen)
summary(Z1.sim.sgopwomen)
summary(Z1.sim.sgopmen)
summary(Z1.sim.sgopgunmencollege)

ggCaterpillar(lme4::ranef(M1,condVar=TRUE), QQ=FALSE, likeDotplot=TRUE, reorder=FALSE)[["regioncondensed"]]

M1ranef <- lme4::ranef(M1,condVar=TRUE)
colnames(M1ranef[["year"]]) <- c("Intercept","Party ID (D to R)") 
ggCaterpillar(M1ranef, QQ=FALSE, likeDotplot=TRUE, reorder=FALSE)[["year"]]
ggsave(file="gss-m1-ranef.png")

M1aranef <- lme4::ranef(M1a,condVar=TRUE)
colnames(M1aranef[["regioncondensed"]]) <- c("Intercept","Party ID (D to R)") 
ggCaterpillar(M1aranef, QQ=FALSE, likeDotplot=TRUE, reorder=FALSE)[["regioncondensed"]]
ggsave(file="gss-m1a-ranef.png")

M1branef <- lme4::ranef(M1b,condVar=TRUE)
colnames(M1branef[["regioncondensed"]]) <- c("Intercept","Party ID (D to R)") 
ggCaterpillar(M1branef, QQ=FALSE, likeDotplot=TRUE, reorder=FALSE)[["regioncondensed"]]
ggsave(file="gss-m1b-ranef.png")


summary(M2 <- glmer(gunsalesd ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | regioncondensed), data=GSS, family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
summary(M3 <- glmer(gunsdrugd ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | regioncondensed), data=GSS, family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
summary(M4 <- glmer(semiguns ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | regioncondensed), data=GSS, family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
summary(M5 <- glmer(gunsdrnk ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | regioncondensed), data=GSS, family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))
summary(M6 <- glmer(guns911 ~ z.age + collegeed + female + black + otherrace + z.partyid + owngund + (1 + z.partyid | regioncondensed), data=GSS, family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))))

Table <- stargazer(M1, M2, M3, M4, M5, M6, type="html", style = "ajps", title="Mixed Effects Models of Attitudes toward Gun Control",
	covariate.labels=c("Age", "College Educated", "Female", "Black", "Other Race (Not White)", "Party ID (D to R)", "Gun in the Household"),
	dep.var.labels=c("Require Police Permit?", "Background Check for Private Sales?", "Tougher Penalties than Selling Drugs?", "Limit Semi-Automatics to Police/Military?", "Illegal to Carry a Gun while Drunk?", "Tougher Gun Control Laws after 9/11?"), dep.var.labels.include = TRUE, 
	model.names=FALSE, omit.stat=c("aic","ll","bic")
)
