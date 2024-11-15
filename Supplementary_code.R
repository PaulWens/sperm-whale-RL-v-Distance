Sys.setenv(TZ='GMT') 

## Load R packages and data 
library(mgcv)
library(scales)

load("Supplementary_data.Rdata")


# Variable names in code vs. manuscript (see also Table 2)
# mtab$NF_prob          # pNFA      
# mtab$BUZZ_pres        # buzz          
# mtab$fluken           # fluke
# mtab$state            # state
# mtab$solarnoon        # solarnoon
# mtab$depth_first      # depth
# mtab$Bathy_First      # depth_bathy
# mtab$DW_2h            # BF_2h       
# mtab$UPAS_2h          # UPAS_2h       
# mtab$TS               # TS
# mtab$NS               # NS
# mtab$SEL_max_incr     # SELmax
# mtab$SL_hrange_min    # DISTmin  
# mtab$prev_hrange      # DISTmin_prev     
# mtab$prev_hrange_post # DISTmin_post   
# mtab$DW_2h_dur        # BF_2h_dur


## Fit GAMMs
# pNFA

fit_baseline <- gamm(NF_prob ~ s(Bathy_First, bs="ts") + s(solarnoon, bs="cc",k=7) + 
                       DW_2h + UPAS_2h + TS,
                     data=mtab[fitBool_baseline,], random=list(ind=~1), 
                     correlation=corARMA(p=1,q=1), family=quasibinomial(link="logit"))
summary(fit_baseline$lme) # get corARMA coefs
fit_baseline <- gamm(NF_prob ~ s(Bathy_First, bs="ts") + s(solarnoon, bs="cc",k=7) + DW_2h, 
                     data=mtab[fitBool_baseline,], random=list(ind=~1), 
                     correlation=corARMA(value=c(0.7918997, -0.4155180), p=1, q=1, fixed=T),
                     family=quasibinomial(link="logit")) 
fit_exposure <- gamm(NF_prob ~ s(solarnoon, bs="cc", k=7) + DW_2h +
                       NS + DW_2h_dur + s(SEL_max_incr, bs="ts", k=5) + s(SL_hrange_min, bs="ts", k=5) +
                       ti(SEL_max_incr, SL_hrange_min, bs="ts")  +
                       s(prev_hrange, bs="ts", k=5) + s(prev_hrange_post, bs="ts", k=5), 
                     data=mtab[fitBool_exposures,], random=list(ind=~1), 
                     correlation=corARMA(value=c(0.7918997, -0.4155180), p=1, q=1, fixed=T),
                     family=quasibinomial(link="logit"))
fit_exposure <- gamm(NF_prob ~ s(solarnoon, bs="cc", k=7) + DW_2h +
                       s(SEL_max_incr, bs="ts", k=5) + s(SL_hrange_min, bs="ts", k=5) +
                       ti(SEL_max_incr, SL_hrange_min, bs="ts")  +
                       s(prev_hrange, bs="ts", k=5) + s(prev_hrange_post, bs="ts", k=5), 
                     data=mtab[fitBool_exposures,], random=list(ind=~1), 
                     correlation=corARMA(value=c(0.7918997, -0.4155180), p=1, q=1, fixed=T),
                     family=quasibinomial(link="logit"))
anova(fit_exposure$gam)
summary(fit_exposure$gam)
plot(fit_exposure$gam,rug=T,shade=T,scheme=2,too.far=0.1)
acf(residuals(fit_exposure$lme, type="normalized"), lag=20, main="")  

# ROC curve
y <- as.numeric(mtab$NF[fitBool_exposures]) 
preds <- predict(fit_exposure$gam, type="response")
TH_vec <- seq(0,1,0.001)
TPR <- NA # Sensitivity 
TNR <- NA # Specificity
for(j in 1:length(TH_vec)) {
  TP <- sum(as.numeric(y==1 & preds > TH_vec[j]))
  FP <- sum(as.numeric(y==0 & preds > TH_vec[j]))
  FN <- sum(as.numeric(y==1 & preds <= TH_vec[j]))
  TN <- sum(as.numeric(y==0 & preds <= TH_vec[j]))
  TPR[j] <- TP/(TP + FN)
  TNR[j] <- TN/(FP + TN)
}
plot(1-TNR, TPR, type="l", xlim=c(0,1), ylim=c(0,1), lwd=2.5)
abline(0,1,col="grey",lty=2)


# buzz
# only includes foraging states
fitBool_baseline.forage <- fitBool_baseline & mtab$state!="Surface" & mtab$state!="Resting" & mtab$state!="NF Active"
fitBool_exposures.forage <- fitBool_exposures & mtab$state!="Surface" & mtab$state!="Resting" & mtab$state!="NF Active"

fit_baseline <- gamm(BUZZ_pres ~ s(Bathy_First, bs="ts") + s(solarnoon, bs="cc", k=7) +
                       s(depth_first, bs="ts") + DW_2h + UPAS_2h + TS,
                     mtab[fitBool_baseline.forage,], random=list(ind=~1),
                     correlation=corARMA(p=1,q=1), family=quasibinomial(link="logit"))
summary(fit_baseline$lme) # get corARMA coefs
fit_baseline <- gamm(BUZZ_pres ~ s(Bathy_First, bs="ts") + s(solarnoon, bs="cc", k=7) +
                       s(depth_first, bs="ts"),
                     mtab[fitBool_baseline.forage,], random=list(ind=~1), 
                     correlation=corARMA(value=c(0.8497081, -0.7364822), p=1, q=1, fixed=T),
                     family=quasibinomial(link="logit"))
fit_exposure <- gamm(BUZZ_pres ~ s(solarnoon, bs="cc", k=7) + s(depth_first, bs="ts") + 
                       NS + DW_2h_dur + s(SEL_max_incr, bs="ts", k=5) + s(SL_hrange_min, bs="ts", k=5) + 
                       ti(SEL_max_incr, SL_hrange_min, bs="ts")  +
                       s(prev_hrange, bs="ts", k=5) + s(prev_hrange_post, bs="ts", k=5), 
                     data=mtab[fitBool_exposures.forage,], random=list(ind=~1),
                     correlation=corARMA(value=c(0.8572635, -0.7479970), p=1, q=1, fixed=T),
                     family=quasibinomial(link="logit"))
fit_exposure <- gamm(BUZZ_pres ~ s(solarnoon, bs="cc", k=7) + s(depth_first, bs="ts") + 
                       s(SEL_max_incr, bs="ts", k=5) + s(SL_hrange_min, bs="ts", k=5) + 
                       ti(SEL_max_incr, SL_hrange_min, bs="ts")  +
                       s(prev_hrange, bs="ts", k=5) + s(prev_hrange_post, bs="ts", k=5), 
                     data=mtab[fitBool_exposures.forage,], random=list(ind=~1),
                     correlation=corARMA(value=c(0.8572635, -0.7479970), p=1, q=1, fixed=T),
                     family=quasibinomial(link="logit"))
anova(fit_exposure$gam)
summary(fit_exposure$gam)
plot(fit_exposure$gam,rug=T,shade=T,scheme=2,too.far=0.1)
acf(residuals(fit_exposure$lme, type="normalized"), lag=20, main="") 

# ROC curve
y <- as.numeric(fit_exposure$gam$y) 
preds <- predict(fit_exposure$gam, type="response") 
TH_vec <- seq(0,max(y),0.01)
TPR <- NA # Sensitivity 
TNR <- NA # Specificity 
for(j in 1:length(TH_vec)) {
  TP <- sum(as.numeric(y==1 & preds > TH_vec[j]))
  FP <- sum(as.numeric(y==0 & preds > TH_vec[j]))
  FN <- sum(as.numeric(y==1 & preds <= TH_vec[j]))
  TN <- sum(as.numeric(y==0 & preds <= TH_vec[j]))
  TPR[j] <- TP/(TP + FN)
  TNR[j] <- TN/(FP + TN)
}
plot(1-TNR, TPR, type="l", xlim=c(0,1), ylim=c(0,1), lwd=2.5, las=1)
abline(0,1,col="grey",lty=2)

# fluke
# excludes surface state
fitBool_baseline.nosurf <- fitBool_baseline & mtab$state!="Surface"
fitBool_exposures.nosurf <- fitBool_exposures & mtab$state!="Surface"

fit_baseline <- gamm(fluken ~ state + s(Bathy_First, bs="ts") + s(solarnoon, bs="cc",k=7) + 
                       DW_2h + UPAS_2h +TS, 
                     mtab[fitBool_baseline.nosurf,], random=list(ind=~1), 
                     correlation=corAR1(), family=quasipoisson(link="log"))
fit_baseline <- gamm(fluken ~ state + s(Bathy_First, bs="ts") + s(solarnoon, bs="cc",k=7) + DW_2h, 
                     mtab[fitBool_baseline.nosurf,], random=list(ind=~1), 
                     correlation=corAR1(), family=quasipoisson(link="log"))
fit_exposure <- gamm(fluken ~ state + DW_2h + NS + DW_2h_dur + 
                       s(SEL_max_incr, bs="ts", k=5) + s(SL_hrange_min, bs="ts", k=5) +
                       SEL_max_incr:SL_hrange_min +
                       s(prev_hrange, bs="ts", k=5) + s(prev_hrange_post, bs="ts", k=5),
                     data=mtab[fitBool_exposures,], random=list(ind=~1), 
                     correlation=corAR1(), family=quasipoisson(link="log"))
fit_exposure <- gamm(fluken ~ state + DW_2h + 
                       s(SEL_max_incr, bs="ts", k=5) + s(SL_hrange_min, bs="ts", k=5) +
                       SEL_max_incr:SL_hrange_min +
                       s(prev_hrange, bs="ts", k=5) + s(prev_hrange_post, bs="ts", k=5),
                     data=mtab[fitBool_exposures,], random=list(ind=~1), 
                     correlation=corAR1(), family=quasipoisson(link="log"))
anova(fit_exposure$gam)
summary(fit_exposure$gam)
plot(fit_exposure$gam,rug=T,shade=T,scheme=2,too.far=0.1)
acf(residuals(fit_exposure$lme, type="normalized"), lag=20, main="") 

# Predicted vs observed
m <- fit_exposure
obs <- m$gam$model$fluken
pred.mu1 <- rep(0,nrow(m$gam$model)) 
for(i in 1:nrow(m$gam$model)){
  preddata <- data.frame(state=m$gam$model$state[i], DW_2h=m$gam$model$DW_2h[i],
                         SEL_max_incr=0, 
                         SL_hrange_min=0,
                         prev_hrange=m$gam$model$prev_hrange[i],
                         prev_hrange_post=m$gam$model$prev_hrange_post[i])
  tObj <- predict(m$gam, newdata=preddata, type="response", se.fit=T)
  pred.mu1[i] <- tObj$fit
}
sts <- mtab$state[fitBool_exposures.nosurf]
plot(obs[sts=="Resting"], pred.mu1[sts=="Resting"], col = alpha("#F0E442", 0.1), pch=16, 
     xlim=c(0,15), xlab = expression(paste("Observed fluke rate (min"^-1,")")), 
     ylim=c(0,6.5), ylab = expression(paste("Predicted fluke rate (min"^-1,")")))
abline(a=0,b=1,col="grey",lty=2)
points(obs[sts=="LRS"], pred.mu1[sts=="LRS"], col = alpha("#009E73", 0.1), pch=16)
points(obs[sts=="Ascent"]+0.2, pred.mu1[sts=="Ascent"], col = alpha("#56B4E9", 0.1), pch=16)  
points(obs[sts=="Descent"]-0.2, pred.mu1[sts=="Descent"], col = alpha("#E69F00", 0.1), pch=16)
points(obs[sts=="NF Active"], pred.mu1[sts=="NF Active"], col = alpha("#CC79A7", 0.1), pch=16)
