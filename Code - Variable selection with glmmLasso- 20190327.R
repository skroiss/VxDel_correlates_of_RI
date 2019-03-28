# Purpose: Perform variable selection for correlates of RI analysis
# Inputs:  Curated WUENIC DTP3 & potential correlates
# Output:  -
# Author:  Steve Kroiss (skroiss@idmod.org)
# Date:    March 2018

# Clear objects from R workspace
rm(list=ls())  



#############################
# Load libraries
#############################

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(glmmLasso)
options(tibble.width = Inf)



#############################
# Load data
#############################

# Load WUENIC data
load('Curated WUENIC estimates and correlates of RI - 20190115.RData')



#############################
# Prep data for analysis 
#############################

# Select variables with low missingness & eliminate rows with NAs
ri_corr_select <- ri_corr %>%
  select(country,
         iso,
         iso_id,
         year,
         wuenic_dtp3,
         adolescent_fert,
         controlofcorruption,
         fertility,
         governmenteffectiveness,
         ln_gdp,
         ln_healthexp,
         politicalstability,
         regulatoryquality,
         ruleoflaw,
         rural_pop,
         voiceandaccountability) %>%
  # Omit NAs
  na.omit %>%
  mutate(country_fctr = factor(country))


# Center & scale data
ri_corr_select_centered_scaled <- ri_corr_select
ri_corr_select_centered_scaled[,c(6:16)] <- scale(ri_corr_select[,c(6:16)],
                                        center=TRUE,
                                        scale=TRUE)
ri_corr_select_centered_scaled <- data.frame(ri_corr_select_centered_scaled)

# # Export data
# save(ri_corr_select,
#      ri_corr_select_centered_scaled,
#      file = 'RI covariates - selected for low missingness - NAs omitted - 20190327.RData')


#############################
# Check on basic lme4 model
#############################

# Run linear regression model
m_lmer_null <- lmer(wuenic_dtp3 ~ (1|country_fctr), 
                   data=ri_corr_select_centered_scaled)
m_lmer_int <- lmer(wuenic_dtp3 ~ adolescent_fert + 
                     controlofcorruption + 
                     fertility + 
                     governmenteffectiveness + 
                     ln_gdp + 
                     ln_healthexp + 
                     politicalstability + 
                     regulatoryquality + 
                     ruleoflaw + 
                     rural_pop + 
                     voiceandaccountability + 
                     (1|country_fctr), 
                   data=ri_corr_select_centered_scaled)
# Print summaries
summary(m_lmer_null)
summary(m_lmer_int)

# Check BICs
BIC(m_lmer_null)
BIC(m_lmer_int)

# Likelihood ratio test suggests full model is preferred over null model
anova(m_lmer_null, m_lmer_int)





#############################
# Run an unpenalized LASSO
#############################

# Unpenalized LASSO
m_lasso_uni <- glmmLasso(wuenic_dtp3 ~ adolescent_fert,
                         rnd = list(country_fctr = ~1),
                         data = ri_corr_select_centered_scaled,
                         family = gaussian(link="identity"),
                         lambda=0)
m_lasso_int <- glmmLasso(wuenic_dtp3 ~ adolescent_fert + 
                              controlofcorruption + 
                              fertility + 
                              governmenteffectiveness + 
                              ln_gdp + 
                              ln_healthexp + 
                              politicalstability + 
                              regulatoryquality + 
                              ruleoflaw + 
                              rural_pop + 
                              voiceandaccountability,
                            rnd = list(country_fctr = ~1),
                            data = ri_corr_select_centered_scaled,
                            family = gaussian(link="identity"),
                            lambda=0)
summary(m_lasso_int)



m_lasso_uni$bic
m_lasso_int$bic




################## First Simple Method ############################################
## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda

## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
lambda <- seq(500, 0, by=-5)
BIC_vec<-rep(Inf, length(lambda))
family = gaussian(link="identity")

## first fit good starting model
library(MASS)
library(nlme)
PQL<-glmmPQL(wuenic_dtp3~1, 
             random = ~1|country_fctr, 
             family=family, 
             data=ri_corr_select_centered_scaled)
Delta.start<-c(as.numeric(PQL$coef$fixed),
               rep(0,11),
               as.numeric(t(PQL$coef$random$country_fctr)))
Q.start<-as.numeric(VarCorr(PQL)[1,1])

# Run glmmLasso
for(j in 1:length(lambda)){
  print(paste("Iteration ", j,sep=""))
  
  glm1 <- try(glmmLasso(wuenic_dtp3 ~ adolescent_fert + 
                          controlofcorruption + 
                          fertility + 
                          governmenteffectiveness + 
                          ln_gdp + 
                          ln_healthexp + 
                          politicalstability + 
                          regulatoryquality + 
                          ruleoflaw + 
                          rural_pop + 
                          voiceandaccountability,
                        rnd = list(country_fctr = ~1),
                        family = family,
                        data = ri_corr_select_centered_scaled, 
                        lambda=lambda[j],
                        switch.NR=T,
                        final.re=TRUE,
                        control=list(start=Delta.start, q_start=Q.start)), 
              silent=TRUE) 
  
  if(class(glm1)!="try-error")
  {  
    BIC_vec[j]<-glm1$bic
  }
  
}

opt<-which.min(BIC_vec)

glm1_final <- glmmLasso(wuenic_dtp3 ~ adolescent_fert + 
                          controlofcorruption + 
                          fertility + 
                          governmenteffectiveness + 
                          ln_gdp + 
                          ln_healthexp + 
                          politicalstability + 
                          regulatoryquality + 
                          ruleoflaw + 
                          rural_pop + 
                          voiceandaccountability,
                        rnd = list(country_fctr = ~1),
                        family = family,
                        data = ri_corr_select_centered_scaled,
                        lambda=lambda[opt],
                        switch.NR=F,
                        final.re=TRUE,
                        control=list(start=Delta.start,q_start=Q.start))
summary(glm1_final)



#################################
# More elegant method to identify optimal tuning parameter lambda for LASSO
#################################

## Idea: start with big lambda and use the estimates of the previous fit (BUT: before
## the final re-estimation Fisher scoring is performed!) as starting values for the next fit;
## make sure, that your lambda sequence starts at a value big enough such that all covariates are
## shrinked to zero;
## Source: https://rdrr.io/cran/glmmLasso/src/demo/glmmLasso-soccer.r

## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
lambda <- seq(500, 0, by=-5)
BIC_vec<-rep(Inf, length(lambda))
family = gaussian(link="identity")

# Specify starting values for the very first fit; pay attention that Delta.start has suitable length! 
# The argument start of the control list is a vector of starting values for all model parameters, 
# i.e. all fixed effects (usually including the intercept) followed by all random effects 
# (and thus depending on the length of your gouping factor and on the structure of your random terms)
# suppose your grouping factor id contains 80 individuals, then the start vector would have to be of length 
# 1(intercept)+10(fixed effects)+80(random effects)=91.
## More details on specifying starting values: https://stats.stackexchange.com/questions/70753/correct-estimation-of-arguments-for-glmmlasso-function
Delta.start<-as.matrix(t(rep(0, 1+11+length(unique(ri_corr_select_centered_scaled$country_fctr)))))
Q.start<-0.1  

# Run 
for(j in 1:length(lambda)){
  print(paste("Iteration ", j,sep=""))
  
  glmm_lasso_int <- glmmLasso(wuenic_dtp3 ~ adolescent_fert + 
                                controlofcorruption + 
                                fertility + 
                                governmenteffectiveness + 
                                ln_gdp + 
                                ln_healthexp + 
                                politicalstability + 
                                regulatoryquality + 
                                ruleoflaw + 
                                rural_pop + 
                                voiceandaccountability,
                              rnd = list(country_fctr = ~1),
                              data = ri_corr_select_centered_scaled,
                              family = gaussian(link="identity"),
                              lambda=lambda[j], 
                              switch.NR=F, 
                              final.re=TRUE,
                              control = list(start=Delta.start[j,], q_start=Q.start[j]))
  # Check that initial run has a large enough lambda to force estimates to 0
  # summary(glmm_lasso_int)  # when j = 1; yup, looks good

  print(colnames(glmm_lasso_int$Deltamatrix)[2:12][glmm_lasso_int$Deltamatrix[glmm_lasso_int$conv.step,2:12] != 0])
  BIC_vec[j] <- glmm_lasso_int$bic
  Delta.start <- rbind(Delta.start,
                       glmm_lasso_int$Deltamatrix[glmm_lasso_int$conv.step,])
  Q.start <- c(Q.start,
               glmm_lasso_int$Q_long[[glmm_lasso_int$conv.step + 1]])
}
# ID optimal model
opt_int <- which.min(BIC_vec)

# Run final model
glmm_lasso_int_final <- glmmLasso(wuenic_dtp3 ~ 1 + adolescent_fert + 
                                    controlofcorruption + 
                                    fertility + 
                                    governmenteffectiveness + 
                                    ln_gdp + 
                                    ln_healthexp + 
                                    politicalstability + 
                                    regulatoryquality + 
                                    ruleoflaw + 
                                    rural_pop + 
                                    voiceandaccountability,
                                  rnd = list(country_fctr = ~1),
                                  data = ri_corr_select_centered_scaled,
                                  lambda=lambda[opt_int],
                                  switch.NR=F,
                                  final.re=TRUE,
                                  control = list(start = Delta.start[opt_int,], q_start = Q.start[opt_int]))
summary(glmm_lasso_int_final)


