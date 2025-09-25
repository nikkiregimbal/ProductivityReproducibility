library(ggplot2)
library(dplyr)
library(vegan)
library(lme4)
library(DHARMa)

#Read in dataset
df <- read.csv("./02_outdata/revisedSep18-25_backswimmer_development_dispersal.csv")

# Temperature on Survival -------------------------------------------------

#Generalized linear mixed effects model assessing survival probability
survival <- glmer(Survival ~ Treatment + (1|Block/Mesocosm), data = df, family = binomial(link="logit"))
  #getting singularity warning because Mesocosm is nested in Block - must ignore this as both random effects are important to include
  summary(survival)
  drop1(survival, test = 'Chisq') #p < 0.001, chi = 21.01

#testing model assumptions
  #model assumptions met
testOutliers(survival) #outlier test
testDispersion(survival) #overdispersion test

    #Table showing the survival count and probability for each mesocosm 
    mes <- as.data.frame(df %>%
      group_by(Treatment, Mesocosm) %>%
      summarise(surv_count = sum(Survival == 1)) %>%
      mutate(total_ind = 10, survival_prob = surv_count / total_ind))
  

#Bootstrapping to generate confidence intervals (CI) and fitted survival probability
    nd <- data.frame(unique(subset(df, select =c(Treatment))), "Mesocosm"=NA, "Block" = NA)
    
    mySumm_mod <- function(.){
      predict(., newdata=nd, re.form=~0)
    }
    
    boot <- bootMer(survival, mySumm_mod, nsim = 200, use.u=F, type="parametric")
    
    sumBoot <- function(merBoot) {
      return(
        data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))), #predicted fit
                   lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))), #lower CI 
                   upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))) #upper CI
        )
      )
    }
    
    ## Use the function to estimate the model fit and 95% CI
    survival_est <- sumBoot(boot)
    
    # A GLM with a binomial family uses a logit transformation; in order to get probability estimates, we have to untransform the values
    ## This function untransforms the values from logit units to probability
    p_survival <- exp(survival_est) / (1 + exp(survival_est))


survival_result <- cbind(nd, p_survival);survival_result

survival_result$Treatment <- as.integer(as.character(survival_result$Treatment))

#plotting temp on survival
survival_plot <- ggplot(survival_result, aes(x = Treatment, y = fit))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3, col=NA)+
  geom_point(data = mes, aes(x = Treatment, y = survival_prob),  position=position_jitter(h=0.03,w=0.3))+
  scale_x_continuous("Temperature (°C)")+
  scale_y_continuous("Survival Probability", limits = c(0,1,0.1))+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
survival_plot

#save this plot into figs folder
ggsave("./03_figs/survival_plot.png", plot = survival_plot)

# Temperature on Molt Rate ---------------------------------------------------------------
#Does temperature affect the molting rate (total molts/days until adulthood or mortality)

#Use log molt rate because it is normally distributed
hist(df$logMoltRate)

#new df omits NAs from log molt rate (omitting individuals that died before ever molting)
df1 <- df[!is.na(df$logMoltRate),]

#Linear mixed effects model for temperature on molt rate
moltrate <- lmer(logMoltRate ~ Treatment + (1|Block/Mesocosm), data=df1)
drop1(moltrate, test = 'Chisq') #p < 0.001, chi = 55.56
summary(moltrate)

#testing model assumptions
  #model assumptions met
testOutliers(moltrate) #outlier test
testDispersion(moltrate) #overdispersion test

#Confidence interval and predicted fits
    moltc <- babydf %>%
      group_by(Treatment) %>%
      summarise(mean = mean(logMoltRate), sd = sd(logMoltRate), N = n()) %>%
      mutate(se = sd / sqrt(N),
             lwr = mean - qt(1 - (0.05 / 2), N - 1) * se,
             upr = mean + qt(1 - (0.05 / 2), N - 1) * se)
    
    moltc <- moltc %>%
      mutate(meant = exp(mean), sdt = exp(sd), set = exp(se), lwrt = exp(lwr), uprt = exp(upr))


#plotting temperature on molt rate
molt_plot <- ggplot(moltc, aes(x = Treatment, y = meant))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwrt, ymax=uprt), alpha=0.5, col=NA)+
  #labs(tag="A")+
  scale_x_continuous("Temperature (C)")+
  scale_y_continuous("Molt Rate (# molts / day)", limits = c(0.05,0.2,0.05))+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
molt_plot 

#saving molt rate plot to figs folder
ggsave("./03_figs/temp_moltrate_plot.png", plot = molt_plot)


#bootstrapping to find the fitted probability of molt rate as a function of temp 
nd <- df1

mySumm_mod <- function(.){
  predict(., newdata=nd, re.form=~0)
}

boot <- bootMer(moltrate, mySumm_mod, nsim = 200, use.u=F, type="parametric")

sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

## Use the function to estimate the model fit and 95% CI
molt_est <- sumBoot(boot)

mr <- cbind(nd, molt_est);mr
mr$Treatment <- as.integer(as.character(mr$Treatment))

#molt rate predicted fit plot
logmolt_plot <- ggplot(mr, aes(x = Treatment, y = fit))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.5, col=NA)+
  labs(tag="A")+
  scale_x_continuous("Temperature")+
  scale_y_continuous("log(Molt Rate)")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
logmolt_plot 

#saving molt rate predicted fit plot to figs folder
ggsave("./03_figs/temp_logmoltrate_plot.png", plot = logmolt_plot)


# Molt Rate on Survival ---------------------------------------------

#new df omits NAs from log molt rate (omitting individuals that died before ever molting)
df1 <- df[!is.na(df$logMoltRate),] #code for this subsetted dataset is repeated from the last section

#Generalized linear mixed effects model for molt rate on survival probability
molt_survival1<- glmer(Survival ~ MoltRate + (1|Block/Mesocosm), family = binomial(link = 'logit'), data = df1)
  summary(molt_survival1)
  drop1(molt_survival1, test = "Chisq") #p < 0.001, chi = 16.65

#bootstrapping to generate predicted fits and CIs
  nd <- df1
    
    data.frame(unique(subset(df1, select =c(MoltRate))), "Mesocosm"=NA, "Block" = NA)
  
  mySumm_mod <- function(.){
    predict(., newdata=nd, re.form=~0)
  }
  
  boot <- bootMer(molt_survival1, mySumm_mod, nsim = 200, use.u=F, type="parametric")
  
  sumBoot <- function(merBoot) {
    return(
      data.frame(survival_fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))), #predicted fit
                 survival_lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))), #lower CI
                 survival_upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE))) #upper CI
      )
    )
  }
  
  ## Use the function to estimate the model fit and 95% CI
  ms_est <- sumBoot(boot)
  
  pmoltsurv <- exp(ms_est) / (1 + exp(ms_est))
  
  mr_s <- cbind(nd, pmoltsurv);mr_s
  mr_s$Treatment <- as.integer(as.character(mr$Treatment))

#plotting the impact of molt rate on survival  
molt_survival_plot <- ggplot(mr_s, aes(x = MoltRate, y = survival_fit))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=survival_lwr, ymax=survival_upr), alpha=0.5, col=NA)+
  scale_x_continuous("Molt Rate")+
  scale_y_continuous("Survival Probability", limits = c(0,0.6), breaks = seq(0,0.6,0.2))+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
molt_survival_plot

#saving molt rate on survival plot to figs folder
ggsave("./03_figs/moltrate_survival_plot.png", plot = molt_survival_plot)


# Temperature on Body condition -------------------------------------------
#Body condition is only assessed for individuals that survived until adulthood

#creating df of only individuals that survived to adulthood
survived <- subset(df, AdulthoodDay > 0) #omitting NA values where the backswimmer never became an adult
head(survived)

  #Table of number of survivors by temp treatment
  t <- as.data.frame(table(survived$Treatment))

#Body Condition (BC) immediately after adult emergence
survived$BC <- survived$Mass1/survived$Thorax_Width1

#Linear mixed effects model for temperature on body condition 
condition <- lmer(BC ~ Treatment + (1|Block/Mesocosm), data = survived)
  summary(condition)
  drop1(condition, test=  'Chisq') #not significant


