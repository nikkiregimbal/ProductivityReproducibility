library(ggplot2)
library(dplyr)
library(vegan)
library(lme4)
library(partR2)

df <- read.csv("./02_outdata/revisedSep18-25_backswimmer_development_dispersal.csv")

# Temperature on Survival -------------------------------------------------

survival <- glmer(Survival ~ Treatment + (1|Block/Mesocosm), data = df, family = binomial(link="logit"))
  #getting singularity warning bc Mesocosm is nested in Block - must ignore this bc both random effects are important to include
summary(survival)
drop1(survival, test = 'Chisq') #p < 0.001, chi = 19.892


mes <- as.data.frame(babydf %>%
  group_by(Treatment, Mesocosm) %>%
  summarise(surv_count = sum(Survival == 1)) %>%
  mutate(total_ind = 10, survival_prob = surv_count / total_ind))

nd <- data.frame(unique(subset(babydf, select =c(Treatment))), "Mesocosm"=NA, "Block" = NA)

mySumm_mod <- function(.){
  predict(., newdata=nd, re.form=~0)
}

boot <- bootMer(survival2, mySumm_mod, nsim = 200, use.u=F, type="parametric")

sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

## Use the function to estimate the model fit and 95% CI
survival_est <- sumBoot(boot)

# A GLM with a binomial family uses a logit transformation; in order to get probability estimates, we have to untransform the values
## This function untransforms the values from logit units to probability
p_survival2 <- exp(survival_est) / (1 + exp(survival_est))

# Create a dataframe witplot <- ggplot(result, aes(x = Trt, y = fit))+h the results
survival_result2 <- cbind(nd, p_survival2);survival_result2
survival_result$Treatment <- as.integer(as.character(survival_result$Treatment))

baby_survival_plot <- ggplot(survival_result2, aes(x = Treatment, y = fit))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3, col=NA)+
  geom_point(data = mes, aes(x = Treatment, y = survival_prob),  position=position_jitter(h=0.03,w=0.3))+
  #labs(tag="A")+
  scale_x_continuous("Temperature (C)")+
  scale_y_continuous("Survival Probability", limits = c(0,0.6,0.1))+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
baby_survival_plot

library(DHARMa)
#model assumptions met
testOutliers(survival2)
testDispersion(survival2)



# Temperature on Molting Rate ---------------------------------------------

#BABY LONG-TERM EXPOSURE


#Does temperature affect the molting rate (total molts/days until adulthood/mortality)

babydf <- babydf %>%
  mutate(MoltRate = case_when(!is.na(AdulthoodDay) ~ AdulthoodDay/MoltCount, Survival == 0 | is.na(AdulthoodDay)  ~ Death_Day/MoltCount))


hist(babydf$MoltRate) #Molt rate is right skewed

babydf <- babydf %>%
  mutate(logMoltRate = case_when(MoltRate > 0 ~ log(MoltRate)))


hist(log(babydf$MoltRate)) #log transformed is normally distributed 

moltdf <- subset(babydf, MoltCount > 0) #omitting individuals that died before molting 


babymolt <- lmer(logMoltRate ~ Treatment + (1|Block/Mesocosm), data = moltdf)
summary(babymolt)
drop1(babymolt, test = 'Chisq')

babymolt2 <- lmer(logMoltRate ~ poly(Treatment,2) + (1|Block/Mesocosm), data = moltdf)
summary(babymolt2)
drop1(babymolt2, test = 'Chisq')

anova(babymolt, babymolt2)

library(DHARMa)
#model assumptions met
testOutliers(babymolt2)
testDispersion(babymolt2)

conditional <- partR2(babymolt2, data = babydf, R2_type = "conditional", nboot = 5) #fixed + random effects
marginal <- partR2(babymolt2, data = babydf, R2_type = "marginal", nboot = 5) #fixed effects only

#babydf$logMoltRate[is.na(babydf$logMoltRate)] <- 0
babydf <- babydf[!is.na(babydf$logMoltRate),]

moltc <- babydf %>%
  group_by(Treatment) %>%
  summarise(mean = mean(logMoltRate), sd = sd(logMoltRate), N = n()) %>%
  mutate(se = sd / sqrt(N),
         lwr = mean - qt(1 - (0.05 / 2), N - 1) * se,
         upr = mean + qt(1 - (0.05 / 2), N - 1) * se)

moltc <- moltc %>%
  mutate(meant = exp(mean), sdt = exp(sd), set = exp(se), lwrt = exp(lwr), uprt = exp(upr))



baby_molt_plot <- ggplot(moltc, aes(x = Treatment, y = meant))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwrt, ymax=uprt), alpha=0.5, col=NA)+
  #labs(tag="A")+
  scale_x_continuous("Temperature (C)")+
  scale_y_continuous("Mean Molt Rate")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
baby_molt_plot 




nd <- data.frame(unique(subset(moltdf, select =c(Treatment))), "Mesocosm"=NA, "Block" = NA)

mySumm_mod <- function(.){
  predict(., newdata=nd, re.form=~0)
}

boot <- bootMer(babymolt, mySumm_mod, nsim = 200, use.u=F, type="parametric")

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

# Create a dataframe witplot <- ggplot(result, aes(x = Trt, y = fit))+h the results
mr <- cbind(nd, molt_est);mr
mr$Treatment <- as.integer(as.character(mr$Treatment))

baby_molt_plot <- ggplot(mr, aes(x = Treatment, y = fit))+
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
baby_molt_plot 


#predicting the molt rate for the entire df - these predictions will be used in the next model

nd <- moltdf

mySumm_mod <- function(.){
  predict(., newdata=nd, re.form=~0)
}

boot <- bootMer(babymolt, mySumm_mod, nsim = 200, use.u=F, type="parametric")

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

# Create a dataframe witplot <- ggplot(result, aes(x = Trt, y = fit))+h the results
mr <- cbind(nd, molt_est);mr
mr$Treatment <- as.integer(as.character(mr$Treatment))

baby_molt_plot <- ggplot(mr, aes(x = Treatment, y = fit))+
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
baby_molt_plot 


#use the fitted probability for molt rate to predict survival
molt_survival <- glmer(Survival ~ fit + (1|Block/Mesocosm), family = binomial(link = 'logit'), data = mr)
summary(molt_survival)
drop1(molt_survival, test = "Chisq")

nd <- data.frame(unique(subset(mr, select =c(fit))), "Mesocosm"=NA, "Block" = NA)

mySumm_mod <- function(.){
  predict(., newdata=nd, re.form=~0)
}

boot <- bootMer(molt_survival, mySumm_mod, nsim = 200, use.u=F, type="parametric")

sumBoot <- function(merBoot) {
  return(
    data.frame(survival_fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               survival_lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               survival_upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

## Use the function to estimate the model fit and 95% CI
survival_est <- sumBoot(boot)

psurvival <- exp(survival_est) / (1 + exp(survival_est))

# Create a dataframe witplot <- ggplot(result, aes(x = Trt, y = fit))+h the results
mr_s <- cbind(nd, psurvival);mr_s
mr_s$Treatment <- as.integer(as.character(mr$Treatment))

molt_survival_plot <- ggplot(mr_s, aes(x = fit, y = survival_fit))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=survival_lwr, ymax=survival_upr), alpha=0.5, col=NA)+
  scale_x_continuous("log(Molt Rate)")+
  scale_y_continuous("Survival Probability")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
molt_survival_plot 


#5TH INSTAR SHORT-TERM EXPOSURE
#timing until becoming an adult rather than a rate using a single molt point

teendf$Adult_Date <- as.Date(teendf$Adult_Date)

difftime(teendf$Adult_Date,teendf$Start_Date, units = days)





# Temperature on Body condition -------------------------------------------

survived$Age <- "Baby"
survived1$Age <- "Teen"
survived <- subset(survived, select = c("Treatment", "BodyCondition1", "Age"))
survived1 <- subset(survived1, select = c("Treatment", "BodyCondition1", "Age"))

all <- rbind(survived, survived1)
all$Treatment <- as.numeric(as.character(all$Treatment))

#combined plot of baby and teen survival
plotc <- ggplot(all, aes(x = Treatment, y = BodyCondition1, color = Age, fill = Age))+
  geom_point(alpha = 0.5)+
  #geom_line(size=2)+
  geom_smooth()+
  #geom_ribbon(aes(ymin=lwr, ymax=upr, fill = Age), alpha=0.5, col=NA)+
  scale_fill_manual(values = c("skyblue", "coral"), aesthetics = c("color", "fill"))+
  labs(x = "Temperature (C)", y = "Body Condition (g/mm)")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))
plotc


#BABY

survived <- subset(babydf, AdulthoodDay > 0) #omitting NA values where the backswimmer never became an aduwlt
head(survived)

t <- as.data.frame(table(survived$Treatment))
t <- c("Treatment", "Adult Count")
ggplot(data = t, aes(Var1, Freq))+
  geom_bar(stat='identity')+
  labs(x = "Treatment", y = "Number of Adults - Baby Exp")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))


#Body Condition immediately after adult emergence
survived$BodyCondition1 <- survived$Mass1/survived$Thorax_Width1

hist(survived$BodyCondition1)

condition <- lmer(BodyCondition1 ~ poly(Treatment,2) + (1|Block/Mesocosm), data = survived)
summary(condition)
drop1(condition, test=  'Chisq')

anova(c, condition)

conditional <- partR2(babymolt2, data = babydf, R2_type = "conditional", nboot = 5) #fixed + random effects
marginal <- partR2(babymolt2, data = babydf, R2_type = "marginal", nboot = 5) #fixed effects only


library(DHARMa)
#model assumptions met
testOutliers(condition)
testDispersion(condition)

bc <- survived %>%
  group_by(Treatment) %>%
  summarise(BCmean = mean(BodyCondition1), BCsd = sd(BodyCondition1), N = n()) %>%
  mutate(seBC = BCsd / sqrt(N),
         lwr = BCmean - qt(1 - (0.05 / 2), N - 1) * seBC,
         upr = BCmean + qt(1 - (0.05 / 2), N - 1) * seBC)
survived$Timing <- "After emergence"


baby_bc <- ggplot(bc, aes(x = Treatment, y = BCmean))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.5, col=NA)+
  #labs(tag="A")+
   scale_x_continuous("Temperature (C)")+
  scale_y_continuous("Mean body condition (g/mm)")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
baby_bc


#Thorax Width immediately after adult emergence
hist(survived$Thorax_Width1)
survived$Treatment <- as.factor(survived$Treatment)

plot(survived$Treatment, survived$Mass1)

thorax.baby <- lmer(Thorax_Width1 ~ poly(Treatment,2) + (1|Block/Mesocosm), data = survived)
summary(thorax.baby)
drop1(thorax.baby, test=  'Chisq')

#Mass immediately after adult emergence
hist(survived$Mass1)
survived$Treatment <- as.factor(survived$Treatment)

mass.baby <- lmer(Mass1 ~ poly(Treatment,2) + (1|Block/Mesocosm), data = survived)
summary(mass.baby)
drop1(mass.baby, test=  'Chisq')


#Now for individuals that made it to the very end of the exp 
s2 <- subset(survived, Thorax._Width2 >0) #selecting individuals that made it to second measurement

#body condition

s2$BodyCondition2 <- s2$Mass2 / s2$Thorax._Width2
condition2 <- lmer(BodyCondition2 ~ poly(Treatment,2) + (1|Block/Mesocosm), data = s2)
summary(condition2)
drop1(condition2, test=  'Chisq')

library(DHARMa)
#model assumptions met
testOutliers(condition2)
testDispersion(condition2)

bc2 <- s2 %>%
  group_by(Treatment) %>%
  summarise(BCmean = mean(BodyCondition2), BCsd = sd(BodyCondition2), N = n()) %>%
  mutate(seBC = BCsd / sqrt(N),
         lwr = BCmean - qt(1 - (0.05 / 2), N - 1) * seBC,
         upr = BCmean + qt(1 - (0.05 / 2), N - 1) * seBC)
bc2$Exposure <- "Long-term"

s2$Timing <- "Before dispersal assay"

babyafter <- ggplot(bc2, aes(x = Treatment, y = BCmean))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.5, col=NA)+
  #labs(tag="A")+
  scale_x_continuous("Temperature (C)")+
  scale_y_continuous("Mean body condition (g/mm)")+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
babyafter

bc$Timing <- "After emergence"
bc2$Timing <- "Before dispersal assay"
bc_combined <- rbind(bc, bc2)


combined_bc <- ggplot(bc_combined, aes(x = Treatment, y = BCmean, col = Timing, fill = Timing))+
  geom_line(size = 2)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.25, col=NA)+
  geom_point(data = s2, aes(x = Treatment, y = BodyCondition2),  position=position_jitter(h=0.001,w=0.3), col = "coral", size = 2)+
  geom_point(data = survived, aes(x = Treatment, y = BodyCondition1),  position=position_jitter(h=0.001,w=0.3), col = "skyblue", size = 2)+
 # geom_point(data = s2, aes(x = Treatment, y = BodyCondition2), col = 'coral', position=position_jitter(h=0.001,w=0.2),size = 2)+
  scale_fill_manual(values = c("After emergence" = "skyblue","Before dispersal assay" = "coral"), aesthetics = c("color", "fill"))+
labs(tag="A")+
  scale_x_continuous("Temperature (C)")+
  scale_y_continuous("Mean body condition (g/mm)", limits = c(0.014,0.024, 0.02))+
  theme_bw(base_size=12)+ 
  theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
  theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
  theme(plot.tag = element_text(size = rel(2)))
combined_bc


#how did the composition of individuals change between before and after bc measurements
stable <- table(survived$Treatment)
s2table <- table(s2$Treatment)


#Mass
mass.baby2 <- lmer(Mass2 ~ Treatment + (1|Block/Mesocosm), data = s2)
summary(mass.baby2)
drop1(mass.baby2, test=  'Chisq')


#Thorax width
thorax.baby2 <- lmer(Thorax._Width2 ~ Treatment + (1|Block/Mesocosm), data = s2)
summary(thorax.baby2)
drop1(thorax.baby2, test=  'Chisq')


#Cor between thorax and mass
cor(survived$Thorax_Width1, survived$Mass1, method = 'pearson')
plot(survived$Thorax_Width1, survived$Mass1)




#TEEN
  survived1 <- subset(teendf, Ambient_Days > 0) #omitting NA values where the backswimmer never became an aduwlt
  head(survived)
  survived1$Treatment <- as.factor(survived$Treatment)
  survived1$BodyCondition1 <- survived1$Mass1/survived1$Thorax_Width1
  
  
  
  t <- as.data.frame(table(survived$Treatment))
  ggplot(data = t, aes(Var1, Freq))+
    geom_bar(stat='identity')+
    labs(x = "Treatment", y = "Number of Adults - Teen Exp")+
    theme_bw(base_size=12)+ 
    theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
    theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
          axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))
  
  hist(survived1$BodyCondition1)
  
  teenbc <- lmer(BodyCondition1 ~ Treatment + (1|Block/Mesocosm), data = survived1)
  summary(condition)
  drop1(teenbc, test=  'Chisq')
  teenbc2 <- lmer(BodyCondition1 ~ poly(Treatment,2) + (1|Block/Mesocosm), data = survived1)
  summary(condition)
  drop1(teenbc2, test=  'Chisq')
  
  library(DHARMa)
  #model assumptions met
  testOutliers(teenbc2)
  testDispersion(teenbc2)
  
  survived1$Timing <- "After emergence"
  
  anova(teenbc, teenbc2)
  
#Thorax width after emergence
  teen_thorax <- lmer(Thorax_Width1 ~ Treatment + (1|Block/Mesocosm), data = survived1)
  summary(teen_thorax)
  drop1(teen_thorax, test=  'Chisq')
  
  
#Mass after emergence
  teen_mass <- lmer(Mass1 ~ Treatment + (1|Block/Mesocosm), data = survived1)
  summary(teen_mass)
  drop1(teen_mass, test=  'Chisq')
  
  
  #Teen - before dispersal assay
  teen2 <- subset(survived1, Thorax_Width2 > 0) #omitting adults that died before dispersal assess day
  
  #Body condition after emergence
  teen2$BodyCondition2 <- teen2$Mass2 / teen2$Thorax_Width2
  teen_condition2 <- lmer(BodyCondition2 ~ Treatment + (1|Block/Mesocosm), data = teen2)
  summary(teen_condition2)
  drop1(teen_condition2, test=  'Chisq')
  
  teenbc2 <- lmer(BodyCondition2 ~ poly(Treatment,2) + (1|Block/Mesocosm), data = teen2)
  drop1(teenbc2, test=  'Chisq')
  
  library(DHARMa)
  #model assumptions met
  testOutliers(teenbc2)
  testDispersion(teenbc2)
  
  anova(teen_condition2, teenbc2)
  
  teen2$Timing <- "Before dispersal assay"
  
  
  t1 <- survived1 %>%
    group_by(Treatment) %>%
    summarise(BCmean = mean(BodyCondition1), BCsd = sd(BodyCondition1), N = n()) %>%
    mutate(seBC = BCsd / sqrt(N),
           lwr = BCmean - qt(1 - (0.05 / 2), N - 1) * seBC,
           upr = BCmean + qt(1 - (0.05 / 2), N - 1) * seBC)
  t1$Timing <- "After emergence"
  
  t2 <- teen2 %>%
    group_by(Treatment) %>%
    summarise(BCmean = mean(BodyCondition2), BCsd = sd(BodyCondition2), N = n()) %>%
    mutate(seBC = BCsd / sqrt(N),
           lwr = BCmean - qt(1 - (0.05 / 2), N - 1) * seBC,
           upr = BCmean + qt(1 - (0.05 / 2), N - 1) * seBC)
  
  t2$Timing <- "Before dispersal assay"
  
  teenbc_combined <- rbind(t1,t2)
  
  
  
  teencombined_bc <- ggplot(teenbc_combined, aes(x = Treatment, y = BCmean, col = Timing, fill = Timing))+
    geom_line(size = 2)+
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.25, col=NA)+
    geom_point(data = teen2, aes(x = Treatment, y = BodyCondition2),  position=position_jitter(h=0.001,w=0.3), col = "coral", size = 2)+
    geom_point(data = survived, aes(x = Treatment, y = BodyCondition1),  position=position_jitter(h=0.001,w=0.3), col = "skyblue", size = 2)+
    # geom_point(data = s2, aes(x = Treatment, y = BodyCondition2), col = 'coral', position=position_jitter(h=0.001,w=0.2),size = 2)+
    scale_fill_manual(values = c("After emergence" = "skyblue","Before dispersal assay" = "coral"), aesthetics = c("color", "fill"))+
    labs(tag="B")+
    scale_x_continuous("Temperature (C)")+
    scale_y_continuous("Mean body condition (g/mm)", limits = c(0.014,0.024, 0.02))+
    theme_bw(base_size=12)+ 
    theme(panel.border = element_rect(color="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))+
    theme(axis.text.x = element_text(size=rel(1.5)), axis.text.y = element_text(size=rel(1.5)),
          axis.title.x = element_text(size=rel(2)),axis.title.y = element_text(size=rel(2)))+
    theme(plot.tag = element_text(size = rel(2)))
  teencombined_bc
  
  
  #how did the composition of individuals change between before and after bc measurements
  ttable <- table(survived1$Treatment)
  t2table <- table(teen2$Treatment)
  
  
  ggarrange(combined_bc, teencombined_bc, ncol = 2, legend =  "right", common.legend = TRUE)
  
  
  #Thorax width after emergence
  teen_thorax2 <- lmer(Thorax_Width2 ~ Treatment + (1|Block/Mesocosm), data = teen2)
  summary(teen_thorax2)
  drop1(teen_thorax2, test=  'Chisq')
  
  
  #Mass after emergence
  teen_mass2 <- lmer(Mass2 ~ Treatment + (1|Block/Mesocosm), data = teen2)
  summary(teen_mass2)
  drop1(teen_mass2, test=  'Chisq')
  
  

# Temperature on thorax width ---------------------------------------------

  

