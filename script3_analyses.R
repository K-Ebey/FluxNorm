#Script 3: R code to analyze data on egg viability, assign sex, generate reaction norms,
#compare morphology to gene expression and generate all other figures and table values.

library(tidyverse)
library(lmerTest)
library(lme4)
library(glmmTMB)
library(car)
library(MASS)
library(embryogrowth)

load("gene_expression_data.RData") #Gene expression data with cleaning 
load("gene_expression_no_omit_data.RData") #Gene expression data without cleaning

eggs <- read.csv('2024FN_COMBINED.csv')

#for embryos scored with Y assign oviduct score of 3
eggs[eggs$Oviduct.=='Y', 'Oviduct.'] <- 3

#Table 1 sample numbers

#filter out staged/discarded eggs not included in the sample
eggs1 <- eggs %>%
  filter(!(Comments %in% c("Staged","DISCARDED",
                           "Staged 20240625 (stage 13.0)",
                           "SHELL CRUSHED ON SIDE/DISCARDED"))) %>%
  dplyr::select(Egg.ID, Clutch.ID, Treatment, Mass..g., 
                Dissection.Egg.Mass..g., Yolk.Mass..g., Embryo.Mass..g.,
                Oviduct., Survive, Dissection.Tube.ID)

#make survival binary values
eggs1[eggs1$Survive!='D', 'Survive'] <- 1
eggs1[eggs1$Survive=='D', 'Survive'] <- 0
eggs1$Survive <- as.numeric(eggs1$Survive)

#Make table with all eggs incubated by clutch and treatment
eggs1 %>% group_by(Clutch.ID,Treatment) %>% filter(!Treatment=='33CON') %>%
  filter(Dissection.Tube.ID!='FN82') %>% #one of the eggs had twins
  #this table is by number of eggs, so removing one of the twins
  summarise(count = n()) %>%
  pivot_wider(names_from=Treatment,
              values_from = count)

#Make table with only embryos sampled by clutch and treatment
eggs1 %>% filter(Survive==1 & grepl("FN", Dissection.Tube.ID,fixed=T)) %>%
  #filter only samples that survived and were sampled
  group_by(Clutch.ID,Treatment) %>%
  #filter by clutch and treatment
  summarise(count = n()) %>% #get number of samples in the groups
  pivot_wider(names_from=Treatment,
              values_from = count)

#Make table with nonviable eggs by clutch and treatment
eggs1 %>% filter(Survive==0) %>%
  #filter only samples did not survive
  group_by(Clutch.ID,Treatment) %>%
  #filter by clutch and treatment
  summarise(count = n()) %>% #get number of samples in the groups
  pivot_wider(names_from=Treatment,
              values_from = count)


#graph survival by temperature treatment and clutch ID
#Figure 1====
eggs1 %>% group_by(Treatment) %>% filter(Treatment!='33CON') %>%
  filter(Dissection.Tube.ID!='FN82') %>% #one of the eggs had twins
  #this analysis is by number of eggs, so removing one of the twins
  summarize(survival_rate=mean(Survive)) %>%
  ggplot(aes(x=Treatment, y=survival_rate))+geom_bar(stat='identity')+
  theme_classic()+ylab('Proportion Viable') + ylim(0,1)

summary(glmer(Survive~Treatment+(1|Clutch.ID),
              data=eggs1[!(eggs1$Treatment=='33CON')&eggs1$Dissection.Tube.ID!='FN82',], family=binomial(link="logit")))
#one of the eggs had twins; this analysis is by number of eggs, so removing one of the twins

Anova(glmer(Survive~Treatment+(1|Clutch.ID),
            data=eggs1[!(eggs1$Treatment=='33CON')&eggs1$Dissection.Tube.ID!='FN82',], family=binomial(link="logit")))
#one of the eggs had twins; this analysis is by number of eggs, so removing one of the twins

#make dataframe with egg data and gene expression data
comb <- merge(eggs, all_genes, by='Dissection.Tube.ID') #cleaned gene expression data
comb_no_omit <-merge(eggs, all_no_omit, by='Dissection.Tube.ID') #uncleaned gene expression data

comb_33con <- comb %>% filter(Treatment=='33CON') #33 constant treatment

#only include fluctuating temperature treatments
comb <- comb %>% filter(!(Treatment=='33CON'))
comb_no_omit <- comb_no_omit %>% filter(!(Treatment=='33CON'))

#figure A2
#visualize how gene expression overall-RPL8 and standardized-Cyp19 and AMH varies by treatment group
ggplot(comb, aes(x=Treatment, y=log(rpl8_SQ)))+geom_boxplot()+theme_classic()

#test if rpl8 levels depend on temperature treatment
shapiro.test(log(comb[!(comb$Treatment=='33CON'),]$rpl8_SQ))
#rpl8 expression is not normal p<0.001, so use Kruskal-Wallis test
kruskal.test(rpl8_SQ~Treatment, data=comb[!(comb$Treatment=='33CON'),])

#Assign sex====
comb$Sex <- NA
comb_no_omit$Sex <- NA
#1=F
#0=M
#sex is assigned based on log(St_cyp19)>-6 based on two clumps in the graph

#Assign sex for flux treatments with cleaned gene expression data
comb[!is.na(comb$st_cyp19) & is.na(comb$st_amh)& log(comb$st_cyp19)> -6,]$Sex <- 1
comb[is.na(comb$st_cyp19) & !is.na(comb$st_amh),]$Sex <- 0
comb[!is.na(comb$st_cyp19) & !is.na(comb$st_amh) & log(comb$st_cyp19)> -6,]$Sex <- 1
comb[!is.na(comb$st_cyp19) & !is.na(comb$st_amh) & log(comb$st_cyp19)< -6,]$Sex <- 0


#Assign sex for 33 constant treatment with cleaned gene expression data
comb_33con$Sex <- NA
#sex is assigned based on log(St_cyp19)>-6 based on two clumps in the graph
comb_33con[!is.na(comb_33con$st_cyp19) & is.na(comb_33con$st_amh)& log(comb_33con$st_cyp19)> -6,]$Sex <- 1
comb_33con[is.na(comb_33con$st_cyp19) & !is.na(comb_33con$st_amh),]$Sex <- 0
comb_33con[!is.na(comb_33con$st_cyp19) & !is.na(comb_33con$st_amh) & log(comb_33con$st_cyp19)> -6,]$Sex <- 1
comb_33con[!is.na(comb_33con$st_cyp19) & !is.na(comb_33con$st_amh) & log(comb_33con$st_cyp19)< -6,]$Sex <- 0

#Assign sex for flux treatments with uncleaned gene expression data
comb_no_omit[!is.na(comb_no_omit$st_cyp19) & is.na(comb_no_omit$st_amh)& log(comb_no_omit$st_cyp19)> -6,]$Sex <- 1
comb_no_omit[is.na(comb_no_omit$st_cyp19) & !is.na(comb_no_omit$st_amh),]$Sex <- 0
comb_no_omit[!is.na(comb_no_omit$st_cyp19) & !is.na(comb_no_omit$st_amh) & log(comb_no_omit$st_cyp19)> -6,]$Sex <- 1
comb_no_omit[!is.na(comb_no_omit$st_cyp19) & !is.na(comb_no_omit$st_amh) & log(comb_no_omit$st_cyp19)< -6,]$Sex <- 0

#check that sex is the same for cleaned and uncleaned gene expression data
comb_no_omit$Sex==comb$Sex

#get numbers of males and females by clutch and treatment
comb %>%
  group_by(Clutch.ID,Treatment) %>% filter(Sex==0) %>% summarise(count = n()) %>% 
  pivot_wider(names_from=Treatment, values_from = count)
comb %>%
  group_by(Clutch.ID,Treatment) %>% filter(Sex==1) %>% summarise(count = n()) %>% 
  pivot_wider(names_from=Treatment, values_from = count)

comb_33con %>%
  group_by(Clutch.ID) %>% filter(Sex==0) %>%summarise(count = n())
comb_33con %>%
  group_by(Clutch.ID) %>% filter(Sex==1) %>%
  summarise(count = n())

#figure 2 a b and c====

#for visualization purposes assign Na values
comb2 <- comb 
comb2$st_cyp19 <- replace_na(comb$st_cyp19,0.000001)
comb2$st_amh <- replace_na(comb$st_amh,0.0000001)

ggplot(comb2, aes(x=log(st_amh), y=log(st_cyp19)))+
  geom_point(aes(color=as.factor(Sex)))+theme_classic()+geom_hline(yintercept = -6)+
  guides(fill="none")+guides(color="none")+scale_fill_manual(values=c('darkgreen', 'orange'))+scale_color_manual(values=c('darkgreen', 'orange'))

ggplot(comb[!is.na(comb$Sex),], aes(x=Treatment, y=log(st_amh), fill=as.factor(Sex)))+
  geom_boxplot(alpha=0.6)+theme_classic()+geom_jitter(width=0.15, aes(color=as.factor(Sex)))+ (fill="none")+guides(color="none")+ scale_fill_manual(values=c('darkgreen', 'orange'))+ scale_color_manual(values=c('darkgreen', 'orange'))

ggplot(comb[!is.na(comb$Sex) & !(comb$Treatment=='33CON'),], aes(x=Treatment, y=log(st_cyp19), fill=as.factor(Sex)))+
  geom_boxplot(alpha=0.6)+theme_classic()+geom_jitter(width=0.15, aes(color=as.factor(Sex)))+ (fill="none")+guides(color="none")+ scale_fill_manual(values=c('darkgreen', 'orange'))+scale_color_manual(values=c('darkgreen', 'orange'))

#Test if cyp19 and amh expression varies by treatment
#replace nondetect amh and cyp19 expression with 0
comb2$st_amh <- replace_na(comb$st_amh,0)
comb2$st_cyp19 <- replace_na(comb$st_cyp19,0)

#cyp19 and amh distriubtions are non-normal
shapiro.test(comb2[!is.na(comb2$Sex)&comb2$Sex==0,]$st_amh)
shapiro.test(comb2[!is.na(comb2$Sex)&comb2$Sex==1,]$st_amh)
shapiro.test(comb2[!is.na(comb2$Sex)&comb2$Sex==0,]$st_cyp19)
shapiro.test(comb2[!is.na(comb2$Sex)&comb2$Sex==1,]$st_cyp19)

kruskal.test(st_cyp19~Treatment, data=comb2[!is.na(comb2$Sex)&comb2$Sex==1,])
kruskal.test(st_cyp19~Treatment, data=comb2[!is.na(comb2$Sex)&comb2$Sex==0,])
pairwise.wilcox.test(comb2[!is.na(comb2$Sex)&comb2$Sex==1,]$st_cyp19,comb2[!is.na(comb2$Sex)&comb2$Sex==1,]$Treatment,p.adjust.method = 'bonferroni')
pairwise.wilcox.test(comb2[!is.na(comb2$Sex)&comb2$Sex==0,]$st_cyp19,comb2[!is.na(comb2$Sex)&comb2$Sex==0,]$Treatment,p.adjust.method = 'bonferroni')

kruskal.test(st_amh~Treatment, data=comb2[!is.na(comb2$Sex)&comb2$Sex==0,])
kruskal.test(st_amh~Treatment, data=comb2[!is.na(comb2$Sex)&comb2$Sex==1,])



#treat temperature as continuous variable 
#Temperature values are empirical means from temps.R results
comb[comb$Treatment=='30FLUX', 'Temperature'] <- 30
comb[comb$Treatment=='31FLUX', 'Temperature'] <- 31
comb[comb$Treatment=='32FLUX', 'Temperature'] <- 32
comb[comb$Treatment=='33FLUX', 'Temperature'] <- 33.1
comb[comb$Treatment=='34FLUX', 'Temperature'] <- 33.9

#CTEs for each treatment from temps.R results
comb[comb$Treatment=='30FLUX', 'CTE'] <- 30.9
comb[comb$Treatment=='31FLUX', 'CTE'] <- 31.7
comb[comb$Treatment=='32FLUX', 'CTE'] <- 32.5
comb[comb$Treatment=='33FLUX', 'CTE'] <- 33.6
comb[comb$Treatment=='34FLUX', 'CTE'] <- 34.4

#Morphological Sex====

comb$Morph_Sex <- NA
#assign morphological sex based on oviduct score 0 and 1 are male; 2 and 3 are female
comb[comb$Oviduct.%in%c(0,1), 'Morph_Sex'] <- 0
comb[comb$Oviduct.%in%c(2,3), 'Morph_Sex'] <- 1

#figure 4 a====
eggs1 %>% filter(Treatment!='33CON') %>% filter(Survive==1) %>% filter(Oviduct.%in%c(0,1,2,3)) %>%ggplot(aes(x=Treatment, fill=as.factor(Oviduct.), y=Mass..g.))+ geom_bar(position='fill', stat='identity')+ theme_classic()+ylab('Proportion')


#compare morphological sex results to gene expression sex
nrow(comb[comb$Treatment!='33CON'&!is.na(comb$Morph_Sex)&!is.na(comb$Sex),])
View(comb %>% filter(Treatment!='33CON') %>%filter(!is.na(Morph_Sex)) %>%
       filter(!is.na(Sex)) %>% filter(Morph_Sex!=Sex) %>% 
       dplyr::select(c(Dissection.Tube.ID, Treatment, Oviduct., st_cyp19,st_amh, Morph_Sex, Sex)) %>%mutate(log_st_cyp19=log(st_cyp19)))

#figure 4 b ====
ggplot(comb[comb$Oviduct.%in%c('0','1','2','3')&!comb$Treatment=='33CON',], aes(y=log(st_amh), x=Oviduct.))+
  geom_boxplot()+geom_jitter(width=0.2, alpha=0.9, height=0, aes(color=as.factor(Sex)))+
  geom_smooth()+theme_classic()+scale_color_manual(values=c('darkgreen', 'orange'))+
  xlab("Oviduct score")
       
       kruskal.test(st_amh~as.factor(Oviduct.),data=comb2[comb2$Oviduct.%in%c('0','1','2','3')&!comb2$Treatment=='33CON'&!is.na(comb$Sex),])
       pairwise.wilcox.test(comb2[comb2$Oviduct.%in%c('0','1','2','3')&!comb2$Treatment=='33CON'&!is.na(comb$Sex),]$st_amh,                 as.factor(comb2[comb2$Oviduct.%in%c('0','1','2','3')&!comb2$Treatment=='33CON'&!is.na(comb$Sex),]$Oviduct.),p.adjust.method = 'bonferroni')             
       
       #embryogrowth package reaction norm====
       constant_temp <- read_csv("Constant_temp_RxnNorm_data.csv")
       cons_temp <- na.omit(constant_temp %>%dplyr::select(`Temperatures`,males,females))
       
       
       cons_temp <- cons_temp %>% rename(temperatures=`Temperatures`)
       
       
       ct <-pivot_longer(cons_temp,cols=c(males,females),names_to='Sex')
       
       Temperature <- c()
       Sex <- c()
       
       for (n in (1:nrow(ct))){
         Tmp <- as.numeric(ct[n,'temperatures'])
         Temperature <- c(Temperature, rep(Tmp,ct[n,'value']))
         if (ct[n,'Sex']=='males'){
           Sex <- c(Sex, rep(0,ct[n,'value']))
         }
         if (ct[n,'Sex']=='females'){
           Sex <- c(Sex, rep(1,ct[n,'value']))
         }
       }
       
       cst <- data.frame(Temperature=Temperature,Sex=Sex)
       
       
       #constatn temperature reaction norm
       model_TSDII <- tsd(males=cons_temp$males,females=cons_temp$females,temperatures=cons_temp$temperatures, males.freq=FALSE, parameters.initial=c(P_low=31.5, S_low=0.3, P_high=34, S_high=-0.4),equation='logistic')
       
       priors <- tsd_MHmcmc_p(result=model_TSDII, accept=TRUE)
       out_mcmc <- tsd_MHmcmc(result=model_TSDII, n.iter=10000, parametersMCMC=priors)
       plot(model_TSDII, resultmcmc=out_mcmc, lab.TRT = "TRT l = 5 %",males.freq = F)
       
       #print upper and lower pivotal temperature credibility interval
       data.frame(summary(out_mcmc)[2])[c(1,3),]
       
       #Girondot 1999 defines TRT as S*2*ln(l/(1-l))
       #Define TRT as 5-95% sex ratio, so l=0.05
       #print |upper and lower TRT|
       data.frame(summary(out_mcmc)[2])[c(2,4),]*2*log(1/19)
       
       #Figure 3A
       #plot reaction norm with input data
       plot(model_TSDII, resultmcmc=out_mcmc, lab.TRT = "TRT l = 5 %",males.freq = F,
            show.observations = F,show.PTRT = F, ylab='Proportion female', xlab='Temperature (°C)')+
         geom_jitter(data=cst, aes(x=Temperature, y=Sex),width = 0.2, height = 0.05,alpha=0.2)
       
       #Set up Flux dataframe by copying values from comb
       comb %>%
         group_by(Treatment) %>% filter(Sex==0) %>%
         summarise(count = n())
       
       comb %>%
         group_by(Treatment) %>% filter(Sex==1) %>%
         summarise(count = n())
       unique(comb$Temperature)
       unique(comb$CTE)
       
       flux_avgtemp <- data.frame(temperatures=c(30,31,32,33.1,33.9),
                                  males=c(0,0,23,22,4),females=c(25,30,3,2,23))
       flux_cte <- data.frame(temperatures=c(30.9,31.7,32.5,33.6,34.4),
                              males=c(0,0,23,22,4),females=c(25,30,3,2,23))
       
       #fluctuating temperature reaction norm based on average temperatures
       model_TSDII_at <- tsd(males=flux_avgtemp$males,females=flux_avgtemp$females,temperatures=flux_avgtemp$temperatures, males.freq=FALSE,parameters.initial=c(P_low=31.5, S_low=0.3, P_high=34, S_high=-0.4),equation='logistic')
       priors_at <- tsd_MHmcmc_p(result=model_TSDII_at, accept=TRUE)
       out_mcmc_at <- tsd_MHmcmc(result=model_TSDII_at, n.iter=10000, parametersMCMC=priors_at)
       plot(model_TSDII_at, resultmcmc=out_mcmc_at, lab.TRT = "TRT l = 5 %",males.freq = F)
       summary(out_mcmc_at)[2]
       
       #Figure 3B
       plot(model_TSDII_at, resultmcmc=out_mcmc_at, lab.TRT = "TRT l = 5 %",males.freq = F, show.observations = F, show.PTRT = F,ylab='Proportion female', xlab='Average Temperature (°C)',xlim=c(29.5,34.5))+ geom_jitter(data=comb, aes(x=Temperature, y=Sex),width = 0.2, height = 0.05,alpha=0.5)
       
       #print Tpiv
       data.frame(summary(out_mcmc_at)[2])[c(1,3),]
       #print |upper and lower TRT|
       data.frame(summary(out_mcmc_at)[2])[c(2,4),]*2*log(1/19)
       
       #fluctuating temperature reaction norm based on CTE
       model_TSDII_cte <- tsd(males=flux_cte$males,females=flux_cte$females,temperatures= flux_cte$temperatures, males.freq=FALSE, parameters.initial=c(P_low=31.5, S_low=0.3, P_high=34, S_high=-0.4), equation='logistic')
       priors_cte <- tsd_MHmcmc_p(result=model_TSDII_cte, accept=TRUE)
       out_mcmc_cte <- tsd_MHmcmc(result=model_TSDII_cte, n.iter=10000, parametersMCMC=priors_cte)
       
       #figure 3C
       plot(model_TSDII_cte, resultmcmc=out_mcmc_cte, lab.TRT = "TRT l = 5 %",males.freq = F, show.observations = F, show.PTRT = F,ylab='Proportion female', xlab='CTE (°C)',xlim=c(30,34.8))+ geom_jitter(data=comb, aes(x=CTE, y=Sex),width = 0.2, height = 0.05,alpha=0.5)
       
       #print Tpiv
       data.frame(summary(out_mcmc_cte)[2])[c(1,3),]
       #print |upper and lower TRT|
       data.frame(summary(out_mcmc_cte)[2])[c(2,4),]*2*log(1/19)
       
       #make dataframe with pivotal temperatures and CI copied from outputs
       PT_data <- data.frame(Reaction_Norm = c(rep('Constant',3), rep("Flux-Avg Temp",3), rep("Flux-CTE",3)),
                             L_Tpiv = c(31.90,31.88,31.96,
                                        31.5,31.92,31.73,
                                        32.29,32.06,32.43),
                             U_Tpiv=c(34.92,34.69,35.27,
                                      33.56,33.36,33.75,
                                      34.06,33.84,34.24))
       
       #Figure 3D
       ggplot(PT_data, aes(y=Reaction_Norm,x=L_Tpiv))+ geom_point(color='turquoise')+ stat_summary(fun=median,fun.min=min,fun.max=max,color='turquoise')+ geom_point(aes(y=Reaction_Norm,x=U_Tpiv),color='salmon')+
         stat_summary(aes(y=Reaction_Norm,x=U_Tpiv),fun=median,fun.min=min,fun.max=max,color='salmon')+ theme_classic()+ylab("")+xlab("Pivotal Temperature (°C)")
       