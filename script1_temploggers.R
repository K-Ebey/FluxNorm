#Script 1: R code to analyze incubator temperature loggers and generate Figure 1A and
#information in Table 3.
library(tidyverse)
library(lubridate)
library(hms)
#read in data
all_temps <- read_csv("incubator_temperature_loggers.csv")
#filter data to only include times when the logger was in the incubator
#30Flux filter
all_temps <-all_temps[!(all_temps$Treatment=='30FLUX')|(!((all_temps$Treatment=='30FLUX')&(all_temps$DateTime<mdy("6/27/2024")))&!((all_temps$Treatment=='30FLUX')&(all_temps$DateTime>mdy("8/20/2024")))),]
#31Flux filter
all_temps <-all_temps[!(all_temps$Treatment=='31FLUX')|(!((all_temps$Treatment=='31FLUX')&(all_temps$DateTime<mdy("6/27/2024")))&!((all_temps$Treatment=='31FLUX')&(all_temps$DateTime>mdy("8/20/2024")))),]
#32Flux filter
all_temps <-all_temps[!(all_temps$Treatment=='32FLUX')|(!((all_temps$Treatment=='32FLUX')&(all_temps$DateTime<mdy("6/27/2024")))&!((all_temps$Treatment=='32FLUX')&(all_temps$DateTime>mdy("8/17/2024")))),]
#33Flux filter
all_temps <-all_temps[!(all_temps$Treatment=='33FLUX')|(!((all_temps$Treatment=='33FLUX')&(all_temps$DateTime<mdy("6/27/2024")))&!((all_temps$Treatment=='33FLUX')&(all_temps$DateTime>mdy("8/11/2024")))),]
#34Flux filter
all_temps <-all_temps[!(all_temps$Treatment=='34FLUX')|(!((all_temps$Treatment=='34FLUX')&(all_temps$DateTime<mdy("6/27/2024")))&!((all_temps$Treatment=='34FLUX')&(all_temps$DateTime>mdy("8/15/2024")))),]
all_temps$Time <- as_hms(all_temps$DateTime)
all_temps$hour <- hour(all_temps$Time)
#Make a figure of average temperatures +/- standard deviation for each temperature treatment and time
all_temps %>%
  group_by(Treatment, hour) %>%
  summarize(meant=mean(Temp),sdt=sd(Temp),ci_low=meant-sdt,ci_high=meant+sdt)%>%
  ggplot(aes(x=hour, y=meant))+geom_line(aes(color=Treatment))+theme_classic()+
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill=Treatment),
              alpha=0.2)+ylab("Temperature °C")+xlab("Hour")
#Table data
all_temps$day <- date(all_temps$DateTime)
all_temps %>% group_by(Treatment, day) %>%
  #calculate the mean, min, max temps each day in each treatment group
  summarize(meant=mean(Temp),
            min=min(Temp),
            max=max(Temp),
            pc95=quantile(Temp,0.95),
            pc5=quantile(Temp,0.05)) %>%
  #calculate the average daily mean,min,max temps in each treatment group
  group_by(Treatment) %>%
  summarize(meand=mean(meant),maxd=mean(max),mind=mean(pc5))
#Calculate average CTE
#using method described in Georges 2004, assuming each day is roughly sinusoidal
#calculate the CTE every day and then average them
#Georges defines CTE as the temp where median development occurs
#calculation method is to find the time (x) when this median temperature is reached and then solve for the temperature
#0=pi/2-(Rsin(x)/(M-T0))-x
#M is the mean temperature and R is the amplitude
solvethis <- function(M, R, x){
  answer=(3.1415926358979323846264/2)-(R*sin(x)/(M-28))-x
}
#x has a range 0 to 2pi, solve for x by cycling through
calc_CTE <- function(M, R){
  answer=1
  x=0
  while(answer>0.01){
    answer=solvethis(M,R,x)
    x=x+0.01
  }
  CTE=R*cos(x)+M
}
#I am using 5th percentile in place of min for emect of door opening
all_temps %>%
  group_by(Treatment, day) %>%
  summarise(M=mean(Temp), R=max(Temp)-quantile(Temp,0.05)) %>%
  group_by(Treatment, day) %>%
  summarise(CTE=calc_CTE(M, R/2)) %>%
  group_by(Treatment) %>%
  summarise(mean_CTE=mean(CTE))
#calculate nominal CTE for experimental treatments
print(calc_CTE(30, 1.5))
print(calc_CTE(31, 1.5))
print(calc_CTE(32, 1.5))
print(calc_CTE(33, 1.5))
print(calc_CTE(34, 1.5))
