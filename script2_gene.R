#Script 2: R code to clean and process qPCR data for each gene.

library(tidyverse)
library(lmerTest)
library(lme4)
library(glmmTMB)
library(car)
#function to clean the data
#if CV is above threshold (currently set at 0.25), examine all data and remove a
#point if dimerence between it and the other 2 points exceeds the sd
clean_data <- function(gene_data) {
  samples <- unique(gene_data$Sample)
  result <- data.frame(Sample=c(), Well=c(), Plate=c(), Cq=c(),
                       Starting.Quantity..SQ.=c(), mean_SQ=c(), SD_SQ=c(),
                       CV_SQ=c(), comments=c())
  for (s in samples){
    temp <- gene_data[gene_data$Sample==s,] #go through sample by sample
    temp<- temp %>% mutate(comments=0)
    temp$CV_SQ <- replace_na(temp$CV_SQ, 0)
    if (unique(temp$CV_SQ)>0.25) { #check if sample CV is above threshold 0.25
      sd <- unique(temp$SD_SQ)
      #take the dimerence between each pair of replicates
      if(length(temp$Sample)==3){
        d1 <- abs(temp$Starting.Quantity..SQ.[1]-temp$Starting.Quantity..SQ.[2])
        d2 <- abs(temp$Starting.Quantity..SQ.[1]-temp$Starting.Quantity..SQ.[3])
        d3 <- abs(temp$Starting.Quantity..SQ.[2]-temp$Starting.Quantity..SQ.[3])
        #
        if (d1>sd & d2>sd & d3>sd){ #case where all three dimerences exceed SD
          temp$comments <- "all spread, none removed"
        }
        else if (d1>sd & d2>sd &d3<sd){ #case where 1-2 and 1-3 exceed SD
          temp$comments <- paste(temp$Well[1], "Spread, Removed")
          temp$Starting.Quantity..SQ.[1]<- NA
          temp$mean_SQ=mean(temp$Starting.Quantity..SQ., na.rm=T)
          temp$SD_SQ=sd(temp$Starting.Quantity..SQ., na.rm=T)
          temp$CV_SQ=sd(temp$Starting.Quantity..SQ.,
                        na.rm=T)/mean(temp$Starting.Quantity..SQ., na.rm=T)
        }
        else if (d1>sd & d2<sd &d3>sd){ #case where 1-2 and 2-3 exceed SD
          temp$comments <- paste(temp$Well[2], "Spread, Removed")
          temp$Starting.Quantity..SQ.[2]<- NA
          temp$mean_SQ=mean(temp$Starting.Quantity..SQ., na.rm=T)
          temp$SD_SQ=sd(temp$Starting.Quantity..SQ., na.rm=T)
          temp$CV_SQ=sd(temp$Starting.Quantity..SQ.,
                        na.rm=T)/mean(temp$Starting.Quantity..SQ., na.rm=T)
        }
        else if (d1<sd & d2>sd &d3>sd){ #case where 1-3 and 2-3 exceed SD
          temp$comments <- paste(temp$Well[3], "Spread, Removed")
          temp$Starting.Quantity..SQ.[3]<- NA
          temp$mean_SQ=mean(temp$Starting.Quantity..SQ., na.rm=T)
          temp$SD_SQ=sd(temp$Starting.Quantity..SQ., na.rm=T)
          temp$CV_SQ=sd(temp$Starting.Quantity..SQ.,
                        na.rm=T)/mean(temp$Starting.Quantity..SQ., na.rm=T)
        }
        else { #case where one or no dimerences exceed SD
          temp$comments <- "none spread, none removed"
        }
      }
    }
    temp <- temp %>% select(Sample, Well, Plate, Cq, Starting.Quantity..SQ., mean_SQ,
                            SD_SQ, CV_SQ, comments)
    result <- rbind(result, temp)
  }
  result
}
#RPL8====
rpl8_raw<-read.csv('RPL8_All_Plates.csv')
rpl8_raw <- rpl8_raw %>% filter(!(is.na(Sample))) #remove standard and ntc data
#group by sample because 3 technical replicates
rpl8<-rpl8_raw %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(Starting.Quantity..SQ.), #mean starting quantity
  SD_SQ=sd(Starting.Quantity..SQ.), #standard deviation starting quantity
  CV_SQ=sd(Starting.Quantity..SQ.)/mean(Starting.Quantity..SQ.), #coemicient of variation starting quantity
  mean_cQ=mean(Cq), #mean cycle threshold
  SD_cQ=sd(Cq), #stardard deviation cycle threshold
  CV_cQ=sd(Cq)/mean(Cq)) #coemicient of variation cycle threshold
#No NA.omit here because samples with less than 3 reads will not have any points
removed
ggplot(rpl8, aes(x=CV_SQ))+geom_histogram()
ggplot(rpl8, aes(x=log(mean_SQ)))+geom_histogram(bins=100)
rpl8_full <- merge(rpl8_raw, rpl8, by='Sample') %>% select(Sample, Well, Plate,
                                                           Cq, Starting.Quantity..SQ.,
                                                           mean_SQ, SD_SQ, CV_SQ, mean_cQ, SD_cQ, CV_cQ)
rpl8_new <- clean_data(rpl8_full)
rpl8_new_sum<-rpl8_new %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(na.omit(Starting.Quantity..SQ.)), #mean starting quantity
  SD_SQ=sd(na.omit(Starting.Quantity..SQ.)), #standard deviation starting quantity
  CV_SQ=sd(na.omit(Starting.Quantity..SQ.))/mean(na.omit(Starting.Quantity..SQ.)),
  #coemicient of variation starting quantity
  mean_cQ=mean(na.omit(Cq)), #mean cycle threshold
  SD_cQ=sd(na.omit(Cq)), #stardard deviation cycle threshold
  CV_cQ=sd(na.omit(Cq))/mean(na.omit(Cq))) #coemicient of variation cycle threshold
#Na.Omit here so samples with outlier removed or no reads will be included
ggplot(rpl8_new_sum, aes(x=CV_SQ))+geom_histogram()
ggplot(rpl8_new_sum, aes(x=log(mean_SQ)))+geom_histogram(bins=100)
#AMH====
amh_raw<-read.csv('AMH_All_Plates.csv')
#group by sample because 3 technical replicates
amh<-amh_raw %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(Starting.Quantity..SQ.), #mean starting quantity
  SD_SQ=sd(Starting.Quantity..SQ.), #standard deviation starting quantity
  CV_SQ=sd(Starting.Quantity..SQ.)/mean(Starting.Quantity..SQ.), #coemicient of variation starting quantity
  mean_cQ=mean(Cq), #mean cycle threshold
  SD_cQ=sd(Cq), #stardard deviation cycle threshold
  CV_cQ=sd(Cq)/mean(Cq)) #coemicient of variation cycle threshold
ggplot(amh, aes(x=CV_SQ))+geom_histogram()
ggplot(amh, aes(x=log(mean_SQ)))+geom_histogram(bins=100)
amh_full <- merge(amh_raw, amh, by='Sample') %>% select(Sample, Well, Plate,
                                                        Target, Cq, Starting.Quantity..SQ.,
                                                        mean_SQ, SD_SQ, CV_SQ, mean_cQ, SD_cQ, CV_cQ)
amh_new <- clean_data(amh_full)
amh_new_sum<-amh_new %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(na.omit(Starting.Quantity..SQ.)), #mean starting quantity
  SD_SQ=sd(na.omit(Starting.Quantity..SQ.)), #standard deviation starting quantity
  CV_SQ=sd(na.omit(Starting.Quantity..SQ.))/mean(na.omit(Starting.Quantity..SQ.)),
  #coemicient of variation starting quantity
  mean_cQ=mean(na.omit(Cq)), #mean cycle threshold
  SD_cQ=sd(na.omit(Cq)), #stardard deviation cycle threshold
  CV_cQ=sd(na.omit(Cq))/mean(na.omit(Cq))) #coemicient of variation cycle threshold
ggplot(amh_new_sum, aes(x=CV_SQ))+geom_histogram()
ggplot(amh_new_sum, aes(x=log(mean_SQ)))+geom_histogram(bins=100)
#Cyp19====
cyp19_raw<-read.csv('CYP19_All_Plates.csv')
cyp19_raw <- cyp19_raw %>% filter(!(is.na(Sample))) #remove standard and ntc data
#group by sample because 3 technical replicates
cyp19<-cyp19_raw %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(Starting.Quantity..SQ.), #mean starting quantity
  SD_SQ=sd(Starting.Quantity..SQ.), #standard deviation starting quantity
  CV_SQ=sd(Starting.Quantity..SQ.)/mean(Starting.Quantity..SQ.), #coemicient of variation starting quantity
  mean_cQ=mean(Cq), #mean cycle threshold
  SD_cQ=sd(Cq), #stardard deviation cycle threshold
  CV_cQ=sd(Cq)/mean(Cq)) #coemicient of variation cycle threshold
ggplot(cyp19, aes(x=CV_SQ))+geom_histogram()
ggplot(cyp19, aes(x=log(mean_SQ)))+geom_histogram(bins=100)
cyp19_full <- merge(cyp19_raw, cyp19, by='Sample') %>% select(Sample, Well, Plate,
                                                              Cq, Starting.Quantity..SQ.,
                                                              mean_SQ, SD_SQ, CV_SQ, mean_cQ, SD_cQ, CV_cQ)
cyp19_new <- clean_data(cyp19_full)
cyp19_new_sum<-cyp19_new %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(na.omit(Starting.Quantity..SQ.)), #mean starting quantity
  SD_SQ=sd(na.omit(Starting.Quantity..SQ.)), #standard deviation starting quantity
  CV_SQ=sd(na.omit(Starting.Quantity..SQ.))/mean(na.omit(Starting.Quantity..SQ.)),
  #coemicient of variation starting quantity
  mean_cQ=mean(na.omit(Cq)), #mean cycle threshold
  SD_cQ=sd(na.omit(Cq)), #stardard deviation cycle threshold
  CV_cQ=sd(na.omit(Cq))/mean(na.omit(Cq))) #coemicient of variation cycle threshold
ggplot(cyp19_new_sum, aes(x=CV_SQ))+geom_histogram()
ggplot(cyp19_new_sum, aes(x=log(mean_SQ)))+geom_histogram(bins=100)
#all genes====
#with cleaned data recalculate SQ for each sample
amh_new_sum <- amh_new_sum %>% rename(amh_SQ = mean_SQ) %>% select(Sample,
                                                                   amh_SQ)
cyp19_new_sum <- cyp19_new_sum %>% rename(cyp19_SQ = mean_SQ) %>%
  select(Sample, cyp19_SQ)
rpl8_new_sum <- rpl8_new_sum %>% rename(rpl8_SQ = mean_SQ) %>% select(Sample,
                                                                      rpl8_SQ)
#merge data files
all_genes <- merge(amh_new_sum, rpl8_new_sum, by='Sample')
all_genes <- merge(all_genes, cyp19_new_sum, by='Sample')
#standardize AMH and CYP19A1 gene expression by that of housekeeping RPL8 gene
all_genes$st_amh <- all_genes$amh_SQ/all_genes$rpl8_SQ
all_genes$st_cyp19 <- all_genes$cyp19_SQ/all_genes$rpl8_SQ
ggplot(all_genes, aes(x=st_amh, y=st_cyp19))+geom_point()
ggplot(all_genes, aes(x=st_amh))+geom_histogram()
ggplot(all_genes, aes(x=st_cyp19))+geom_histogram()
all_genes<-all_genes%>%mutate(Dissection.Tube.ID=paste("FN", Sample, sep=""))
#write file with these values
save(all_genes, file="gene_expression_data.RData")

#Without cleaning====
#calculate the gene expression values without data cleaning to see if results are robust
#calculate mean with na.omit to replicate what happens to cleaned data
rpl8<-rpl8_raw %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(na.omit(Starting.Quantity..SQ.)))  #mean starting quantity
amh<-amh_raw %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(na.omit(Starting.Quantity..SQ.)))  #mean starting quantity
cyp19<-cyp19_raw %>% group_by(Sample) %>% summarise(
  mean_SQ=mean(na.omit(Starting.Quantity..SQ.)))  #mean starting quantity



#with cleaned data recalculate SQ for each sample
amh_no_omit <- amh %>% rename(amh_SQ = mean_SQ) %>% dplyr::select(Sample, amh_SQ)
cyp19_no_omit <- cyp19 %>% rename(cyp19_SQ = mean_SQ) %>% dplyr::select(Sample, cyp19_SQ)
rpl8_no_omit <- rpl8 %>% rename(rpl8_SQ = mean_SQ) %>% dplyr::select(Sample, rpl8_SQ)
#merge data files
all_no_omit <- merge(amh_no_omit, rpl8_no_omit, by='Sample')
all_no_omit <- merge(all_no_omit, cyp19_no_omit, by='Sample')
#standardize AMH and CYP19A1 gene expression by that of housekeeping RPL8 gene
all_no_omit$st_amh <- all_no_omit$amh_SQ/all_no_omit$rpl8_SQ
all_no_omit$st_cyp19 <- all_no_omit$cyp19_SQ/all_no_omit$rpl8_SQ
all_no_omit<-all_no_omit%>%mutate(Dissection.Tube.ID=paste("FN", Sample, sep=""))
#write file with these values
save(all_no_omit, file="gene_expression_no_omit_data.RData")
