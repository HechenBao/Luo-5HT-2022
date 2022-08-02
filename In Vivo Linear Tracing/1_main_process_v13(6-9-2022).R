## This script starts with data cleaning and re-define of catagories, followed by normalization calculation for genotype or treatment.

rm(list = ls())    #clear environment
#sert_data <- read.csv(file.choose())
{
  #load library
  {
    library(oro.nifti)
    library(neurobase)
    library(RColorBrewer)
    library(wesanderson)
    library(scales)
    library(ggplot2)
    library(gplots)
    library(reshape)
    library(reshape2)
    library(ppcor)
    library(plyr)
    library(dplyr)
    library(ggpubr)
    library(factoextra)
    library(corrplot)
  }
  
  # set working directory and load data
  {
    setwd('D:\\5HT_Data\\summarized_data\\5HT_data')
    sert_data <- read.csv('5HT_data_v5.csv')
  }
  
  #rename variables
  {
    # sert_data<-rename(sert_data,c("Gender"="Gender..1.M.0.F."))
    # sert_data<-rename(sert_data,c("Genotype"="Geno.code..5.WT.3.KO." ))
    # sert_data<-rename(sert_data,c("Treatment"="Treatment..1.Baseline.2.Stress.3.PBS..4.FLX." ))
    
    sert_data<-rename(sert_data,c("RFP"="RFP" )) # needs to update by Yanjia 6/9/22 
    
    sert_data<-rename(sert_data,c("Sox2"="sox2." ))
    sert_data<-rename(sert_data,c("Ki67"="Total..Ki67.RFP.." ))

    sert_data<-rename(sert_data,c("K.Sox2"="Ki67.Sox2." ))
    sert_data<-rename(sert_data,c("K.Sox2neg"="ki67.sox2..neg" ))

    sert_data<-rename(sert_data,c("Type1"="type1" ))
    sert_data<-rename(sert_data,c("Type2a"="type2a" ))
        
    sert_data<-rename(sert_data,c("K.Type1"="Ki67.Type1..ki67.sox.2.GFAP.." ))
    sert_data<-rename(sert_data,c("K.Type2a"="Ki67.Type2a..ki67.sox2.GFAP.." ))

    sert_data<-rename(sert_data,c("K.Type1.p"="X.ki67type1.type1" ))
    sert_data<-rename(sert_data,c("K.Type2a.p"="X.ki67type2a.type2a" ))
    sert_data<-rename(sert_data,c("K.Sox2.p"="X.ki67sox2.sox2" ))

    # sert_data<-rename(sert_data,c("K.Type1.p"="X.Ki67.type.1" ))
    # sert_data<-rename(sert_data,c("K.Type2.p"="X.Ki67.type.2" ))
    # sert_data<-rename(sert_data,c("K.Type3.p"="X.Ki67.type.3" ))
    # 
    # sert_data<-rename(sert_data,c("pyk.type1.p"="pyk.type.1" ))
    # sert_data<-rename(sert_data,c("pyk.type2.p"="pyk.type.2" ))
    # sert_data<-rename(sert_data,c("pyk.type3.p"="pyk.type.3" ))

  }
  
  #define factors
  {
    sert_data$Gender<-as.factor(sert_data$Gender)
    sert_data$Gender<-revalue(sert_data$Gender, c("Female" = "F"))
    sert_data$Gender<-revalue(sert_data$Gender, c("Male" = "M"))
    
    sert_data$Genotype<-as.factor(sert_data$Genotype)
    # sert_data$Genotype<-revalue(sert_data$Genotype, c("5" = "WT"))
    # sert_data$Genotype<-revalue(sert_data$Genotype, c("3" = "KO"))
    
    sert_data$Treatment<-as.factor(sert_data$Treatment)
    # sert_data$Treatment<-revalue(sert_data$Treatment, c("1" = "Baseline"))
    # sert_data$Treatment<-revalue(sert_data$Treatment, c("2" = "Stress"))
    # sert_data$Treatment<-revalue(sert_data$Treatment, c("3" = "PBS"))
    # sert_data$Treatment<-revalue(sert_data$Treatment, c("4" = "FLX"))
  }
  
  #define groups
  {
    sert_data$Group=as.factor(paste(paste(as.character(sert_data$Gender),
                                          as.character(sert_data$Genotype),sep="-"),
                                    as.character(sert_data$Treatment), sep ="-"))
    head(sert_data)
  }
  
  #calculate new variables
  {
    sert_data$Sox2.neg<-sert_data$RFP-sert_data$Sox2
    sert_data$Sox2.neg.p<-sert_data$Sox2.neg/sert_data$RFP*100

  }
  
  #subtract useful variables
  # {
  #   sert_data<-select(sert_data,Group,Gender,Genotype,Treatment,Animal.ID,
  #                     RFP,Sox2,Type1,Type2a,Sox2.neg,Ki67,
  #                     Type1.p,Type2a.p,Sox2.neg.p,Ki67.p,Sox2.p,sertVol.p2)  
  # }
  {
    sert_data<-select(sert_data,Group,Gender,Genotype,Treatment,Animal.ID,
                      RFP,Sox2,Sox2.neg,Ki67,Type1,Type2a,K.Type1,K.Type2a,K.Sox2,K.Sox2neg,K.Type1.p,K.Type2a.p,K.Sox2.p)  
  
    var_n = ncol(sert_data)-5
    
  # Add a new column as Gender+Geno (GGType)
  {
    sert_data$GGType <- paste(sert_data$Gender, sert_data$Genotype, sep="_")
    unique(sert_data$Group)
  }
  
  var_agg_n = var_n+5
  var_mean_n = var_n+1
  
  # Get mean value for each group for Baseline and PBS
  {
    baseline_data <- subset(sert_data,  Treatment == 'Baseline')
    baseline_mean <- aggregate(baseline_data[,6:var_agg_n], list(baseline_data$GGType), mean)
    
    # PBS_data <- subset(sert_data,  Treatment == 'PBS')
    # PBS_mean <- aggregate(PBS_data[,6:10], list(PBS_data$GGType), mean)
  }
  
  # Seperate into 2 data sets: Basline/Stress and PBS/FLX
  {
    stress_data <- subset(sert_data,Treatment == "Baseline" | Treatment == "Stress")
    # FLX_data    <- subset(sert_data,Treatment == "PBS" | Treatment == "FLX")
  }
  
  #reorder by group level
  {
    stress_data$Group<-factor(stress_data$Group,levels=c("F-WT-Baseline","F-WT-Stress","F-KO-Baseline","F-KO-Stress",
                                                         "M-WT-Baseline","M-WT-Stress","M-KO-Baseline","M-KO-Stress"))
    # FLX_data$Group<-factor(FLX_data$Group,levels=c("F-WT-PBS","F-WT-FLX","F-KO-PBS","F-KO-FLX",
    #                                                "M-WT-PBS","M-WT-FLX","M-KO-PBS","M-KO-FLX"))
  }
  
  
  ## normalized to individual conditions for the same genotype (regress treatment)
  # Define Functions for normalization
    {
    # function: Normalized score
    norm.sc <- function(newData, inData, meanTable) {
      newData <- inData
      for (i in 1:nrow(inData)){
        for (j in 1:nrow(meanTable)){
          if (inData[i,"GGType"] == meanTable[j,1]) {
            newData[i,6:var_agg_n] <- (inData[i,6:var_agg_n] - meanTable[j,2:var_mean_n]) / meanTable[j,2:var_mean_n]
          }
        }
      }
      return(newData)
    }
    
  }
  
  # Call the functions
  stress_norm.sc <- norm.sc(stress_norm.sc, stress_data, baseline_mean)
  # FLX_norm.sc <- norm.sc(FLX_norm.sc, FLX_data, PBS_mean)
  #sert_norm.sc <- rbind(stress_norm.sc,FLX_norm.sc)
    
  
  # divide raw dataset by 10000 for export plots
  {
  stress_data_2<-stress_data
  stress_data_2[,6:15]<-stress_data_2[,6:15]/10000
  
  # FLX_data_2<-FLX_data
  # FLX_data_2[,6:9]<-FLX_data_2[,6:9]/10000
  }
}
}

## normalized score for baseline data (KO normalized to WT)
# Get mean value for each group for Baseline and PBS
{
  WT_data <- subset(sert_data,  Treatment == 'Baseline'& Genotype == 'WT')
  WT_mean <- aggregate(WT_data[,6:var_agg_n], list(WT_data$Gender), mean)
}

WTKO_data<-subset(stress_data, Treatment == 'Baseline')

# function: Normalized score
base.norm.sc <- function(newData, inData, meanTable) {
  newData <- inData
  for (i in 1:nrow(inData)){
    for (j in 1:nrow(meanTable)){
      if (inData[i,"Gender"] == meanTable[j,1]) {
        newData[i,6:var_agg_n] <- (inData[i,6:var_agg_n] - meanTable[j,2:var_mean_n])/meanTable[j,2:var_mean_n]
      }
    }
  }
  return(newData)
}
{
WTKO_norm.sc <- base.norm.sc(WTKO_norm.sc, WTKO_data, WT_mean)
WTKO_norm.sc$Group<-factor(WTKO_norm.sc$Group,levels=c("F-WT-Baseline","F-KO-Baseline",
                                                       "M-WT-Baseline","M-KO-Baseline"))
}


write.csv(stress_data_2, 'stress_data_v13.csv')
write.csv(stress_norm.sc, 'stress_norm_v13.sc.csv')
write.csv(WTKO_norm.sc, 'stress_WTKOnorm.sc_v13.csv')



######################################################################################
## Exploratory heatmap plots for the dataset

## measure mean of stress/FLX dataset 
stress_norm.sc_ratio <- subset(stress_norm.sc,Treatment == "Stress")
stress.norm.sc_mean <- aggregate(stress_norm.sc_ratio[,4:14], list(stress_norm.sc_ratio$Group), mean)

FLX_norm.sc_ratio <- subset(FLX_norm.sc,Treatment == "FLX")
FLX.norm.sc_mean <- aggregate(FLX_norm.sc_ratio[,4:14], list(FLX_norm.sc_ratio$Group), mean)

## heatmap for combined dataset
rbinded_mean <- rbind(stress.norm.sc_mean,FLX.norm.sc_mean)
row.names(rbinded_mean) <- c("F-WT-Stress","F-KO-Stress","M-WT-Stress","M-KO-Stress",
                             "F-WT-FLX","F-KO-FLX","M-WT-FLX","M-KO-FLX")
rbinded_mean <- rbinded_mean[,2:12]
h <- melt(as.matrix(rbinded_mean))
colnames(h) <- c("Group","var","value")
h$Group <- factor(h$Group,levels=c("M-KO-FLX","M-WT-FLX","F-KO-FLX","F-WT-FLX",
                                   "M-KO-Stress","M-WT-Stress","F-KO-Stress","F-WT-Stress"))
h1 <- subset(h,var == "Type1"|var == "Type2"|var == "Type3")
#h$Group <- factor(h$Group,levels=c("F-WT-Stress","F-KO-Stress","M-WT-Stress","M-KO-Stress",
#                                  "F-WT-FLX","F-KO-FLX","M-WT-FLX","M-KO-FLX")) 
ggplot()+geom_tile(aes(x=h1$var,y=h1$Group,fill=h1$value))+
  scale_fill_gradient2(low=("deepskyblue4"),mid="white",high=("firebrick1"),midpoint=0)
ggplot()+geom_tile(aes(x=h1$var,y=h1$Group,fill=h1$value))+
  scale_fill_gradient2(low=("deepskyblue4"),mid="white",high=("firebrick1"),midpoint=0)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme_classic() +
  coord_equal() 

h2 <- subset(h,var == "RFP"|var == "Type1"|var == "Type2"|var == "Type3"|var == "Ki67")
ggplot()+geom_tile(aes(x=h2$var,y=h2$Group,fill=h2$value))+
  scale_fill_gradient2(low=("deepskyblue4"),mid="white",high=("firebrick1"),midpoint=0)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme_classic()+
  coord_equal() 


a<-as.matrix(rbinded_mean)
heatmap.2(a,trace="none")
