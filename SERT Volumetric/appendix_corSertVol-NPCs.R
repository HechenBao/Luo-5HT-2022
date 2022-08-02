### This script is used to measure the correlation of SERT volumetric and density of different neural precursors.

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
  sert_data <- read.csv('5HT_data.csv')
}

#rename variables
{
  sert_data<-rename(sert_data,c("Gender"="Gender..1.M.0.F."))
  sert_data<-rename(sert_data,c("Genotype"="Geno.code..5.WT.3.KO." ))
  sert_data<-rename(sert_data,c("Treatment"="Treatment..1.Baseline.2.Stress.3.PBS..4.FLX." ))
  
  sert_data<-rename(sert_data,c("RFP"="RFP." ))
  sert_data<-rename(sert_data,c("Ki67"="Ki67.RFP." ))
  sert_data<-rename(sert_data,c("Sox2"="Sox2.RFP." ))
  
  sert_data<-rename(sert_data,c("Type1"="Type1..GFAP.Sox2.RFP.." ))
  sert_data<-rename(sert_data,c("Type2"="Type2..GFAP.Sox2.RFP.." ))
  sert_data<-rename(sert_data,c("Type3"="Type3..Sox2.." ))
  
  sert_data<-rename(sert_data,c("Ki67.p"="X..Ki67.RFP.." ))
  sert_data<-rename(sert_data,c("Sox2.p"="X..Sox2.RFP.." ))
  sert_data<-rename(sert_data,c("K.Sox2"="Sox2.Ki67." ))
  sert_data<-rename(sert_data,c("K.Sox2.p"="X.Ki67Sox2.Sox2" ))
  
  sert_data<-rename(sert_data,c("Type1.p"="X.Type1" ))
  sert_data<-rename(sert_data,c("Type2.p"="X.Type2" ))
  sert_data<-rename(sert_data,c("Type3.p"="X.Type3" ))
  
  sert_data<-rename(sert_data,c("K.Type1.p"="X.Ki67.type.1" ))
  sert_data<-rename(sert_data,c("K.Type2.p"="X.Ki67.type.2" ))
  sert_data<-rename(sert_data,c("K.Type3.p"="X.Ki67.type.3" ))
  
  sert_data<-rename(sert_data,c("pyk.type1.p"="pyk.type.1" ))
  sert_data<-rename(sert_data,c("pyk.type2.p"="pyk.type.2" ))
  sert_data<-rename(sert_data,c("pyk.type3.p"="pyk.type.3" ))
  
}

#define factors
{
  sert_data$Gender<-as.factor(sert_data$Gender)
  sert_data$Gender<-revalue(sert_data$Gender, c("0" = "F"))
  sert_data$Gender<-revalue(sert_data$Gender, c("1" = "M"))
  
  sert_data$Genotype<-as.factor(sert_data$Genotype)
  sert_data$Genotype<-revalue(sert_data$Genotype, c("5" = "WT"))
  sert_data$Genotype<-revalue(sert_data$Genotype, c("3" = "KO"))
  
  sert_data$Treatment<-as.factor(sert_data$Treatment)
  sert_data$Treatment<-revalue(sert_data$Treatment, c("1" = "Baseline"))
  sert_data$Treatment<-revalue(sert_data$Treatment, c("2" = "Stress"))
  sert_data$Treatment<-revalue(sert_data$Treatment, c("3" = "PBS"))
  sert_data$Treatment<-revalue(sert_data$Treatment, c("4" = "FLX"))
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
  sert_data$Type2a<-sert_data$Sox2-sert_data$Type1
  sert_data$Sox2.neg<-sert_data$RFP-sert_data$Sox2
  sert_data$Sox2.neg.p<-sert_data$Sox2.neg/sert_data$RFP*100
  sert_data$Type2a.p<-sert_data$Type2a/sert_data$RFP*100
}

#subtract useful variables
# {
#   sert_data<-select(sert_data,Group,Gender,Genotype,Treatment,Animal.ID,
#                     RFP,Sox2,Type1,Type2a,Sox2.neg,Ki67,
#                     Type1.p,Type2a.p,Sox2.neg.p,Ki67.p,Sox2.p,sertVol.p2)  
# }
{
  sert_data<-select(sert_data,Group,Gender,Genotype,Treatment,Animal.ID,
                    Sox2,Type1,Type2a,K.Sox2,K.Sox2.p,sertVol.p2)  
  
  # Add a new column as Gender+Geno (GGType)
  {
    sert_data$GGType <- paste(sert_data$Gender, sert_data$Genotype, sep="_")
    unique(sert_data$Group)
  }
  
  
  
  # Get mean value for each group for Baseline and PBS
  {
    baseline_data <- subset(sert_data,  Treatment == 'Baseline')
    baseline_mean <- aggregate(baseline_data[,6:10], list(baseline_data$GGType), mean)
    
    PBS_data <- subset(sert_data,  Treatment == 'PBS')
    PBS_mean <- aggregate(PBS_data[,6:10], list(PBS_data$GGType), mean)
  }
  
}
}


## cell density
# baseline by group
baseline_data_2<-baseline_data
baseline_data_2[,6:9]<-baseline_data_2[,6:9]/10000
{
g1<-ggscatter(baseline_data_2, x = "Sox2", y = "sertVol.p2",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Sox2", ylab = "sertVol%")+ facet_grid(~Group)+ 
          theme_classic()
g2<-ggscatter(baseline_data_2, x = "Type1", y = "sertVol.p2",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Type1", ylab = "sertVol%")+ facet_grid(~Group)+ 
          theme_classic()
g3<-ggscatter(baseline_data_2, x = "Type2a", y = "sertVol.p2",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Type2a", ylab = "sertVol%")+ facet_grid(~Group)+ 
          theme_classic()
g4<-ggscatter(baseline_data_2, x = "K.Sox2", y = "sertVol.p2",
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Ki67 Sox2", ylab = "sertVol%")+ facet_grid(~Group)+ 
  theme_classic()
g5<-ggscatter(baseline_data_2, x = "K.Sox2.p", y = "sertVol.p2",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Ki67 Sox2.p", ylab = "sertVol%")+ facet_grid(~Group)+ 
          theme_classic()


SertVol.group.plot<-ggarrange(g1, g2, g3, g4, g5,
                           labels = c("A", "B", "C", "D","E"),
                           ncol = 1, nrow = 5)
SertVol.group.plot
}

# baseline by genotype
{
  g1<-ggscatter(baseline_data_2, x = "RFP", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "RFP", ylab = "sertVol%")+ facet_grid(~Genotype)+ 
    theme_classic()
  g2<-ggscatter(baseline_data_2, x = "Type1", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Type1", ylab = "sertVol%")+ facet_grid(~Genotype)+ 
    theme_classic()
  g3<-ggscatter(baseline_data_2, x = "Type2a", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Type2a", ylab = "sertVol%")+ facet_grid(~Genotype)+ 
    theme_classic()
  g4<-ggscatter(baseline_data_2, x = "Ki67", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Ki67", ylab = "sertVol%")+ facet_grid(~Genotype)+ 
    theme_classic()
  SertVol.group.plot<-ggarrange(g1, g2, g3, g4,
                                labels = c("A", "B", "C", "D"),
                                ncol = 1, nrow = 4)
  SertVol.group.plot
}

# baseline by sex
{
  g1<-ggscatter(baseline_data_2, x = "RFP", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "RFP", ylab = "sertVol%")+ facet_grid(~Gender)+ 
    theme_classic()
  g2<-ggscatter(baseline_data_2, x = "Type1", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Type1", ylab = "sertVol%")+ facet_grid(~Gender)+ 
    theme_classic()
  g3<-ggscatter(baseline_data_2, x = "Type2a", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Type2a", ylab = "sertVol%")+ facet_grid(~Gender)+ 
    theme_classic()
  g4<-ggscatter(baseline_data_2, x = "Ki67", y = "sertVol.p2",
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Ki67", ylab = "sertVol%")+ facet_grid(~Gender)+ 
    theme_classic()
  SertVol.group.plot<-ggarrange(g1, g2, g3, g4,
                                labels = c("A", "B", "C", "D"),
                                ncol = 1, nrow = 4)
  SertVol.group.plot
}

