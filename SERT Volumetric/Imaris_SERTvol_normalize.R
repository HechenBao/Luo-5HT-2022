rm(list = ls())    #clear environment

{
  library("pdftools")
  library("ggplot2")
  library("tidyverse")
  library("data.table")
}

# set working directory and load data
{
  setwd("D:\\5HT_Data\\summarized_data\\SERT volumetric\\IMARIS-SERT\\SERT_analysis\\SERT_v2\\Baseline_control")
  sert_vol <- read.csv('Females_measures.csv')
  sert_ROI <- read.csv('Females_ROIs.csv')
}

{
rownames(sert_vol) <- sert_vol[,1] #change row names as animal ID
rownames(sert_ROI) <- sert_ROI[,1] #change row names as animal ID

sert_vol[,1] <- NULL #remove the first column
sert_ROI[,1] <- NULL #remove the first column
}

norm_matrix <- sert_vol/sert_ROI*100

norm_matrix$means <- rowMeans(norm_matrix,na.rm = TRUE)

write.csv(norm_matrix, 'norm_matrix2.csv')


