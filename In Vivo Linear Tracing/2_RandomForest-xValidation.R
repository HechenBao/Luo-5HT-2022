rm(list = ls())    #clear environment
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
    library(rstatix)
    library(ggpubr)
    library(heatmap.plus)
    library(randomForest)
    library(grid)
    library(ggplotify)
    library(randomForest)
    library(pROC)
  }
  
  # set working directory and load data
  {
    setwd('D:\\5HT_Data\\2020\\3way_test\\R_data_test')
    sert_data <- read.csv('5HT_data.csv')
    sert_data=sert_data[1:115,1:18]
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
    
    sert_data<-rename(sert_data,c("Type1.p"="X.Type1" ))
    sert_data<-rename(sert_data,c("Type2.p"="X.Type2" ))
    sert_data<-rename(sert_data,c("Type3.p"="X.Type3" ))
    
    sert_data<-rename(sert_data,c("K.Type1.p"="X.Ki67.type.1" ))
    sert_data<-rename(sert_data,c("K.Type2.p"="X.Ki67.type.2" ))
    sert_data<-rename(sert_data,c("K.Type3.p"="X.Ki67.type.3" ))
  }
  #str(sert_data)
  
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
  
  # normality test  
  
  
  
  # Seperate into 2 data sets: Basline/Stress and PBS/FLX
  {
    stress_data <- subset(sert_data,Treatment == "Baseline" | Treatment == "Stress")
    FLX_data    <- subset(sert_data,Treatment == "PBS" | Treatment == "FLX")
  }
  
  #reorder by group level
  {
    stress_data$Group<-factor(stress_data$Group,levels=c("F-WT-Baseline","F-WT-Stress","F-KO-Baseline","F-KO-Stress",
                                                         "M-WT-Baseline","M-WT-Stress","M-KO-Baseline","M-KO-Stress"))
    FLX_data$Group<-factor(FLX_data$Group,levels=c("F-WT-PBS","F-WT-FLX","F-KO-PBS","F-KO-FLX",
                                                   "M-WT-PBS","M-WT-FLX","M-KO-PBS","M-KO-FLX"))
  }
}  


## select data for random forest
data.rf<-select(sert_data,Gender,Genotype,Treatment,RFP,Ki67,Ki67.p,
                Type1,Type2,Type3,Type1.p,Type2.p,Type3.p)
# data.rf<-select(sert_data,Gender,Genotype,Treatment,
#                 Type1.p,Type2.p,Type3.p)
# Gender factor to numeric values ("F" to "0", "M" to "1")  
data.rf$Gender<-as.numeric(data.rf$Gender)
data.rf$Gender<-data.rf$Gender-1

# generate output array for RF-importance values and ROC values
{
  testing.n = 10  
  feature.n = ncol(data.rf)-1  
  ## for ROC curves
  roc.rf.matrix <- as.data.frame(matrix(nrow=1,ncol=3))
  colnames(roc.rf.matrix) <- c("tpp","fpp","rep")
  
  # auc
  auc.matrix <- as.data.frame(matrix(nrow=1,ncol=testing.n))
  rownames(auc.matrix) <- c("AUC")
  
  ## for featured variables
  imp.MSE.matrix <- as.data.frame(matrix(nrow=feature.n,ncol=testing.n))
  imp.INP.matrix <- as.data.frame(matrix(nrow=feature.n,ncol=testing.n))
  
}


### randomForest smoothing ROC
{
  for (i in 1:testing.n){
    data_set_size=floor(nrow(data.rf)*(5/10))
    index <- sample(1:nrow(data.rf),size=data_set_size)
    training <- data.rf[index,]
    testing <- data.rf[-index,]
    
    # bestmtry <- tuneRF(training,training$Genotype,
    #                    stepFactor = 1.2,improve=0.01,trace=T,plot=T)
    rf <- randomForest(Gender~.,data=training,mtry=3,ntree=1001,
                       importance=TRUE,proximity=TRUE)
    
    # rf.imp.p <- varImpPlot(rf,pch = 20, main = "Importance of Variables")
    imp <- importance(rf)
    imp.MSE.matrix[,i] <- imp[,-2]
    imp.INP.matrix[,i] <- imp[,-1]
    
    result <- data.frame(testing$Gender,predict(rf,testing[,2:ncol(testing)],type="response"))
    roc.rf.info <- roc(result$testing.Gender,result$predict.rf..testing...2.ncol.testing....type....response..,
                       ci=TRUE,plot=FALSE,legacy.axes=TRUE,auc=TRUE,smooth=TRUE)
    roc.rf.df <- data.frame(
      tpp=roc.rf.info$sensitivities*100, ## tpp = true positive percentage
      fpp=(1 - roc.rf.info$specificities)*100 ## fpp = false positive precentage
    )
    rep <- rep(i,nrow(roc.rf.df))
    roc.rf.df$rep <- rep
    roc.rf.matrix <- rbind(roc.rf.matrix,roc.rf.df)
    auc.matrix[,i] <- auc(roc.rf.info) 
  }
}

# organize output data     
{
  rownames(imp.MSE.matrix) <- rownames(imp)
  rownames(imp.INP.matrix) <- rownames(imp)
  
  roc.rf.matrix.p <- roc.rf.matrix[-1,]
  roc.rf.matrix.p$rep <- as.factor(roc.rf.matrix.p$rep)
  roc.rf.matrix.p$line <- rep(seq(1,nrow(roc.rf.df),1),testing.n)
}

# mean curve for roc
{
  roc.rf.mean <- as.data.frame(matrix(nrow=nrow(roc.rf.df),ncol=3))
  colnames(roc.rf.mean) <- c("tpp","fpp","rep")
  roc.rf.mean.tpp <- aggregate(roc.rf.matrix.p$tpp,list(roc.rf.matrix.p$line),mean)
  roc.rf.mean.fpp <- aggregate(roc.rf.matrix.p$fpp,list(roc.rf.matrix.p$line),mean)
  roc.rf.mean$tpp <- roc.rf.mean.tpp$x
  roc.rf.mean$fpp <- roc.rf.mean.fpp$x
  roc.rf.mean$rep <- as.factor(rep(100,nrow(roc.rf.df)))
}


roc.p <- ggplot(roc.rf.matrix.p,aes(x=fpp,y=tpp,group=rep))+
  geom_line(size=1,colour = 'steelblue',alpha=0.15)+
  labs(x="False Positive Rate", y="True Postive Rate")+
  ggtitle("ROC")+theme_classic()+theme(text=element_text(size=20))+
  coord_equal()+
  geom_line(data=roc.rf.mean,aes(x=fpp,y=tpp),size=1.5,colour='steelblue')+
  geom_abline(slope=1,size=1.2, color='red',linetype = "dashed")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))

roc.p
#####################################################################





### randomForest non-smoothing ROC
{
  {
    for (i in 1:10){  
      data_set_size=floor(nrow(data.rf)*(5/10))
      index <- sample(1:nrow(data.rf),size=data_set_size)
      training <- data.rf[index,]
      testing <- data.rf[-index,]
      
      rf <- randomForest(Gender~.,data=training,mtry=3,ntree=1001,
                         importance=TRUE,proximity=TRUE)
      rf.imp.p <- varImpPlot(rf,pch = 20, main = "Importance of Variables")
      imp <- importance(rf)
      imp.MSE.matrix[,i] <- imp[,-2]
      imp.INP.matrix[,i] <- imp[,-1]
      
      
      result <- data.frame(testing$Gender,predict(rf,testing[,2:ncol(testing)],type="response"))
      roc.rf.info <- roc(result$testing.Gender,result$predict.rf..testing...2.ncol.testing....type....response..,
                         plot=TRUE,legacy.axes=TRUE,
                         xlab="False Positive Percentage", ylab="True Postive Percentage",
                         col="#377eb8", lwd=2, auc=TRUE, print.auc=TRUE)
      roc.rf.df <- data.frame(
        tpp=roc.rf.info$sensitivities*100, ## tpp = true positive percentage
        fpp=(1 - roc.rf.info$specificities)*100 ## fpp = false positive precentage
      )
      rep <- rep(i,nrow(roc.rf.df))
      roc.rf.df$rep <- rep
      roc.rf.matrix <- rbind(roc.rf.matrix,roc.rf.df)
      auc.matrix[,i] <- auc(roc.rf.info) 
    }  
  }
  
  # plot all non-smoothing ROC
  roc.rf.matrix.p <- roc.rf.matrix[-1,]
  roc.rf.matrix.p$rep <- as.factor(roc.rf.matrix$rep)
  
  
  ggplot(roc.rf.matrix.p,aes(x=fpp,y=tpp,color=rep))+geom_line(size=1)+
    scale_color_brewer(palette="Blues")+
    labs(x="False Positive Percentage", y="True Postive Percentage")+
    ggtitle("ROC")+theme_classic()+theme(text=element_text(size=20))+
    coord_equal()
  
}




### modeling with logistic regression to predict Gender
glm.fit=glm(Gender~ .,data=data.rf, family="binomial")
par(pty = "s")
roc.info <- roc(data.rf$Gender, glm.fit$fitted.values, plot=TRUE,legacy.axes=TRUE,
                xlab="False Positive Percentage", ylab="True Postive Percentage",
                col="#377eb8", lwd=2)
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)
roc(data.rf$Gender, glm.fit$fitted.values, plot=TRUE,legacy.axes=TRUE,
    xlab="False Positive Percentage", ylab="True Postive Percentage",
    col="#377eb8", lwd=2, auc=TRUE, print.auc=TRUE) 

# roc(data.rf$Gender, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, 
#     print.auc=TRUE, print.auc.x=45, partial.auc=c(100, 90), 
#     auc.polygon = TRUE, auc.polygon.col = "#377eb822")

### overlay plot
roc(data.rf$Gender, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

plot.roc(result$testing.Gender,result$predict.rf..testing...2.ncol.testing....type....response.., percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Logisitic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)
