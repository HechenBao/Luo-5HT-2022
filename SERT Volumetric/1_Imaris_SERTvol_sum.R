{
library("pdftools")
library("ggplot2")
library("tidyverse")
library("data.table")
}

#Before using this code make sure your files are named ("A# section #") so the code can read the excel files. Also make sure that all the volume files are within the working directory.
#set working directory to where you want your data to be saved edit for your own data and create a folder named "ROI Folder" that has all the ROI's.
#At the end there will be a csv table output with the sums of each object for each section named the same as the file selected.
#There should be no 0 values in the output table unless a section didn't have data to input.
{
datadirectory <- ("D:/5HT_Data/summarized_data/SERT volumetric/IMARIS-SERT/SERT_analysis/SERT_v2") ##edit directory accordingly requires single backslash for windows directory change accordingly
setwd(datadirectory)

#read in excel file by selecting csv file to read in data this will open a new window
path = choose.files()
paths <- unlist(strsplit(path,"\\\\"))
paths <- tail(paths,1)

#If then statement to determine if file is csv or excel.
csvread <- fread(file=path, skip = 2, header = TRUE,fill = F)

#make sure that excel file from above has removed the name columns if excel file is unedited from imaris, this part should work fine
volumecol <- csvread$Volume
sectionnamecol <- csvread$`Original Image Name`
}

{
#Input will be requested in the console for 2 things, total animals and total sections. Input is required in console for code to work!
totalanimals <- "Enter the total number of animals to analyze:" %>% readline() %>% as.integer()
totalsections <- "Enter the total number of sections to analyze per animal:" %>% readline() %>% as.integer()
print(paste("You have",paste(totalanimals,paste("animals and",paste(totalsections,"sections per animal")))))


#function making animal number and pasting ".*" so it can search all text until it reaches section #
#make sure animals in data set are labeled as animal # section # ex:("a1 section 1"...) since it will search for that text

animallist <- array(1:totalanimals)
animalID <- paste(unique(paste(substr(sectionnamecol,1,5))),".*",sep = "")

#forloop making list of sections from entered section data
sectionlist <- array()
for (i in 1:totalsections){
   sectionlist[i] <- (paste("G00",i,sep = ""))


}

#create empty matrix for forloop
meanmatrix <- matrix(nrow = length(animallist),ncol = length(sectionlist))

#make large forloop that goes through animallist and sectionlist and sums up the values for each animal from all section numbers
for (i in 1: length(animallist))
  {
  for (j in 1:length(sectionlist)){
  animalsecpaste<- paste(animalID[i],sectionlist[j],sep="")
  searchpastedtext<- grepl(animalsecpaste,sectionnamecol,ignore.case = TRUE)
  meanmatrix[i,j] <- sum(as.numeric(csvread$Volume[searchpastedtext]))
    }
}

#transpose matrix to have columns be animals and sections be rows
meanmatrix <- as.data.frame((meanmatrix))

#make array to only have animal number in name list
splitname <- function(x){
  (strsplit(x,".*", fixed=TRUE))}

animalnames <- unlist(lapply(animalID, splitname))

#name mean matrix columns and rows
colnames(meanmatrix) <- sectionlist
rownames(meanmatrix) <- animalnames
}
#visualize data, look for empty bars, indicates missing data
#pdf(file = paste0(datadirectory,"Check Volumes.pdf"))
barplot(as.matrix(t(meanmatrix)), xlab = "Animals", ylab = "Volume of section",main = "Check that volumes of animals are present", col = c("azure","azure1","azure2","azure3","lightblue","pink","yellow" ))
#dev.off()

#prints and exports mean matrix table that sums up all the objects, the file has the prefix matrix output added to the file name
print(meanmatrix)

#change name for up/low ROIs
write.csv(meanmatrix, file = paste("object_matrix_output",paths))

########## ROI IMPORT ############
# To use this section select a folder that has all the ROI's csv files it should only have the ROI's for data set you are analyzing. Any other files will disrupt the process. 
# The ROI volumes should be named "A# section #"
# Furthermore the ROI files should have 3 lines as header and then the ROI volume for the section. i.e( row1 =blank, row2=blank, row3= "Volume" row4=(Volume of ROI))

#create a list of the files from your target directory and sort by number
{
roidirectory <- choose.dir() %>% setwd()

#create a list of the files from your target directory and sort by number
file_list <- list.files() %>% str_sort(numeric = TRUE)

#confirm files list is in order from decending to ascending
print(file_list)

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
roimatrix <- matrix(nrow = length((animallist)),ncol = length(sectionlist))

#Function to read roi file and extract sum of volume only
readroi <- function(roi){
  tempread2 <- fread(file_list[roi],header = T,skip = 3)
  value <- as.numeric(sum(tempread2$Volume))
}

#Loop through animals and sections checking for a value in original data file. If no value enter "0" in roimatrix.
#If value exist then enter value in roimatrix.

for (i in 1: length(animallist)){
  for(j in 1: length(sectionlist)){
    animalsecpaste <- paste(animalID[i],sectionlist[j],sep="")
    searchpastedtext <-  grep(animalsecpaste ,csvread$`Original Image Name`,ignore.case = T)
      if(sum(searchpastedtext) > 0){
        searchanimallist<- grep(animalsecpaste,file_list,ignore.case =T)
        tempread <- readroi(searchanimallist)
        roimatrix[i,j] <- as.numeric(tempread)
      }
      else{
        roimatrix[i,j] <- as.numeric(0)
      }
    }
  }
}


# Change dimensions of Roimatrix
{
roimatrix <- as.data.frame(roimatrix)

#name ROI Matrix
colnames(roimatrix) <- sectionlist
rownames(roimatrix) <- animalnames

#Print ROI matrix and verify that values are within the matrix
print(roimatrix)

write.csv(roimatrix, file = paste("ROImatrix_output",paths))
}

### Matrix divisions to normalize object volume ROI by volume of ROI ###


normalizedmatrix <- (meanmatrix/roimatrix)*100

#set NA to 0 values in matrix
normalizedmatrix[is.na(normalizedmatrix)] <- 0

#view values in normalized matrix
print(normalizedmatrix)


#name normalized ROI
colnames(normalizedmatrix) <- animalnames
rownames(normalizedmatrix) <- sectionlist

sumnormalizedmatrix <- sapply(normalizedmatrix, sum)

#pdf(paste0(directory,"Density of Measured Objects.pdf"))

barplot(as.matrix(normalizedmatrix), xlab = "Animals", ylab = "% Object volume/ROI Volume",main = "Density of Measured Objects", col = c("azure","azure1","azure2","azure3","lightblue","pink","yellow" ))

normalizedmatrix <- rbind(normalizedmatrix,sumnormalizedmatrix)

#dev.off()

print(normalizedmatrix)

#output as csv file Normalized matrix file = object_volume/ROI_volume this can be graphed directly
write.csv(normalizedmatrix, file = paste0(datadirectory,paste("normalizedmatrix_output",paths)))

