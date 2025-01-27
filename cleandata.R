# We need to pull the data from folders
# It needs to 
# 1. Go to the Label file as denoted by the group and datametric
# 2. For each label in there, look for a file that fits that label, datametric, and form
# 3. Pull the data from that file
# 4. Run a transpose on that data 
# 5. Append it as the next index of the output list
# 6. Once all files have been looped through, save the list as a file 
#Load Functions
library(MASS)
library(readr) # Used for reading CSV file faster than the built in r function
library(dplyr) # used for selecting columns in the read functions
library(caTools)
library(readxl)
remakeIDs<- function(){#Function to remake the IDs for use in the final function
  
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/")
  file_name <- paste0("iADRCautoPtxth_0.005_quant_OD.csv")
  # Read the data
  data <- read_delim(file_name, delim = ',', col_select = c(2,3,6)) # Change the delimiter to space
  # Define empty lists for each group
  fcnIDs <- list()
  fmciIDs <- list()
  fscdIDs <- list()
  mcnIDs <- list()
  mmciIDs <- list()
  mscdIDs <- list()
  combinedIDs <- list()
  
  # Loop to sort people into respective groups
  for (i in 1:nrow(data)) {
    person <- data[i, ]
    # Append to combined list
    combinedIDs <- append(combinedIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
    # Sort based on SEX and diag
    if (person$SEX == 2) {  # Female
      if (person$diag == 'CN' || person$diag == 'CN (IO)') {
        fcnIDs <- append(fcnIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
      } else if (person$diag == 'MCI') {
        fmciIDs <- append(fmciIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
      } else {
        fscdIDs <- append(fscdIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
      }
    } else {  # Male
      if (person$diag == 'CN' || person$diag == 'CN (IO)') {
        mcnIDs <- append(mcnIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
      } else if (person$diag == 'MCI') {
        mmciIDs <- append(mmciIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
      } else {
        mscdIDs <- append(mscdIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag)))
      }
    }
  }
  if (!dir.exists("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs")) dir.create("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs",recursive = TRUE)
  # When running, we want to change this to a DIR we can access
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs")
  
  saveRDS(fcnIDs,file = "fcnIDs.rds")
  saveRDS(fmciIDs,file = "fmciIDs.rds")
  saveRDS(fscdIDs,file = "fscdIDs.rds")
  saveRDS(mcnIDs,file = "mcnIDs.rds")
  saveRDS(mscdIDs,file = "mscdIDs.rds")
  saveRDS(mmciIDs,file = "mmciIDs.rds")
  saveRDS(combinedIDs,file = "combinedIDs.rds")
  
  
} 
sdaf<- function(group,datametric,form){ #save data as file, Used for processing raw dMRI data for use in the ABC model
  output <- list()
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs")
  file_name <- paste0(group,'IDs.rds') 
  IDs <- readRDS(file_name) 
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/SC_quant")
  for (i in IDs){
    filename <- paste0(i$`Subject ID`, "_SC_", datametric, "_", form, ".csv")
    if (file.exists(filename)) {
      mat <- read_delim(filename, delim = ' ', col_names = FALSE, show_col_types = FALSE)
      # Run the Transpose 
      for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
          if (i != j) {
            # Check if mat[i, j] is zero but mat[j, i] is non-zero
            if (mat[i, j] == 0 && mat[j, i] != 0) {
              mat[i, j] <- mat[j, i]
            }
            # Check if mat[j, i] is zero but mat[i, j] is non-zero
            else if (mat[j, i] == 0 && mat[i, j] != 0) {
              mat[j, i] <- mat[i, j]
            }
          }
        }
      }
    
      
      # Convert the list to a matrix
      m <- matrix(unlist(mat),ncol =(length(mat)))
      output[[length(output)+1]] = m
      }
  }
  final_list <- list()
  final_list[[1]] <-output
  o <- "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features"
  if (!dir.exists(o)) dir.create(o,recursive = TRUE)
  setwd(o)
  final_name = paste0(group,"_",datametric,"_",form,".rds")
  saveRDS(final_list,final_name)
}
pullfibers<- function(group,datametric){ #save data as file, Used for processing raw dMRI data for use in the ABC model
  output <- list()
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs")
  file_name <- paste0(group,'IDs.rds') 
  IDs <- readRDS(file_name) 
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/SC_quant")
  first <- IDs[[1]]
  for(i in IDs){
    filename <- paste0(IDs[[1]]$`Subject ID`, "_SC_", datametric, ".csv")
    if (file.exists(filename)) mat <- read_delim(filename, delim = ' ', col_names = FALSE, show_col_types = FALSE)
    else next
    
    # Run the Transpose 
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        if (i != j) {
          # Check if mat[i, j] is zero but mat[j, i] is non-zero
          if (mat[i, j] == 0 && mat[j, i] != 0) {
            mat[i, j] <- mat[j, i]
          }
          # Check if mat[j, i] is zero but mat[i, j] is non-zero
          else if (mat[j, i] == 0 && mat[i, j] != 0) {
            mat[j, i] <- mat[i, j]
          }
        }
      }
    }
    
    
    # Convert the list to a matrix
    m <- matrix(unlist(mat),ncol =(length(mat)))
    output[[length(output)+1]] = m
  }
  
  final_list <- list()
  final_list[[1]] <-output
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features")
  final_name = paste0(group,"_",datametric,".rds")
  saveRDS(final_list,final_name)
}
makeFiles <- function(grouplist,datametriclist,formlist){

for (i in grouplist){
  for (d in datametriclist){
    for (f in formlist){
      sdaf(i,d,f)
      }
    }
  }
}
# Pull the Labels
pullYs<-function(group){
  # Our Data is held within the iADRC_Struture_Diffusion_Tau_Abeta_84ROIcombo.xlsx file
  # The data has the following Columns we care about
  # 1. Study_ID: Holds the ID for that person
  # 2. Tau : 84 different columns that hold the Tau data, will be Col 1 of our final file
  # 3. Amyloid: Same Structure as Tau, will be col 2
  # Make what will become our final output list
  # For each label, make a new index in that list that is a 84x2 matrix 
  #For each Label we have, pull all the Tau, then all of the Amyloid Data 
  # Do this for each label in the group
  
  #Read The IDs for our Group, we can use the old deprecated IDs files
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs")
  filename <- paste0(group,'IDs.rds')
  Ids <- readRDS(filename) # Only Column 1 is important for this 
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/")
  fulldata <- read_xlsx("iADRC_Struture_Diffusion_Tau_Abeta_84ROIcombo.xlsx")
  # Drop the Unneeded columns 
  # We only need Study_ID, SCroi(1-84)_Tau, SCroi(1-84)_Amyloid (col 1, and 1147 - 1314)
  data <- fulldata[,c(1,1147:1314)]
  rm(fulldata) # To save memory
  # In Theory, this should have loaded just the columns for IDs, Tau, and Amyloid Data 
  output <- list()
  for(i in Ids){
    # Find if the current ID is in the excel data
    curID <- i$`Subject ID`
    index <- which(data$Study_ID ==curID)
    # If it is
    if (length(index) == 0) {
      next
    } else {
      tau_data <- data[index,2:85]
      amy_data <- data[index,86:169]
      # Convert these to a 84x1 matrix each
      tau_matrix <- as.matrix(t(tau_data))
      amy_matrix <- as.matrix(t(amy_data))
      # combine them together into 1 84x2 matrix
      final_matrix <- cbind(tau_matrix,amy_matrix)
      # append them to the output list 
      output[[length(output) + 1]] <- final_matrix
    }
  }
  o <- "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/"
  if (!dir.exists(o)) dir.create(o,recursive = TRUE)
  setwd(o)
  finalname <- paste0(group,'label.rds')
  saveRDS(output,finalname)
}

grouplist <- list('fcn','fmci','fscd','mcn','mmci','mscd') # Before we run the whole group, I need to fix the names, might mean rerunning the data creation scripts 
datametriclist <- list('Da','Dr','MD','FA','ICVF','OD') # We need to make a new version of the script for mean length
formlist <-list('mean','median')
makeFiles(grouplist,datametriclist,formlist)


for(g in grouplist){ 
  sdaf(g,'OD','mean')
  pullYs(g)
}


