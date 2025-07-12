#cleandata.R is a colleciton of four functions that handle data processing to be used in the ABC model 
#makeIDs - Used to create files that store the IDs for each of our 7 groups 
#savedata - For all datametrics, pulls the raw dMRI data for each participant, and saves it as a file
#pullYs - pulls the Tau and Amyloid data for use in model 
library(MASS)
library(readr) # Used for reading CSV file faster than the built in r function
library(dplyr) # used for selecting columns in the read functions
library(caTools)
library(readxl)
makeIDs<- function(dataloc="/N/u/conlcorn/BigRed200/SexLinkedProject/data/",
                   outputloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs",
                   file_name = "iADRCautoPtxth_0.005_quant_OD.csv" ){#Function to make the IDs for use in the final function
  setwd(dataloc)
  # Read the data
  data <- read_delim(file_name, delim = ',', col_select = c(2,3,6)) # Read the collums, only the ones that we need 
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

  if (!dir.exists(outputloc)) dir.create(outputloc,recursive = TRUE)
  setwd(outputloc)

  saveRDS(fcnIDs,file = "fcnIDs.rds")
  saveRDS(fmciIDs,file = "fmciIDs.rds")
  saveRDS(fscdIDs,file = "fscdIDs.rds")
  saveRDS(mcnIDs,file = "mcnIDs.rds")
  saveRDS(mscdIDs,file = "mscdIDs.rds")
  saveRDS(mmciIDs,file = "mmciIDs.rds")
  saveRDS(combinedIDs,file = "combinedIDs.rds")
} 
savedata<- function(group = 'combined',
                    datametric = 'OD',
                    form = 'mean', 
                    IDLoc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs" ,
                    dataloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/SC_quant",
                    outputloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features" ){ #save data as file, Used for processing raw dMRI data for use in the ABC model
  output <- list()
  setwd(IDLoc)
  file_name <- paste0(group,'IDs.rds') 
  IDs <- readRDS(file_name) 
  setwd(dataloc)
  for (i in IDs){
    if(datametric == 'meanlength' || datametric == 'numoffibers') filename <- paste0(i$`Subject ID`, "_SC_", datametric, ".csv")
    else filename <- paste0(i$`Subject ID`, "_SC_", datametric, "_", form, ".csv")

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
  o <- outputloc 
  if (!dir.exists(o)) dir.create(o,recursive = TRUE)
  setwd(o)
  final_name = paste0(group,"_",datametric,"_",form,".rds")
  saveRDS(final_list,final_name)
}
# Pull the Labels
pullYs<-function(group='fcn',
                 IDLoc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/IDs",
                 dataloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/",
                 dataname = "iADRC_Struture_Diffusion_Tau_Abeta_84ROIcombo.xlsx",
                 outputloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels/"){
  #Read The IDs for our Group, we can use the old deprecated IDs files
  setwd(IDLoc)
  filename <- paste0(group,'IDs.rds')
  Ids <- readRDS(filename) # Only Column 1 is important for this 
  setwd(dataloc)
  fulldata <- read_xlsx(dataname)
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
  o <- outputloc
  if (!dir.exists(o)) dir.create(o,recursive = TRUE)
  setwd(o)
  finalname <- paste0(group,'label.rds')
  saveRDS(output,finalname)
}



remakeAtlas <- function()
{
  library(dplyr)
  library(stringr)
  # Remakes the brain atlas file to fix formating 
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/")
  #pull the two atlases we need to merge 
  
  rightindex <- read_excel("iADRC_Struture_Diffusion_Tau_Abeta_84ROIcombo.xlsx", sheet = 2)
  right <- rightindex[10:nrow(rightindex),]
  regions <- read_xlsx("Table_84ROI_Yeo7_relabeled_withlobes.xlsx")
  
  #Clean the names 
  
  ROI <- right %>% mutate(ROI = str_replace(ICV, "^ctx-", ""))  # Removes "left_" or "right_" from the start of the string
  
  combined_df <- ROI %>%
    mutate(row_number = row_number()) %>%  # Add row number to 'right'
    left_join(regions, by = "ROI")  # Join on cleaned ICV names
  
  final_df <- combined_df %>%
    dplyr::select(ROI, row_number, LOBE)
  # save the df as a file to use later 
  save(final_df, file = "finalAtlas.rds")
  
}
grouplist <- list("fcn", "fmci", 
                   "mcn", "mmci"
)
metriclist <- list("Da","Dr","FA","ICVF","MD")

for(g in grouplist) for (m in metriclist) savedata(group = g, datametric = m)


UVC2UVPM <- function(group,limit = NA, modeldir, outputdir){
# Recalculates the UVPM from a given UVC file, more robust than original UVPM

setwd(modeldir)

#load the UVC file
uvc_file <- paste0('Average_',group,'_UVC.rds')
if (!file.exists(uvc_file)) {
  stop(paste("UVC file for group", group, "does not exist."))
}
uvc <- readRDS(uvc_file)

if(is.na(limit)) {
  # If limit is not specified, use the full UVC
  UVC <- uvc
} else{ 
  UVC <- uvc[limit,]
}

# create a 84x84 matrix to hold the UVPM
uvpm <- matrix(0, nrow = 84, ncol = 84)

  loc_matrix <- matrix(0, nrow = 84, ncol = 84)
  index <- 1
  for (j in 2:84) { # Start from the second column
    for (i in 1:(j - 1)) { # Only fill below the diagonal
      loc_matrix[i, j] <- index
      index <- index + 1
    }
  }

  # using the loc_matrix for indexing, fill out the uvpm matrix
  for (i in 1:84) {
    for (j in 1:84) {
      if (i != j || loc_matrix[i,j] != 0) { # Skip the diagonal
        uvpm[i, j] <- mean(UVC[loc_matrix[i, j]])
      }
    }
  }

  #place the data in the upper tri into the lower tri
  for (i in 1:84) {
    for (j in 1:84) {
      if (i < j) {
        uvpm[j, i] <- uvpm[i, j]
      }
    }
  }

  #save the UVPM matrix
  if (!dir.exists(outputdir)) dir.create(outputdir, recursive = TRUE)
  setwd(outputdir)
  uvpm_file <- paste0('Average_', group, '_UVPM.rds')
  saveRDS(uvpm, uvpm_file)
  message(paste("UVPM for group", group, "has been recalculated and saved to", uvpm_file))
}


for (g in grouplist)
{
  if (g == 'fmci') limit = 4500:5500
  else limit = NA
  UVC2UVPM(group = g, limit = limit, modeldir = '~/Documents/Work/FinalFiles', outputdir = '~/Documents/Work/FinalFiles/fixed')
  
  
}
