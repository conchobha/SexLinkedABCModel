

#Load Functions
library(MASS)
library(readr) # Used for reading CSV file faster than the built in r function
library(dplyr) # used for selecting collums in the read functions
library(caTools) #Used for splitting data 
# Source Entire Folder function***
sourceEntireFolder <- function(folderName, verbose=FALSE, showWarnings=TRUE) { 
  # Takes in Name, Verbosity, and Warnings. Sources all files within a folder, one of which being the abc.R file. 
  #Frankly we only need to source the needed, but the model works so why change it
  files <- list.files(folderName, full.names=TRUE) # Lists all files within the folder
  files <- files[ grepl("\\.[rR]$", files) ]# Grab only R files 
  if (!length(files) && showWarnings) # Activates if there are no R files in the specified folder 
    warning("No R files in ", folderName)
  
  for (f in files) { #Runs through every R file in the Folder
    if (verbose) 
      cat("sourcing: ", f, "\n") 
    ## TODO:  add caught whether error or not and return that
    try(source(f, local=FALSE, echo=FALSE), silent=!verbose) 
    # Attempts to Read and execute all R files in the folder, giving us those variables that exist  
  }
  return(invisible(NULL))
}



replace_nan_with_NA <- function(list_of_matrices) { #Model needs NA values instead of nan in which our dataset gives, this fixes the data
  lapply(list_of_matrices, function(matrix) {
    # Convert the matrix to numeric, which will turn any non-numeric values to NA (including 'NaN' strings)
    numeric_matrix <- suppressWarnings(as.numeric(matrix))
    
    # Reshape the numeric vector back to the original matrix dimensions
    numeric_matrix <- matrix(numeric_matrix, nrow = nrow(matrix), ncol = ncol(matrix))
    
    # Replace NA values (originally 'NaN' strings) with 0
    numeric_matrix[is.na(numeric_matrix)] <- NA
    
    return(numeric_matrix)
  })
}
RunModel <- function(dataname='OD',
                     dn=1,sn,fn=100,num=42,burns=1000,
                     group='combined',form='mean', 
                     outputDIR="/N/u/conlcorn/BigRed200/SexLinkedProject/output",
                     dataloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData")
{ 
# Main function that is called once data is processed
#' @param dataname Defines the metric used 
#' @param dn Defines the dimension used for the ABC model
#' @param sn Defines the number of MCMC itterations
#' @param num Defines the seed, and can also be used as an itteration marker
#' @param burns The number of itterations that are discarded
#' @param group The group used in the definition 
#' @param form Defines if it is the mean or median data 
#' @param outputDIR Declares where the final model will be saved 
#' @param dataloc defines where the preprocessed data is stored 
  cpath <- getwd()
  start <- Sys.time()
	set.seed(dn*fn*num)
  setwd(cpath) # Set the working directory to the location of the script
  sourceEntireFolder("code_behavior_updated",verbose = FALSE, showWarnings=TRUE) # Source the entire folder containing the functions
  # This is where the ABC function is located

#Read the data 
  features <- paste0(group,"_",dataname,"_",form,".rds")
  labels <- paste0(group,"label.rds")
  setwd(paste0(dataloc,"/Features"))
  X <-readRDS(features)
  X_list <- X[[1]]
  setwd(paste0(dataloc,"/Labels"))
  labeldata <- readRDS(labels) #Pull the 
# We need to clean the data, as the model does not like NaN values
  labeldata_cleaned <- replace_nan_with_NA(labeldata)
  rm(X,labeldata) # Cleaning variables no longer needed
# Split for training and Testing sets 
#Get Validation, Test and Train ID
  train_ratio <- 0.8
# We are having an issue with train_ratio for some reason, only when ran using slurm
  n <- length(X_list)
# Generate training indices by sampling.
  train_indices <- sample(seq_len(n), size = floor(0.8 * n))

# Get the test indices as the ones not in train_indices.
  test_indices <- setdiff(seq_len(n), train_indices)

  X_train <- X_list[train_indices]
  X_test <- X_list[test_indices]
  Y_train <- labeldata_cleaned[train_indices]
  Y_test <- labeldata_cleaned[test_indices]

# Run ABC Model
  md=abc(X=X_train, Y=Y_train,W=NULL, H=NULL, K=dn,
       indices = NULL, indices_irt = NULL,
       seed = num, nscan = sn, burn = burns,odens = 10,
       print = TRUE, gof=FALSE, plot=FALSE,
       prior=list())
  end <- Sys.time()
  comptime <- as.numeric(difftime(end,start,units = "secs"))
# Return Model as file for future validation
  res=list("model"=md,"dn"=dn,"sn"=sn, "fn"=fn,"num"=num 
         ,'testX' = X_test,'testY' = Y_test, 'Comptime' = as.numeric(difftime(end,start,units = "secs"))
         )
  if (!dir.exists(outputDIR)) dir.create(outputDIR,recursive = TRUE) #should have been made before running the script, but just to be sure
  setwd(outputDIR)  # When running, we want to change this to a DIR we can access
  saveRDS(res,paste0("ADNI_",dataname,"_",group,"_",form,"_",sn,"_",burns,".rdata")) # Saves the Model as an R Object in the Above DIR
}
