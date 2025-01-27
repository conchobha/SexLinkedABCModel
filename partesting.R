#Load Functions
library(MASS)
library(foreach)
library(doParallel)
library(parallel)
library(readr) # Used for reading CSV file faster than the built in r function
library(dplyr) # used for selecting collums in the read functions
library(caTools) #Used for splitting data 
library(doParallel)
library(foreach)
library(MASS)
library(readr)
library(dplyr)
# Source Entire Folder function***
sourceEntireFolder <- function(folderName, verbose=FALSE, showWarnings=TRUE) { # Takes in Name, Verbosity, and Warnings
  files <- list.files(folderName, full.names=TRUE) # Lists all files within the folder
  files <- files[ grepl("\\.[rR]$", files) ]# Grab only R files 
  if (!length(files) && showWarnings) # Activates if there are no R files in the specified folder 
    warning("No R files in ", folderName)
  
  for (f in files) { #Runs through every R file in the Folder
    if (verbose) 
      cat("sourcing: ", f, "\n") 
    ## TODO:  add caught whether error or not and return that
    try(source(f, local=FALSE, echo=FALSE), silent=!verbose) # Attempts to Read and execute all R files in the folder, giving us those variables that exist  
  }
  return(invisible(NULL))
}
replace_nan_with_NA <- function(list_of_matrices) {
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
RunModel <- function(dataname,dn,sn,fn,num,burns,group,form)
{
  start <- Sys.time()
  node_name <- Sys.info()["nodename"]  # Get the node name
  print(paste0("Starting model on node: ", node_name, 
               " | Group: ", group, 
               " | Datametric: ", dataname, 
               " | Dimensionality: ", dn, 
               " | Replicate: ", num))
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/methods")
  sourceEntireFolder("code_behavior_updated",verbose = FALSE, showWarnings=TRUE)
  set.seed(num)
  # Read the data
  if(dataname == 'meanlength' || dataname == 'numoffibers') features <- paste0(group,"_",dataname,".rds")
  else features <- paste0(group,"_",dataname,"_",form,".rds")
  labels <- paste0(group,"label.rds")
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features")
  X <-readRDS(features)
  X_list <- X[[1]]
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Labels")
  labeldata <- readRDS(labels)
  # Normalize the X
  labeldata_cleaned <- replace_nan_with_NA(labeldata)
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
  print("Begining Fitting Model") #Used for tracking the step
  # Run ABC Model
  md=abc(X=X_train, Y=Y_train,W=NULL, H=NULL, K=dn,
         indices = NULL, indices_irt = NULL,
         seed = num, nscan = sn, burn = burns,odens = 10,
         print = TRUE, gof=FALSE, plot=FALSE,
         prior=list())
  
  
  # If you name the model variable df, it fails
  
  #corr <- cor(c(md$EFLPM),c(X_test))
  end <- Sys.time()
  comptime <- as.numeric(difftime(end,start,units = "secs"))
  # Return Model as file for future validation
  res=list("model"=md,"dn"=dn,"sn"=sn, "fn"=fn,"num"=num 
           ,'testX' = X_test,'testY' = Y_test, 'Comptime' = comptime
           #,"validation.id"=validation.id, "test.id"=test.id,"train.id"=train.id #-> Removed for now, will return later 
  )
  
 
  outfile <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/model_outputs/Dimensonality_testing/DM-",dn,"/Rep_",num)
  if (!dir.exists(outfile)) dir.create(outfile,recursive = TRUE)
  
  print("Finished Fitting, Saving Model")
  
  
  setwd(outfile)  # When running, we want to change this to a DIR we can access
  saveRDS(res,paste0("ADNI_",dataname,"_",group,"_",form,"_",sn,"_",burns,".rdata")) # Saves the Model as an R Object in the Above DIR
}

# For testing
dn <- 1:8 
sn <- 50000
fn <- 100
burns <- 1000
grouplist <- list('fcn', 'fmci', 'mcn', 'mmci', 'combined') # Ensure correct names for groups
datametriclist <- list(
  'FA', 'OD', 'meanlength', 'numoffibers' # Simplified for testing
)

# Load required libraries
library(doParallel)
library(foreach)
library(MASS)
library(readr)
library(dplyr)

# Load the Arguments for the Job Array 
args = commandArgs(trailingOnly = TRUE)
rep <-as.numeric(args[1])
dim <-as.numeric(args[2])

# Initialize and register the cluster
# Start the parallel backend with the desired number of cores
message("Intitalized Cluster")         

for (group in grouplist) {
  tryCatch({
    RunModel(datametriclist[4], dim, sn, fn, rep, burns, group, 'mean')
  }, error = function(e) {
    message(sprintf("An error occurred for group '%s': %s", group, e$message))
    # Continue to the next iteration
  })
}


# Load required library
