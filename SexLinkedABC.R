
#Load Functions
library(MASS)
library(readr) # Used for reading CSV file faster than the built in r function
library(dplyr) # used for selecting collums in the read functions
library(caTools) #Used for splitting data 
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


#This function parses through a folder, and for each file that has a matching entry in the IDs list
# will add that data to the output list
# IDs: A list of every person whose lable matches which group we want to test 
# type: Allows us to only pull the data for the model we want to run
#       final file  Will be in the structure of [ID]_SC_[type]_[form].csv
#form: Will allow us to pick if we want the mean, median, or normalized values 
# output a list that will hold our output
pull_data <- function(IDs, type, form, output) {#TODO Test function
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/SC_quant")
  counter <- 1
  
  for (i in IDs) {
    filename <- paste0(i$`Subject ID`, "_SC_", type, "_", form, ".csv")
    
    if (file.exists(filename)) {
      IDdata <- read_delim(filename, delim = ' ')
      
      # Store data as a list containing the ID and the matrix
      output[[counter]] <- list(ID = i$`Subject ID`, Data = IDdata)
      counter <- counter + 1
    } else {
      warning("Error: File: ", filename, " does not exist")
    }
  }
  
  return(output)
}


# Sets the DIR to the code needed of sourcing
setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/methods")
sourceEntireFolder("code_behavior_updated",verbose = FALSE, showWarnings=TRUE)

#Input Perameters for Execution

args <- commandArgs(trailingOnly = TRUE)
#We need ones to define
dataname <-strtoi(args[0]) #Which data are we  fitting model to, will be a string 
# Can be , Da, Dr, MD,FA,ICVF,OD,meanlength,numoffibers
#Number of Dimensions
dn <- strtoi(args[1])
#Number of MCMC Interations, needs to be large
sn <- strtoi(args[2])
#fn
fn <- strtoi(args[3])
#the Seed
num <- strtoi(args[4]) 
burns <-strtoi(args[5]) #Defines how many iterations we will remove in our model
group <- strtoi(args[6]) # Defines which of 7 groups we will fit our model to 
form <- strtoi(args[7]) #Defines if we are using mean, median, or normalized.

set.seed(num)


# In the future we may want to directly access the data storage, but this will be here until 
# We know the code wont mess with the source data 
# Load Data into X and Y
setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/") #Sets working DIR to location of data

# Remember each person falls into one of six places
# MCI, Normal cognition, or AD
# Male or Female

# First Step, get a list of which subject IDs fall within each group
# This is needed so we can choose which groups we will have
# We will have seven groups, the above groups plus a combined group
# sex can be either 2 or 1, the paper states that 2 is female
if(dataname != 'meanlength' || dataname != 'numoffibers'){ # We only have labels for the others
file_name <- paste0("iADRCautoPtxth_0.005_quant_", dataname, ".csv")
} else { # in case we are doing one of those, unsure at this moment
  #TODO Find out what to do in the case of these two
  stop("No Labels File exists currently")
} # Structure is: ID,Sex,diag,Tau,Amyloid
# Read the data
data <- read_delim(file_name, delim = ',', col_select = c(2,3,6,22,23)) # Change the delimiter to space

# Define empty lists for each group
fcnIDs <- list()
fmciIDs <- list()
fadIDs <- list()
mcnIDs <- list()
mmciIDs <- list()
madIDs <- list()
combinedIDs <- list()

# Loop to sort people into respective groups
for (i in 1:nrow(data)) {
  person <- data[i, ]
  
  # Check if the data is valid
  if (is.na(person$TauPos) || is.na(person$AmyPos)) { #Might need to change, since some have very few non-nan values
    next # Skip to the next iteration
  }
  
  # Append to combined list
  combinedIDs <- append(combinedIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
  
  # Sort based on SEX and diag
  if (person$SEX == 2) {  # Female
    if (person$diag == 'CN' || person$diag == 'CN (IO)') {
      fcnIDs <- append(fcnIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
    } else if (person$diag == 'MCI') {
      fmciIDs <- append(fmciIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
    } else {
      fadIDs <- append(fadIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
    }
  } else {  # Male
    if (person$diag == 'CN' || person$diag == 'CN (IO)') {
      mcnIDs <- append(mcnIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
    } else if (person$diag == 'MCI') {
      mmciIDs <- append(mmciIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
    } else {
      madIDs <- append(madIDs, list(list(`Subject ID` = person$`Subject ID`, SEX = person$SEX, diag = person$diag, TauPos = person$TauPos, AmyPos = person$AmyPos)))
    }
  }
}


# Check which group of IDs to pull from 
# Create a named list to map groups to their respective IDs
id_list <- list(
  fcn = fcnIDs,
  fmci = fmciIDs,
  fad = fadIDs,
  mcn = mcnIDs,
  mmci = mmciIDs,
  mad = madIDs,
  combined = combinedIDs
)

# Use the group variable to extract the corresponding IDs
IDs <- id_list[[group]]

# Check if the group is valid
if (is.null(IDs)) {
  stop("Invalid group specified.")
}

raw_input_data <- list()
raw_input_data <- pull_data(IDs, dataname,'mean',raw_input_data) 
# IDs: The list of IDs
#type
#form
# Assuming no errors in our pulling of the data, each index of the data array will have it's tau and amyloid
# as the same index in the IDs array (the 4th and 5 column)

#TODO Scale data for a mean of 0 and a sd of 1  

# Separate training and test sets 

set.seed(num) # For predictability, the final version might randomly set the seed 
# Current structure of our data, Two lists
#Raw_input_data: a list with the first column being the IDs, and the second being a VxV matrix 
#IDs: Holds our labels to match with our Data
#    Structure is: ID,Sex,diag,Tau,Amyloid 
# Take the Combined data, and split them into test and training sets

split <-sample.split(IDs, SplitRatio = 0.7)



#Pull the Train data, Train labels, and their testing counterparts
X_train <- raw_input_data[split]
Y_train <- IDs[split]
X_test <- raw_input_data[!split]
Y_test <- IDs[!split]

#Remove IDs from the data, since it is no longer needed to match the data 

#X should only be the data
#Y should only be the tau and Amy values 

# Run ABC Model
df=abc(X=X_train, Y=Y_train,W=NULL, H=NULL, K=dn,
       indices = NULL, indices_irt = NULL,
       seed = num, nscan = sn, burn = burns, odens = 10,
       print = FALSE, gof=FALSE, plot=FALSE,
       prior=list())

# Return Model as file for future validation
res=list("model"=df,"dn"=dn,"sn"=sn, "fn"=fn,"num"=num, "validation.id"=validation.id, "test.id"=test.id,"train.id"=train.id)



setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/output")  # When running, we want to change this to a DIR we can access
saveRDS(res,paste0("ADNI",dn,sn,fn,"_",num,".rdata")) # Saves the Model as an R Object in the Above DIR