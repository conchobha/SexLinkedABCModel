# This script is what is called to define the parameters, and run the model

setwd(getwd()) #All needed functions should be in the same directory as this script
source("dMRIABC.R") # Source our model script

#Set up The variables
fn <- 100
burns <- 1000 #how many itterations to burn
itterations <- 50000 #how many itterations to run
grouplist <- list('fcn','fmci','fscd','mcn','mmci','mscd') # Define the list of our groups

args = commandArgs(trailingOnly = TRUE) #pull the arguments from the command line call
#We are going to be passing the group and replication, group
gr <- as.numeric(args[1]) #Define which group we are using out of the six total
group <- grouplist[gr]


DM <- as.numeric(args[2]) #Define which dim we are using out of the 1-8
rep <- as.numeric(args[3]) #Define which replication we doing. Impacts both seed and output location



Metric <- 'OD' 

output_loc <- paste0("/N/slate/conlcorn/SexLinkedProject/DimTesting/",Metric,"/",DM,"/Rep-",rep,"/") #make sure to make this directory before running the script to avoid errors
if (!dir.exists(output_loc)) dir.create(output_loc,recursive = TRUE)


set.seed(rep*DM*as.numeric(args[1])) #sets the seed in a way every run is unique, but repeatable

RunModel(dataname=Metric, #sourced from dMRIABC.R
         dn=DM,
         num=rep,
         sn=itterations,
         group=group,
         outputDIR=output_loc
)
