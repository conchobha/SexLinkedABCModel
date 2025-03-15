#Used to Run the Final Model after tweaking
# Source the Needed Functions
setwd(getwd())
source("dMRIABC.R") # Adds the Functions we need in order to run the model
#Set up The variables

fn <- 100
burns <- 1000
grouplist <- list('fcn','fmci','fscd','mcn','mmci','mscd') # Ensure correct names for groups

args = commandArgs(trailingOnly = TRUE)
#We are going to be passing the group and replication, group
group <- grouplist[1] # 6 groups
DM <- as.numeric(args[1])
rep <- as.numeric(args[2])
# Pull Which dim we are using 
#for OD
#dm 1 for FCN, 2 for others 

dm <- 'OD'

output_loc <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/DimTesting/",DM,"/Rep-",rep,"/")
if (!dir.exists(output_loc)) dir.create(output_loc,recursive = TRUE)
set.seed(rep)

RunModel(dataname=dm,
         dn=DM,
         num=rep,
         sn=50000,
         group=group,
         outputDIR=output_loc
         )

