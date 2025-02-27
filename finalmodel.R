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
group <- grouplist[as.numeric(args[1])] # 6 groups
dm <- 'OD'
rep <- 42
# Pull Which dim we are using 
#for OD
#dm 1 for FCN, 2 for others 
d <- 2
if(group == "fcn"){ d <- 1}


output_loc <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/",dm,"/Rep-",rep,"/")
if (!dir.exists(output_loc)) dir.create(output_loc,recursive = TRUE)

RunModel(dataname=dm,
         dn=d,
         num=rep,
         sn=200000,
         group=group,
         outputDIR=output_loc
         )

