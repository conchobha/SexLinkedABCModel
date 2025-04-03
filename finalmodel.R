#Used to Run the Final Model after tweaking
# Source the Needed Functions
setwd(getwd())
source("dMRIABC.R") # Adds the Functions we need in order to run the model
#Set up The variables

fn <- 100
burns <- 1000
grouplist <- list('fcn','fmci',
                  'fscd','mcn','mmci','mscd'
                  ) # Ensure correct names for groups

args = commandArgs(trailingOnly = TRUE)
#We are going to be passing the group and replication, group
gr <- as.numeric(args[1])
if (gr != 7){ 
  group <- grouplist[gr]
  flag <- FALSE
}else{
    flag <- TRUE
    }# 6 groups
DM <- as.numeric(args[2])
rep <- as.numeric(args[3])
# Pull Which dim we are using 
#for OD
#dm 1 for FCN, 2 for others 

Metric <- 'meanlength'

output_loc <- paste0("/N/slate/conlcorn/SexLinkedProject/DimTesting/",Metric,"/",DM,"/Rep-",rep,"/")
if (!dir.exists(output_loc)) dir.create(output_loc,recursive = TRUE)


set.seed(rep*DM*as.numeric(args[1]))
if(!flag){ 
  RunModel(dataname=Metric,
         dn=DM,
         num=rep,
         sn=50000,
         group=group,
         outputDIR=output_loc
         )
}else{
  for (g in grouplist){
    
    RunModel(dataname=Metric,
             dn=DM,
             num=rep,
             sn=50000,
             group=g,
             outputDIR=output_loc
    )
    
    
  } 
  
  
}

