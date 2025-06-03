#Used to Run the Final Model after tweaking
# Source the Needed Functions
setwd(getwd())
source("dMRIABC.R") # Adds the Functions we need in order to run the model
#Set up The variables

fn <- 100
burns <- 1000
sn <- 200000    
grouplist <- list('fcn','fmci','fscd','mcn','mmci','mscd') # Ensure correct names for groups
#metriclist <- list("Da","Dr","FA","ICVF","MD","meanlength")
metriclist <- list("Da","Dr","MD") # Temp while I re run failed models 
#args = commandArgs(trailingOnly = TRUE)
#We are going to be passing the group and replication, group

#group <- grouplist[as.numeric(args[1])] # 6 groups
#dm <- metriclist[as.numeric(args[2])] #6

#Temp assignments
#dm <- metriclist[as.numeric(args[2])]
#grouptemp <- as.numeric(args[1])

group <- 'mmci'
dm <- 'Da'
#if (grouptemp == 2) {
 # group <- 'mmci'
#} else if (dm == 'Da' && grouptemp == 1) {
 # group <- 'fmci'
#} else {
 # group <- 'mcn'
#}

# Pull Which dim we are using 
#for OD
#dm 1 for FCN, 2 for others 

#dm <- 'OD'
rep <- 10
DM <- 1
output_loc <- paste0("/N/slate/conlcorn/SexLinkedProject/FinalModels/SuppModels/",dm,"/")
if (!dir.exists(output_loc)) dir.create(output_loc,recursive = TRUE)

print(paste0("Running with ",group, " - ", dm))
        # Starting rep
max_retries <- 10


for (attempt in 1:max_retries) {
  cat("Trying with rep =", rep, "and seed =", sn, "\n")
  
  result <- tryCatch({
    RunModel(
      dataname = dm,
      dn = DM,
      num = rep,
      sn = sn,
      group = group,
      outputDIR = output_loc
    )
    return(TRUE)  # Success
  }, error = function(e) {
    message(paste("Error with rep", rep, ":", e$message))
    rep <<- rep + 1  # Increment rep on error
    return(FALSE)
  })
  
  if (result) break  # Exit loop on success
}

if (!result) {
  cat("All attempts failed.\n")
}


