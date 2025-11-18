



MakeSigEdge <- function(g1, g2 , metric = "OD", av = FALSE, range = NA,
                    atlasloc = '~/Documents/Work/FinalFiles/finalAtlas.rds',
                    modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                    outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps_Final/", within = TRUE
                    
) {


 
  setwd(modeldir)
  o <- outputdir 
  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
 
  if (!is.na(g2)) {
    fname <- paste0(g1, "_vs_", g2, "_", metric, "_SigEdges")
    # Modify the code to calc the 95% CI for these. If g2 is much less, set it to -1, elif it's more, 1, else 0
    CI <- function(x) {
      quantile(x, probs = c(0.025, 0.975))
    }
    # load the MCMC data for both groups
    if(av){
      model1name <- paste0("Average_", g1, "_UVC.rds")
      model2name <- paste0("Average_", g2, "_UVC.rds")
      UVC1_1 <- readRDS(model1name)
      UVC2_1 <- readRDS(model2name)
    }else{
      g1name <- paste0("ADNI_",metric,"_", g1, "_mean_5e+05_10000.rdata")
      data <- readRDS(g1name)
      model <- data$model
      UVC1_1 <- model$UVC
      
      g2name <- paste0("ADNI_",metric,"_", g2, "_mean_5e+05_10000.rdata")
      data <- readRDS(g2name)
      model <- data$model
      UVC2_1<- model$UVC
    }
    if(g2 == 'fmci') UVC2 <- UVC2_1
    else UVC2 <- UVC2_1
    if(g1 == 'fmci') UVC1 <- UVC1_1
    else UVC1 <- UVC1_1
    #UVC1 <- UVC1_1
    #UVC2 <- UVC2_1
    
    
    hi.g1 <- t(apply(UVC1, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) # returns a 84x83/2 matrix, which is the CI's for each unique connection
    hi.g2 <- t(apply(UVC2, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
    
    # Create an 84x84 matrix filled with zeros
    loc_matrix <- matrix(0, nrow = 84, ncol = 84)
    
    # Start filling the upper triangle downwards, column by column
    index <- 1
    for (j in 2:84) { # Start from the second column
      for (i in 1:(j - 1)) { # Only fill below the diagonal
        loc_matrix[i, j] <- index
        index <- index + 1
      }
    }
    
    # We can use this matrix to pull the proper CI's for each connection, then return the matrix to be made into the heatmap
    # Use the number stored in the loc_matrix to pull the CI's from the hi.g1 and hi.g2 matrices
    matrix_data <- matrix(0, nrow = 84, ncol = 84)
    for (i in 1:84) {
      for (j in 1:84) {
        n <- loc_matrix[i, j]
        if (n == 0) next
        
        tmp1 <- hi.g1[n, ]
        tmp2 <- hi.g2[n, ]
        if (max(tmp1[1], tmp2[1]) > min(tmp1[2], tmp2[2])) { #checks if they don't overlap
          if (tmp1[2] <= tmp2[1]) { #checks if g2 is larger or smaller than g1
            matrix_data[i, j] <- 1
          } else if (tmp1[1] >= tmp2[2]) {
            matrix_data[i, j] <- -1
          } 
        }
      }
    }
    # We need to make sure the matrix is symmetric
    for (i in 1:nrow(matrix_data)) {
      for (j in 1:ncol(matrix_data)) {
        if (i != j) {
          # Check if matrix_data[i, j] is zero but matrix_data[j, i] is non-zero
          if (matrix_data[i, j] == 0 && matrix_data[j, i] != 0) {
            matrix_data[i, j] <- matrix_data[j, i]
          }
          # Check if matrix_data[j, i] is zero but matrix_data[i, j] is non-zero
          else if (matrix_data[j, i] == 0 && matrix_data[i, j] != 0) {
            matrix_data[j, i] <- matrix_data[i, j]
          }
        }
      }
    }
  } else matrix_data <- center_matrix(g1data) # if we are only looking at the heatmap of one group
  
  # Remove rows 35-49 and 84 from matrix_data
  rows_to_remove <- c(35:49, 84)
  matrix_data <- matrix_data[-rows_to_remove, -rows_to_remove]
  
  fmat <- matrix_data
  fmat[matrix_data==1] <- 2
  fmat[matrix_data == -1] <- 1
 # Save matrix_data as a .edge file
  write.table(fmat, file = paste0(o, fname, ".edge"), row.names = FALSE, col.names = FALSE)

  #Outside of saving the 'Total' significant edges, we also want to save the significant edges for each lobe. We can do this by going row by row, and setting the numbers to 0 if they aren't in the lobe we want to look at, and doing this for all lobes
  #pull the atlas, to see which regions correspond to which lobes
  load(atlasloc)


  atlas <- final_df[-rows_to_remove, -rows_to_remove]
  #going by the $LOBE column of the atlas, we can see which regions correspond to which lobes

  lobelist <- unique(atlas$LOBE)

  for(l in lobelist){
    lobe_matrix <- matrix_data
    regions <- atlas$row_number[which(atlas$LOBE == l)]
    # Remove rows and columns that aren't in regions
    all_regions <- 1:nrow(matrix_data)
    regions_to_remove <- setdiff(all_regions, regions)
  
    if(within) lobe_matrix[regions_to_remove,regions_to_remove] <- 0 #only if we want just within lobe connections
    else  
      lobe_matrix[regions_to_remove, ] <- 0
    # Save the lobe_matrix as a .edge file
    
    mat_new <- lobe_matrix
    
    # Only change original 1s
    mat_new[lobe_matrix == 1] <- 2
    
    # Only change original -1s
    mat_new[lobe_matrix == -1] <- 1
    
    write.table(mat_new, file = paste0(o, fname, "_", l, ".edge"), row.names = FALSE, col.names = FALSE)
  }
}

#run a test of this using the 500k model
MakeSigEdge(g1 = 'mcn', g2 = 'mmci', metric = 'OD', av = FALSE,
            modeldir = '~/Documents/Work/FinalFiles/500k',
            outputdir = '/home/connor/Documents/Work/plots/BrainPlots1/')



g1 <- 'mcn'
g2 <- 'mmci'
metric <- 'OD'
av <- FALSE
modeldir <- '~/Documents/Work/FinalFiles/500k'
outputdir <- '/home/connor/Documents/Work/plots/BrainPlots1/'


