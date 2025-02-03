library(coda)
library(MASS)
library(boot)
library(latentnet)
library(lvm4net)
library(amen)
library(pROC)
library(corrplot)
Heatmap <- function( metric, group, value = TRUE){
  modelname=paste0("ADNI_",metric,"_",group,"_mean_1e+05_1000.rdata")
  #Takes in a Model, and Saves a PDF holding an image of the heat map of the estimated correlation 
  location <- paste0('/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/',metric,'/Rep-1/')
  setwd(location)
  data <- readRDS(modelname)
  model <-data$model
  matrix_data <- model$UVPM
  labels <- seq(1, 84, by=8)
  o <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps/",metric)
  if (!dir.exists(o)) dir.create(o,recursive = TRUE)
  setwd(o)
  
  pdf(file = paste0(metric,group,'_heatmap.pdf'))
  # Plot with custom labels
  if(value){
  corrplot(matrix_data,
           method = "color",                            # Use color method for heatmap
           col = colorRampPalette(c("blue", "white", "red"))(200), # Color gradient
           tl.pos = "lt",                               # Place axis labels on the left and top
           tl.col = "black",                            # Set axis numbers color to black
           tl.cex = 0.3,                                # Adjust the font size of the axis numbers
           addgrid.col = "white",                       # Set grid lines to white
           number.cex = 0.1,
           cl.pos = "b",                                # Place color legend at the bottom
           is.corr = FALSE                              # Indicate this is not a correlation matrix
  )
  }else{
    corrplot(matrix_data,
             method = "color",                            # Use color method for heatmap
             col = c("orange", "cyan"),                   # Use two colors: orange for negative, cyan for positive
             breaks = c(-Inf, 0, Inf),                    # Define breaks to apply colors
             tl.pos = "lt",                               # Place axis labels on the left and top
             tl.col = "black",                            # Set axis numbers color to black
             tl.cex = 0.3,                                # Adjust the font size of the axis numbers
             addgrid.col = "white",                       # Set grid lines to white
             number.cex = 0.1,
             cl.pos = "b",                                # Place color legend at the bottom
             is.corr = FALSE                              # Indicate this is not a correlation matrix
    )
    
  }
  # Add custom stepped labels    
  dev.off()
}
Traceplot <- function(modelname="ADNI_Da_combined_mean_10000_1000.rdata",group,filedir){
  setwd(filedir)
  data <- readRDS(modelname)
  model1=data$model
  UVPM1<-model1$UVPM
  UVC1 <- model1$UVC
  # model1$TAC
  l <- data$sn/10
  mc1 <- mcmc(data=UVC1[,90:100],start = 1, end = l, thin=1)
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/Traceplots")
  pdf(file = paste0(group,'_Traceplot.pdf'))
  traceplot(mc1,ylab = 'Covariance Estimate')
  dev.off()
}
Correlation <- function (modelname="ADNI_Dr_combined_mean_10000_1000.rdata",modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output" ){
  
  setwd(modelloc)
  data <- readRDS(modelname)
  model1=data$model
  x_test <- data$testX
  l <- length(x_test)
  est <- model1$EFlPM
  est_split <- est[1:l]
  vec1 <- unlist(x_test)
  vec2 <- unlist(est_split)
  correlations <- cor(vec1, vec2)
  # View the correlation
  print(paste0("Correlation for ",modelname,": ",  correlations))
  return(correlations)
  
}
histo <- function(modelname,metric, 
                  modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output" , 
                  outputloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots" ){# makes a histogram for the specified model
  setwd(modelloc)
  data <- readRDS(modelname)
  model1=data$model
  dt <- unlist(model1$EFlPM)
  setwd(outputloc)
  pdf(file = paste0(metric,'_histogram.pdf'))
  hist(dt,
       main=paste0("Frequency for values for the EFLPM: ",metric),
       col="darkmagenta"
  )
    
  dev.off()
    
  
}
dimtest<-function(dm){
  start <- Sys.time()
  #Saves a PDF file, that can display the Correclation for each group, accross different
  #Dimensionality of a model. 
  #Can Average accross multiple replications as well 
  output <- matrix(data = NA, nrow = 8, ncol = 5)
  # Defining our 3d Output Matrix
  # These will define our DM, Row, and Column
  toutput <-  array(dim = c(8, 5, 10)) # First index is Dim, Next is group, last is rep
  
  # We need to load each model, get the correlation, and the store it in the right spot
  dim_list <- list(1,2,3,4,5,6,7,8)
  group_list <- list('fcn','fmci', 'mcn','mmci','combined')
  for (i in 1:10){ # loops over each Replication
    for(d in dim_list){
      location <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/model_outputs/Dimensonality_testing/DM-",d,"/Rep_",i)
      setwd(location)
      count <-1
      for(g in group_list){
        #load the file
        filename <- paste0("ADNI_",dm,"_",g,"_mean_50000_1000.rdata")
        
        data <-tryCatch({
          readRDS(filename)
        },
        error = function(e){
          print(paste0("Error reading file: ",filename," for DIM" ,d, " - Setting to NULL"))
          return(NULL)
        })
        if(!is.null(data)){
          model1=data$model
          x_test <- data$testX
          l <- length(x_test)
          est <- model1$EFlPM
          est_split <- est[1:l]
          vec1 <- unlist(x_test)
          vec2 <- unlist(est_split)
          correlations <- cor(vec1, vec2)
          #get the correlation
          
          #store it in the matrix 
          toutput[d,count,i] <- correlations
        }else output[d,count] <- NA
        count <- count + 1
      }
    }
  }
  #Average together all of the last dim together into the output matrix 
  #Average together all of the last dim together into the output matrix 
  for(g in 1:5){ # For each Group set
    for(d in 1:8){ # For each Dim Set
      temp <- numeric(10) # Initialize a numeric vector with length 10
      for(i in 1:10){
        temp[i] <- toutput[d,g,i] # Store values directly as numeric
      }
      # Get the Average of these, excluding NA values
      tav <- mean(temp, na.rm = TRUE)
      output[d,g] <- tav
    }
  }
  
  
  #Call the function to build the table with this data 
  library(grid)
  library(gridExtra)
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/Dim")
  rownames(output) <- c("1","2", "3", "4", "5", "6", "7", "8")  # Dimensions
  colnames(output) <- c("Female Healthy", "Female MCI", "Combined", 
                        "Male Healthy", "Male MCI")
  name <- paste0("Dim_Testing_",dm,"_rep10.pdf")
  pdf(name, width = 12, height = 6)
  grid.newpage() 
  grid.text(paste0("Dim Testing for ", dm,". Averaged Across 10 Replications"), y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Convert the matrix to a table and display it
  grid.table(output)
  
  # Close the PDF device
  dev.off()
  end <- Sys.time()
  s <- as.numeric(difftime(end,start,units = "secs"))
  m <- as.numeric(difftime(end,start,units = "mins"))
  message(paste0("Execution Completed in ",s," Seconds/ ",m," Minutes"))
}
# Function to print values from .rdata files
print_values_from_rdata <- function(directory) {
  # Check if the directory exists
  if (!dir.exists(directory)) {
    stop("Directory does not exist: ", directory)
  }
  
  # Get a list of all .rdata files in the directory
  rds_files <- list.files(directory, pattern = "\\.rdata$", full.names = TRUE)
  
  if (length(rds_files) == 0) {
    message("No .rds files found in the directory: ", directory)
    return(NULL)
  }
  
  # Iterate over each .rds file
  for (file in rds_files) {
    
    # Read the RDS file
    data <- tryCatch(
      {
        readRDS(file)
      },
      error = function(e) {
        cat("Error reading file:", file, "\n", e$message, "\n")
        return(NULL)
      }
    )
    
    if (is.null(data)) {
      next
    }
    sn <- data$sn
    Comptime <- data$Comptime
    print(paste0(
      "For Itterations of ",sn,", It Took "
      ,round(Comptime,0)," Seconds/ "
      ,round(Comptime/60,2)," Minutes/ "
      ,round(Comptime/3600,2)," Hours/ "
      ,round(Comptime/3600/24,2)," Days"
    ))
    
  }
}
#Confidence intervals
# Prints a confidence interval for the given group, comparing the UVPM between the healthy and unhealthy participants in those groups 
CI <- function(group = 'f', metric = 'OD', modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/Rep-1/") {
  library(dplyr)
  setwd(modeldir) # Set the Directory to where we store the models
  cn_name <- paste0("ADNI_", metric, "_", group, "cn_mean_2e+05_1000.rdata") # file name for the cn group
  mci_name <- paste0("ADNI_", metric, "_", group, "mci_mean_2e+05_1000.rdata") # file name for mci group
  
  credible_interval <- function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }
  
  vec1.tmp <- array(rep(0, 84 * 84), dim = c(84, 84))
  vec1.tmp[c(20, 36, 39, 43, 48,62, 69, 74, 80), ] <- 100
  vv <- vec1.tmp[upper.tri(vec1.tmp, diag = FALSE)]
  index <- which(vv != 0)
  
  # Load CN data
  df <- readRDS(cn_name)
  model1 <- df$model
  UVC1 <- model1$UVC
  hi.g1 <- t(apply(UVC1, 2, credible_interval))
  rhi.g1 <- hi.g1[index, ]
  UVPM <- model1$UVPM
  pm <- UVPM[upper.tri(UVPM, diag = FALSE)]
  pm1 <- pm[index]
  order_pm1 <- order(pm1)
  sorted_pm1 <- pm1[order_pm1]
  negpe.g1 <- sorted_pm1[sorted_pm1 < 0]
  pospe.g1 <- sorted_pm1[sorted_pm1 > 0]
  sorted_rhi.g1 <- rhi.g1[order_pm1, ]
  negatives.g1 <- sorted_rhi.g1[sorted_pm1 < 0, ]
  positives.g1 <- sorted_rhi.g1[sorted_pm1 > 0, ]
  order_pos <- order_pm1[sorted_pm1 > 0] # Define order_pos explicitly
  order_neg <- order_pm1[sorted_pm1 < 0] # Define order_neg explicitly
  
  # Load MCI data
  df <- readRDS(mci_name)
  model1 <- df$model
  UVC1 <- model1$UVC
  hi.g5 <- t(apply(UVC1, 2, credible_interval))
  rhi.g5 <- hi.g5[index, ]
  UVPM <- model1$UVPM
  pm <- UVPM[upper.tri(UVPM, diag = FALSE)]
  pm5 <- pm[index]
  sorted_pm5 <- pm5[order_pm1]
  negpe.g5 <- sorted_pm5[sorted_pm1 < 0]
  pospe.g5 <- sorted_pm5[sorted_pm1 > 0]
  sorted_rhi.g5 <- rhi.g5[order_pm1, ]
  negatives.g5 <- sorted_rhi.g5[sorted_pm1 < 0, ]
  positives.g5 <- sorted_rhi.g5[sorted_pm1 > 0, ]
  
  # Safeguard function for length consistency
  ensure_length_match <- function(x, y, target_length) {
    x <- head(x, target_length)
    y <- head(y, target_length)
    list(x = x, y = y)
  }
  
  # Build out the CIs
  plot_CI_95_dem <- function(size, pospe.g1, positives.g1, pospe.g5, positives.g5, order_pos, is_positive = TRUE) {
    offset <- 0.4 # Controls point and line offset from nominal value
    colsuse <- c('steelblue', 'orange2')
    B <- size
    
    for (j in 1:2) {
      if (j == 1) {
        x_1 <- positives.g1[, 2]
        x_0 <- positives.g1[, 1]
        x_vals <- pospe.g1
      } else {
        x_1 <- positives.g5[, 2]
        x_0 <- positives.g5[, 1]
        x_vals <- pospe.g5
      }
      
      # Ensure lengths match
      lengths <- ensure_length_match(x_vals, 1:B + offset * ifelse(j == 1, 1, -1), B)
      x_vals <- lengths$x
      y_vals <- lengths$y
      
      if (j == 1) {
        plot(x_vals, y_vals, pch = 20,
             #xlim = c(-0.4, 0.4),
             #ylim = c(0,B),
             xlab = "Connectivity",
             ylab = "Different positions",
             cex = 1,
             col = colsuse[j],
             yaxt = 'n')
      } else {
        points(x_vals, y_vals, pch = 20,
               col = colsuse[j],
               yaxt = 'n')
      }
      
      for (i in 1:B) {
        if (between(0, x_0[i], x_1[i])) {
          segments(x_0[i], y_vals[i], x_1[i], y_vals[i], col = colsuse[j], lwd = 1.5)
        } else {
          segments(x_0[i], y_vals[i], x_1[i], y_vals[i], col = colsuse[j], lwd = 1.5)
        }
      }
    }
    
    axis(side = 2, at = 1:B, labels = order_pos, cex.axis = 0.5, padj = 0.6, las = 2)
    axis(side = 1)
    legend("bottomright", legend = c('Healthy', 'Dementia'), lty = 1, col = colsuse, cex = 1)
  }
  
  pdf(paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/CIs/",group,metric,"CIPos.pdf"), height = 11, width = 7)
  message("Generating positive connectivity plot...")
  plot_CI_95_dem(min(length(pospe.g1), length(pospe.g5)), pospe.g1, positives.g1, pospe.g5, positives.g5, order_pos)
  dev.off()
  message("Positive connectivity plot saved.")
  
  ### Negative group
  pdf(paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/CIs/",group,metric,"CINeg.pdf"), height = 11, width = 7)
  message("Generating negative connectivity plot...")
  plot_CI_95_dem(min(length(negpe.g1), length(negpe.g5)), negpe.g1, negatives.g1, negpe.g5, negatives.g5, order_neg, is_positive = FALSE)
  dev.off()
  message("Negative connectivity plot saved.")
}

#Finds and Saves the ten regions in the Estimated Connectivity in which has the highest difference between two given models
find_difference <- function(res=0,model1='ADNI_OD_fmci_mean_2e+05_1000.rdata',model2='ADNI_OD_mmci_mean_2e+05_1000.rdata',modeldir='/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/Rep-1'){
  library(stringr)
  #Determine the Two groups we are comparing and Metric  
  matches <- str_match(model1, "^ADNI_([^_]+)_([^_]+)_")[,2:3]
  metric <-matches[1]
  group1 <-matches[2]
  matches <- str_match(model2, "^ADNI_([^_]+)_([^_]+)_")[,3]
  group2 <- matches
  #Load Model1 and Store the UVPM data 
  setwd(modeldir)
  
  data <- readRDS(model1)
  model <-data$model
  matrix_data1 <- model$UVPM
  
  data <- readRDS(model2)
  model <-data$model
  matrix_data2 <- model$UVPM
  
  
  #We need to find both the Positve And Negative Values 
  pos_result <- matrix_data1-matrix_data2
  neg_result <- matrix_data2-matrix_data1
  abs_result <- abs(matrix_data1 - matrix_data2)
  #find the top 10 combinations
  if (res == 1) result_matrix <- pos_result
  else if (res == -1) result_matrix <- neg_result
  else result_matrix <- abs_result
  top_10_values <- sort(as.vector(result_matrix), decreasing = TRUE)[1:10]
  
  # Print the top 10 values
  print(top_10_values)
  # Get the upper triangle of the matrix (excluding diagonal)
  upper_triangle <- upper.tri(result_matrix)
  
  # Extract values from the upper triangle
  upper_values <- result_matrix[upper_triangle]
  
  # Find the indices of the top 10 values in the upper triangle
  sorted_indices <- order(upper_values, decreasing = TRUE)[1:10]
  
  # Get the row and column indices of the top 10 values
  row_col_indices <- which(upper_triangle, arr.ind = TRUE)[sorted_indices, ]
  
  paired_t_test <- t.test(as.vector(matrix_data1), as.vector(matrix_data2), paired = TRUE)
  print(paired_t_test)
  # Combine the results into a data frame
  top_10_with_positions <- data.frame(
    Value = round(upper_values[sorted_indices],4),
    Row = row_col_indices[, 1],
    Column = row_col_indices[, 2]
  )
  print(top_10_with_positions)
  # Print the result
  setwd(paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps/",metric))
  pdf(file = paste0(group1,group2,'_heatmap.pdf'))
  corrplot(result_matrix,
           method = "color",                            # Use color method for heatmap
           col = colorRampPalette(c("white","orange","red","black"))(200), # Color gradient
           tl.pos = "lt",                               # Place axis labels on the left and top
           tl.col = "black",                            # Set axis numbers color to black
           tl.cex = 0.3,                                # Adjust the font size of the axis numbers
           addgrid.col = "white",                       # Set grid lines to white
           number.cex = 0.1,
           cl.pos = "b",                                # Place color legend at the bottom
           is.corr = FALSE  )                            # Indicate this is not a correlation matrix
  dev.off()
  
}

pull_connect <- function(group1,group2,edge1,edge2){
  #Pulls the connectivity data for a given edge, from each group, and is returned as a dataframe
  
  #Make the filenames for the two groups
  group1name <- paste0(group1,"_OD_mean.rds")
  group2name <- paste0(group2,"_OD_mean.rds")
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features/")
  
  #pull the data from both files
  
  data1 <- readRDS(group1name)
  data2 <- readRDS(group2name)
  list1 <- list()
  list2 <- list()
  # Add the data to the list
  for (person in data1[[1]])list1 <- append(list1,person[edge1,edge2])
  for (person in data2[[1]])list2 <- append(list2,person[edge1,edge2])
  listcomb <- list()
  listcomb[[1]] <- list1
  listcomb[[2]] <- list2
  #Create the dataframe
 
  return(listcomb)
}

find_smallest <- function(g1 = 'fcn',g2 = 'fmci'){
  smallest_p_values <- data.frame(p_value = numeric(0), i = integer(0), j = integer(0))
  for (i in 1:84) {
    for (j in 1:84) {
      d <- pull_connect(g1, g2, i, j)
      sublist1 <- d[[1]]  # List of length 13
      sublist2 <- d[[2]]  # List of length 32
      
      # Access individual vectors within the sublists
      sublist1 <- unlist(sublist1)
      sublist2 <- unlist(sublist2)
      
      # Perform a t-test
      t_test_result <- t.test(sublist1, sublist2)
      
      if (!is.na(t_test_result$p.value)) {
        # Add the current p-value and combination to the data frame
        smallest_p_values <- rbind(smallest_p_values, data.frame(p_value = t_test_result$p.value, i = findEdge(i), j = findEdge(j)))
        
        # Keep only the 5 smallest p-values
        smallest_p_values <- smallest_p_values[order(smallest_p_values$p_value), ]
        if (nrow(smallest_p_values) > 5) {
          smallest_p_values <- smallest_p_values[1:5, ]
        }
      }
    }
  }

  # Display the results
  cat("Top 5 smallest p-values and their combinations:\n")
  print(paste0("For ",g1," and ",g2,"\n"))
  print(smallest_p_values)
}


findEdge <- function(number, filename = "iADRC_Struture_Diffusion_Tau_Abeta_84ROIcombo.xlsx", dataloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/") {
  # Function to search for a given edge name, given its number 
  library(readxl)
  
  if (!(number %in% 1:84)) {
    warning("Number not in Range")
    return("Number not in Range")
  }
  
  # Open the file that stores the key 
  setwd(dataloc) # Sets Directory to the data directory 
  
  # Suppress messages and warnings while reading the Excel file
  suppressMessages({ # Currently slow, could improve by optionally giving it an already loaded file 
    suppressWarnings({
      df <- read_excel(filename, sheet = 2)
    })
  })
  
  # Go to the Index
  result <- df$ICV[df$`...5` == number] 
  
  # Return the String 
  rm(df)
  return(result[9]) # Files have a structure R does not like, this works fine however
}



findavSD <- function(group){
  #create the 84x84 matrix 
  SDs <- matrix(NA, nrow = 84, ncol = 84)
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features/")
  
  # Assuming `group` is defined
  
  
  for (i in 1:84) {
    for (j in 1:84) {
      groupname <- paste0(group, "_OD_mean.rds")
      data1 <- readRDS(groupname)
      
      # Initialize an empty numeric vector to store values
      list1 <- numeric()
      
      # Loop through people in data1[[1]] and extract person[i,j]
      for (person in data1[[1]]) {
        # Ensure that person[i,j] is numeric and valid
        val <- person[i, j]
        
        # Only add numeric values to list1, ignoring non-numeric or NA values
        if (!is.na(val) && is.numeric(val) && val != 0) {
          list1 <- c(list1, val)
        }
      }
      
      # Calculate the SD of the filtered list1
      if (length(list1) > 0) {
        s <- sd(list1, na.rm = TRUE)  # Use `na.rm = TRUE` to ignore NAs
        SDs[i, j] <- s
      } else {
        SDs[i, j] <- NA  # Assign NA if list1 is empty (no valid data)
      }
    }
  }
  # Flatten the SDs matrix into a vector, ensuring it's numeric
  flat_SDs <- as.vector(SDs)
  
  # Ensure it's numeric (in case there are any non-numeric values)
  flat_SDs <- as.numeric(flat_SDs)
  
  # Remove any NA or NaN values (including zero if you want)
  flat_SDs <- flat_SDs[!is.na(flat_SDs) & flat_SDs != 0]
  
  # Calculate the overall average SD, if the vector is not empty
  if(length(flat_SDs) > 0) {
    avg_SD <- mean(flat_SDs)
    print(paste0("Average SD for group: ",group))
    print(avg_SD)
  } else {
    print("No valid data to calculate the average SD.")
  }
}

grouplist <- list('fcn','fmci','fscd','mcn','mmci','mscd')

