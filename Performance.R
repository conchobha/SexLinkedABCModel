# This file is planned to be a collection of final versions of functions that are used to gather analysis
# Working versions of functions will be in the MeasurePerformance.R file
# This file should be sourced and used by other analysis scripts

# ~/Documents/Work/FinalFiles/500k/AveragedResults500k
#' @title Measure Performance
#' @description This function measures the performance of a model using various metrics. Most of these functions are the 
#'  final versions of the functions used to measure performance in the project. In progress files are in MeasurePerformance.R



#APM 

disgraph <- function(loc = '~/Downloads/APM') { # violin plot code 
  setwd(loc)
  grouplist <- c("fcn", "fmci", "mcn", "mmci")  # Change this to a vector
  datalist <- list()
  
  for (g in grouplist) {
    data <- readRDS(paste0('Average_', g, "_APM.rds"))
    APM <- data$APM
    matrix_data <- matrix(unlist(APM), nrow = length(APM), byrow = TRUE)
    Av <- colMeans(matrix_data)
    datalist[[g]] <- Av
  }
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Convert list to data frame
  datalist_df <- stack(datalist)
  colnames(datalist_df) <- c("value", "group")
  
  # Recode group names
  group_labels <- c(
    fcn = "Healthy Female",
    fmci = "Female MCI",
    mcn = "Healthy Male",
    mmci = "Male MCI"
  )
  
  # Count number of datapoints per group
  counts <- datalist_df %>%
    group_by(group) %>%
    summarise(n = n(), .groups = 'drop')
  
  # Create group labels with line break before n
  counts$label <- paste0(group_labels[counts$group], "\n(n = ", counts$n, ")")
  
  # Join updated labels into datalist_df
  datalist_df <- datalist_df %>%
    left_join(counts, by = "group")
  
  # Set factor levels in the correct order
  label_levels <- sapply(grouplist, function(g) {
    paste0(group_labels[g], "\n(n = ", counts$n[counts$group == g], ")")
  })
  
  datalist_df$label <- factor(datalist_df$label, levels = label_levels)
  
  # Violin plot with group-specific n values in x-axis labels
  ggplot(datalist_df, aes(x = label, y = value)) +
    geom_violin(fill = "gray85", color = "black", trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
    theme_minimal() +
    labs(y = "Participant Connectivity Density", x = "") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 30, angle = 45, hjust = 1),
      axis.title = element_text(size = 26),
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold")
    )
}

APMTesting <- function(g1 = 'fcn',g2 = 'fmci',loc = "/N/slate/conlcorn/SexLinkedProject/FinalModelStore")
{
  #This function will take two groups, pull the APM data from each group, and then run a paired t-test on them.

  g1name <- paste0('Average_',g1,"_APM.rds")
  g2name <- paste0('Average_',g2,"_APM.rds")
  setwd(loc)
  data1 <- readRDS(g1name)
  data2 <- readRDS(g2name)
  
  APM1 <- data1$APM
  APM2 <- data2$APM

  # Convert the APM data to a matrix
  matrix_data1 <- matrix(unlist(APM1), nrow = length(APM1), byrow = TRUE)
  matrix_data2 <- matrix(unlist(APM2), nrow = length(APM2), byrow = TRUE)

  #Calculate the Average APM for each participant
  RowAv1 <- colMeans(matrix_data1)
  RowAv2 <- colMeans(matrix_data2)

  #do the t-test
  t_test_result <- t.test(RowAv1, RowAv2, paired = TRUE)
  print(paste("For ", g1, " vs ", g2))
  print(t_test_result)  
}


Traceplot <- function(modelname = "ADNI_Da_combined_mean_10000_1000.rdata", 
                      group, av = TRUE, 
                      filedir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                      outputdir) {
  library(coda)
  setwd(filedir)
  if(!av){
  data <- readRDS(modelname)
  model1 <- data$model
  UVPM1 <- model1$UVPM
  l <- data$sn / 10
  UVC1 <- model1$UVC#[4500:7500, 90:100]
  }else{
    UVC_name <- paste0("Average_",group,"_UVC.rds")
    UVPM_name <- paste0("Average_",group,"_UVPM.rds")
    UV <- readRDS(UVC_name)
    if(group == 'fmci') UVC1 <- UV#[4500:5500, 90:100]
    else 
    UVC1 <- UV[, 90:100]
    UVPM1 <- readRDS(UVPM_name)
    l <- 200000/10
  }
  # model1$TAC
  
  mc1 <- mcmc(data = UVC1, start = 1, end = nrow(UVC1), thin = 1)
  setwd(outputdir)
  pdf(file = paste0(group, "_Traceplot.pdf"))
  traceplot(mc1, ylab = "Covariance Estimate")
  dev.off()
}

Heatmap <- function(g1, g2 = NA, metric = "OD", av = TRUE, order = TRUE, range = NA,
                    atlasloc = '~/Documents/Work/ModelFiles/finalAtlas.rds',
                    modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                    outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps_Final/"
                    ) {
  #' @param g1: The first group to be used in the heatmap
  #' @param g2: The second group to be used in the heatmap, if NA, it will only use g1
  #' @param av: If true, it will use the average data, if false, it will use the raw data from a model output. Use Average for most final plots 
  #' @param order: If true, it will order the heatmap by the regions of interest. If false, it will not order the heatmap. Be sure to have an atlas file
  #' @param range: For the subject heatmap, it allows you to give a legend range, if NA it will create one automatically
  #' @param modeldir: The directory where the model outputs are stored. If using av, point it to where the average data is stored
  #' @param outputdir: The directory where the heatmap will be saved\
  
  #' @description A universal function for making heatmaps. Given a group, and if to order, it will generate a heatmap for the model given. If a second group is given, it will give the difference between the two groups using computed CI's
  library(corrplot)
  Flag <- FALSE # Flag to determine if we are doing a comparison map or not
  
  center_matrix <- function(mat) {
    # Calculate the overall mean of all elements
    mat_mean <- mean(mat)
    
    # If the mean is already zero, return the matrix unchanged
    if(abs(mat_mean) <= 0.001)  {
      print("Mean is very close to 0")
      return(mat)
    }
    
    # Subtract the mean from all elements to center the matrix
    centered_mat <- mat - mat_mean
    print("Mean is not equal to 0, scaling matrix")
    print(paste("Old mean was ", mat_mean))
    return(centered_mat)
  }
  
  
  
  reorder <- function(matrix_to_reorder, atlasloc = '~/Documents/Work/ModelFiles/finalAtlas.rds' ) { #function to reorder the matrix based on a given atlas
    # Load the required data
    load(atlasloc)


    # Sort the data frame by the 'LOBE' column
    sorted_df <- final_df[order(final_df$LOBE), ]

    # Reorder the rows based on the sorted indices
    sorted_indices <- sorted_df$row_number

    # Reorder both rows and columns of the matrix
    reordered_mat <- matrix_to_reorder[sorted_indices, sorted_indices]

    return(reordered_mat)
  }

  add_RSN_borders <- function(atlasloc = '~/Documents/Work/ModelFiles/finalAtlas.rds') { #function to add the borders to the heatmap if we are ordering
    load(atlasloc)
    unique_groups <- unique(final_df[order(final_df$LOBE), ])
    group <- unique(unique_groups$LOBE)
    # Manually figure out the rectangles
    group_ranges <- list(
      list(name = "Cingulate", range = c(1, 8)),
      list(name = "Frontal", range = c(9, 30)),
      list(name = "Insular", range = c(31, 32)),
      list(name = "Occipital", range = c(33, 40)),
      list(name = "Parietal", range = c(41, 50)),
      list(name = "Subcortical", range = c(51, 64)),
      list(name = "Temporal", range = c(65, 82))
    )
    ymax <- 84
    ymin <- 1

    for (group in group_ranges) {
      range <- group$range
      flipped_ybottom <- ymax - (range[2] + 0.5 - ymin)
      flipped_ytop <- ymax - (range[1] - 0.5 - ymin)

      rect(
        xleft = range[1] - 0.5,
        ybottom = flipped_ybottom,
        xright = range[2] + 0.5,
        ytop = flipped_ytop,
        border = "black",
        lwd = 2
      )

      text_width <- strwidth(group$name, cex = 1.5) * 1.1 # Slightly larger for padding
      text_height <- strheight(group$name, cex = 1.5) * 1.1

      # Text position
      text_x <- range[2] + 10
      text_y <- .5 + (flipped_ybottom + flipped_ytop) / 2

      # Draw yellow rectangle as a highlight behind text
      rect(
        xleft = text_x - text_width / 2,
        ybottom = text_y - text_height / 2,
        xright = text_x + text_width / 2,
        ytop = text_y + text_height / 2,
        col = "yellow", # Highlight color
        border = NA # No border to keep it clean
      )

      # Draw text on top of the rectangle
      text(
        x = text_x,
        y = text_y,
        labels = group$name,
        cex = 1.5,
        col = "black",
        adj = c(0.5, 0.5) # Centered alignment
      )
    }
  }

  setwd(modeldir)
  o <- outputdir 
  if (is.na(g2)){
  if (av) { 
    # load the data for g1
    g1name <- paste0("Average_", g1, "_UVPM.rds")
    g1data <- readRDS(g1name)
  } else {
    # load the data for g1
    g1name <- paste0("ADNI_",metric,"_", g1, "_mean_5e+05_10000.rdata")
    data <- readRDS(g1name)
    model <- data$model
    g1data <- model$UVPM
    
    }
  }
  fname <- paste0(g1, "_heatmap2.pdf")
  if (!is.na(g2)) {
    Flag <- TRUE
    fname <- paste0(g1, "_vs_", g2, "_heatmap.pdf")
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
  

  if (order) m <- reorder(matrix_data, atlasloc = atlasloc)
  else m <- matrix_data
  
   

  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
  setwd(o)
   #12 if we are doing full size, 8.5 if side by side 
  if(Flag){
    library(corrplot)
    pdf(file = fname, width = 12)
    m[m > 1] <- 1
    m[m < -1] <- -1
    # Define the color palette: blue for -1, white for 0, orange for 1
    color_palette <- colorRampPalette(c("#1F77B4", "white", "#FF7F0E"))(200)
    # Dr. Wang gave the idea to match the color scheme depending on male and female, so I'll test that 
    if(g1 == 'mcn')
    {
      # Change the color palette to the male set that will be used in the overall Heatmap later on (MeasurePerformance.R)
      #08306B #Decrease in males
      #CB181D #Increase in males
      color_palette <- colorRampPalette(c("#08306B", "white", "#CB181D"))(200)
    }
    # Plot the matrix
    suppressWarnings({
    corrplot(m,
             method = "color",
             col = color_palette,
             tl.pos = "n",   # Remove axis labels
             cl.pos = "n",   # Remove color legend
             is.corr = FALSE,
             col.lim = c(-1, 1)   # This is important! Map -1 to +1
    )
    })
    
    # Add a custom legend
    if(g1 =='fcn'){
    legend("right", legend = c("Smaller", "Near no difference", "Larger"), 
           fill = c("#1F77B4", "#FFFFFF", "#FF7F0E"), 
           bty = "n", cex = 1.2)
    }else{
      legend("right", legend = c("Smaller", "Near no difference", "Larger"), 
             fill = c("#08306B", "#FFFFFF", "#CB181D"), 
             bty = "n", cex = 1.2)
    }
    
  }else {
    pdf(file = fname, width =8.5)
    
    if(is.na(range)) {
      range <- c(min(m, na.rm = TRUE), max(m, na.rm = TRUE))
    }
    suppressWarnings({ 
    corrplot(m,
             method = "color",
             col = colorRampPalette(c("#1F77B4","white", "#FF7F0E"))(200),
             tl.pos = "n", # Remove axis labels (numbers)
             # addgrid.col = "white",
             cl.pos = "b",
             is.corr = FALSE,
             col.lim = range,
             cl.cex = 1.25
    )
  })
  }
  
  
  if (order) add_RSN_borders(atlasloc = atlasloc)
  dev.off()
}


AverageCorr <- function(location = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/DimTesting/", group = "fcn",metric = 'OD', quiet = TRUE) {
    #' @title AverageCorr
    #' @description This function calculates the average correlation for a given group and metric. Best for Dimensionality Testing
    #' @param location: The location of the files to be used
    #' @param group: The group to be used in the analysis
    #' @param metric: The metric to be used in the analysis
    #' 
# The file structure that these files need to be in is as follows:
#' Give it a list of folders that are named 1-N,
#'  and then inside each of those folders,
#'  there are folders for each replication ('Rep-1', 'Rep-2', for example)
#' Inside each of those folders, 
#'  there are files named ADNI_[metric]_[GROUP]_mean_50000_1000.rdata, etc.

  
  setwd(location)
  file_name <- paste0("ADNI_",metric,"_", group, "_mean_50000_1000.rdata")

  # Get a matrix to store the correlations, 1-8 for Dim, 1-10 for each rep
  corr_matrix <- matrix(0, nrow = 8, ncol = 10)


  # for each Dim, we need to find the correlation for each rep
  # Get a list of all the folders in the directory
  Dimfolders <- list.dirs(location, full.names = FALSE, recursive = FALSE)
  # Loop through each folder
  for (folder in Dimfolders){
    Repfolders <- list.dirs(paste0(location,'/', folder), full.names = FALSE, recursive = FALSE)
    temp <- 1 #Used to store the index for rep, as the rep folders are not necessarily purely numerical
    for (rep in Repfolders){
      # Get the file name
      filename <- paste0(folder,"/",rep, "/", file_name)
      # Check if the file exists in the folder
      if (file.exists(filename)) {
        if(!quiet) print(paste0(filename, " Exists, computing Corr"))
        # Load the file
        data <- readRDS(filename)
        # Extract the model
        model1 <- data$model
        x_test <- data$testX
        l <- length(x_test)
        est <- model1$EFlPM
        est_split <- est[1:l]
        vec1 <- unlist(x_test)
        vec2 <- unlist(est_split)
        value <- cor(vec1, vec2)
        # Store the value in the matrix
        corr_matrix[as.numeric(folder), temp] <- value
      }else warning(paste("File does not exist: ", filename))
      
        temp <- temp + 1
    }
        }
  # Average the matrix for each Dim
  avg_corr <- rowMeans(corr_matrix, na.rm = TRUE)

  # We also want the SD of the matrix
  sd_corr <- apply(corr_matrix, 1, sd, na.rm = TRUE)
  
  
  # Save both of these as a rdata file, to be read later
  saveRDS(sd_corr, paste0(location, "/SD_", group, "_",metric,"_Corr.rds"))
  saveRDS(avg_corr, paste0(location, "/Average_", group, "_",metric,"_Corr.rds"))
}

AverageAPM <- function(g='fcn',metric='OD',modelloc = "/N/slate/conlcorn/SexLinkedProject/FinalModels/OD")
{
  # for each folder, we need to see if the file exists, if it does, save the APM. Once we have them all
  # We can compute the average and print it 
  # We should also save both the OG data, and the average as a file 
  setwd(modelloc)
  Dimfolders <- list.dirs(modelloc, full.names = FALSE, recursive = FALSE)
  APM <- list()
  for (folder in Dimfolders){
    filename <- paste0(folder,"/ADNI_",metric,"_",g,"_mean_2e+05_1000.rdata")
    if (file.exists(filename)) {
      data <- readRDS(filename)
      model <- data$model
      apm <- model$APM
      APM[[folder]] <- apm 
    }
  }
  #Convert to a matrix
  matrix_data <- matrix(unlist(APM), nrow = length(APM), byrow = TRUE)
  RowAv <- colMeans(matrix_data)
  average <- mean(RowAv)
  print(paste("Average APM for", g))
  print(average)
  
  file <- list("APM" = APM, "Average" = average)
  
  filename <- paste0("Average_",g,"_APM.rds")
  saveRDS(file,file=filename)
}

CI_SpiderPlot <- function(g1 = "fcn", g2 = "fmci",metric = 'OD',
                        modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                        outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/SpiderPlots", 
                        av = TRUE) 
{
  #' @param outputdir defines the output directory for the plots
  #' @param modelloc defines the location of the models
  #' @param g1 defines the first group
  #' @param g2 defines the second group
  #' This measures if we see where g1 is significantly less than g2, more than, or both

  
  reorder <- function(matrix_to_reorder, atlasloc = '~/Documents/Work/FinalFiles/finalAtlas.rds' ) { #function to reorder the matrix based on a given atlas
    # Load the required data
    load(atlasloc)
    
    
    # Sort the data frame by the 'LOBE' column
    sorted_df <- final_df[order(final_df$LOBE), ]
    
    # Reorder the rows based on the sorted indices
    sorted_indices <- sorted_df$row_number
    
    # Reorder both rows and columns of the matrix
    reordered_mat <- matrix_to_reorder[sorted_indices, sorted_indices]
    
    return(reordered_mat)
  }
  
  
  o <- paste0(outputdir)
  CI <- function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }
  # Defines both the lobe names, and where their connections are
  lobe_info <- list(
    list(name = "Cingulate", range = c(1, 8)),
    list(name = "Frontal", range = c(9, 30)),
    list(name = "Insular", range = c(31, 32)),
    list(name = "Occipital", range = c(33, 40)),
    list(name = "Parietal", range = c(41, 50)),
    list(name = "Subcortical", range = c(51, 64)),
    list(name = "Temporal", range = c(65, 82))
  )
  if (av) { # if we are using average data
    setwd(modelloc)
    model1name <- paste0("Average_", g1, "_UVC.rds")
    model2name <- paste0("Average_", g2, "_UVC.rds")
    UVC1 <- readRDS(model1name)
    UVC2_1 <- readRDS(model2name)
    if(g2 == 'fmci') UVC2 <- UVC2_1
    else UVC2 <- UVC2_1
  } else {
    setwd(modelloc)
    model1name <- paste0("ADNI_",metric,"_", g1, "_mean_2e+05_1000.rdata")
    model1 <- readRDS(model1name)
    model1 <- model1$model
    model2name <- paste0("ADNI_",metric,"_", g2, "_mean_2e+05_1000.rdata")
    model2 <- readRDS(model2name)
    model2 <- model2$model
    UVC2 <- model2$UVC
    UVC1 <- model1$UVC
  }
  # Generate the CI's for g1
  hi.g1 <- t(apply(UVC1, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) # returns a 84x83/2 matrix, which is the CI's for each unique connection
  # How to find which one is the unique region????
  # do the same for g2
  
  # Generate the CI's for g2
  hi.g2 <- t(apply(UVC2, 2, function(x) quantile(x, probs = c(0.025, 0.975))))
  
  # Look for regions in each lobe region in where the CI's do not overlap
  #   If they do not, add them to a list of regions that are different
  #   If they do, add them to a list of regions that are the same
  
  # It wouldn't be wise to reconstruct the matrix, it'll be easier to parse manually
  # We know the order of the 84x83/2 goes through each row of the original matrix, until it can't anymore, then it goes to the next column
  # This matrix stores where we are in the original matrix, and which index in the list responds to that location
  
  # Now, we can use this matrix to parse our data
  Greater_Important <- list(
    list(name = "Cingulate", value = NA),
    list(name = "Frontal", value = NA),
    list(name = "Insular", value = NA),
    list(name = "Occipital", value = NA),
    list(name = "Parietal", value = NA),
    list(name = "Subcortical", value = NA),
    list(name = "Temporal", value = NA)
  )
  
  Lower_Important <- list(
    list(name = "Cingulate", value = NA),
    list(name = "Frontal", value = NA),
    list(name = "Insular", value = NA),
    list(name = "Occipital", value = NA),
    list(name = "Parietal", value = NA),
    list(name = "Subcortical", value = NA),
    list(name = "Temporal", value = NA)
  )
  
  
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

  
  # Reorder the matrix based on lobes
  ordered_matrix <- reorder(matrix_data)
  for (lobe in lobe_info) {
    range <- lobe$range
    name <- lobe$name
    sub_matrix <- ordered_matrix[range[1]:range[2], range[1]:range[2]]
    greater <- sum(sub_matrix == 1)
    lesser  <- sum(sub_matrix == -1)
    # Find the index for the matching lobe name
    lobe_index <- which(sapply(Greater_Important, function(x) x$name == name))
    
    # Append values (can be replaced with any structure you want)
    Greater_Important[[lobe_index]]$value <- greater
    Lower_Important[[lobe_index]]$value   <- lesser
  }
  Greater_DF <- data.frame(
    Name = sapply(Greater_Important, function(x) x$name),
    TotalConnections = sapply(Greater_Important, function(x) x$value)
  )
  
  Lower_DF <- data.frame(
    Name = sapply(Lower_Important, function(x) x$name),
    TotalConnections = sapply(Lower_Important, function(x) x$value))
  
  
  max_val <- max(Greater_DF$TotalConnections, Lower_DF$TotalConnections)
  min_val <- 0
  
  # make a radar chart of the data using the fmsb library
  library(fmsb)
  
  radar_data <- as.data.frame(rbind(
    rep(max_val, nrow(Greater_DF)), # Max values
    rep(min_val, nrow(Greater_DF)), 
    Lower_DF$TotalConnections, # Actual values# Min values
    Greater_DF$TotalConnections
   
  ))
  
  colnames(radar_data) <- Greater_DF$Name # Set column names
  #chart_title <- paste0("Signifigant regions for ", g1, " and ", g2)
  # Create Spider Plot
  fname <- paste0(g1, "_", g2, "_spiderplot.pdf")
  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
  setwd(o)
  pdf(file = fname, width = 9)
  step <- ceiling(max_val / 4) # Ensure step is an integer
  # Define colors for colorblind-friendly scheme
  areas <- c(
    rgb(0, 1, 1, 0.25), # Cyan fill with transparency for the larger group
    rgb(1, 0.5, 0, 0.25) # Orange fill with transparency for the smaller group
  )
  
  # Generate correct axis labels
  axis_labels <- seq(0, max_val, by = step) # Generate sequence
  if (tail(axis_labels, 1) != max_val) {
    axis_labels <- c(axis_labels, max_val) # Ensure max_val is explicitly included
  }
  
  # Convert labels to character to prevent automatic formatting issues
  axis_labels <- as.character(axis_labels)
  
  # Create radar chart
  radarchart(radar_data,
             axistype = 1,
             pcol = c("cyan", "orange"), # Outline colors
             pfcol = areas, # Fill the inner area with colors
             plwd = 3,
             #          title = chart_title,
             vlcex = 2,
             cglcol = "gray",
             cglty = 2,
             axislabcol = "black",
             caxislabels = axis_labels # Force correct label sequence
  )

  # Add a legend
  if (g1 == 'fmci' && g2 == 'mmci') leglab = c("FMCI > MMCI", "FMCI < MMCI")
  else if (g1 == 'mmci' && g2 == 'fmci') leglab = c("FMCI < MMCI", "FMCI > MMCI")
  else leglab = c("CU > MCI", "CU < MCI")
  legend(
    x ="bottomright",
    legend = leglab,
    bty = "n",
    pch = 20,
    col = c("cyan", "orange"),
    pt.cex = 3,
    cex = 1.5,
    text.col = "black",
    horiz = FALSE
  )
  dev.off()
}



CI <- function(group = "f", metric = "OD", modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/",atlasloc , av = TRUE, outputdir, lobe = NA) {
  #This modified one instead is used to look at specific lobes, if no lobe is supplied, it will look at all lobes
  library(dplyr)
  setwd(modeldir) # Set the Directory to where we store the models
  
  if (!av) {
    cn_name <- paste0("ADNI_", metric, "_", group, "cn_mean_1e+05_1000.rdata")
    mci_name <- paste0("ADNI_", metric, "_", group, "mci_mean_1e+05_1000.rdata")
  } else {
    cn_UVC_name <- paste0("Average_", group, "cn_UVC.rds")
    cn_UVPM_name <- paste0("Average_", group, "cn_UVPM.rds")
    mci_UVC_name <- paste0("Average_", group, "mci_UVC.rds")
    mci_UVPM_name <- paste0("Average_", group, "mci_UVPM.rds")
  }
  
  credible_interval <- function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }
  
  # New helper: Adjust point estimate if it's outside its CI
  adjust_point_estimates <- function(estimates, intervals) {
    adjusted <- mapply(function(est, ci) {
      lower <- ci[1]
      upper <- ci[2]
      if (est < lower || est > upper) {
        return(mean(ci))
      } else {
        return(est)
      }
    }, estimates, split(intervals, row(intervals)))
    return(adjusted)
  }
  
  vec1.tmp <- array(rep(0, 84 * 84), dim = c(84, 84))

  
#if lobe is NA, we can use the bellow indeces
if(is.na(lobe)) {
  vec1.tmp[c(20, 36, 39, 43, 48, 69, 74, 80), ] <- 100
}else{
  load(atlasloc)
  sorted_df <- final_df[order(final_df$LOBE), ]
  # since we pulled the final_df out, it's still in memory, so we can use it to filter the rows and columns of the matrix
  lobe_indices <- which(final_df$LOBE == lobe) # Get the indices of the lobe we want to analyze
  vec1.tmp[lobe_indices,lobe_indices ] <- 100 # Set the rows of the lobe to 100
}
  vv <- vec1.tmp[upper.tri(vec1.tmp, diag = FALSE)]
  index <- which(vv != 0)
  
  # we can use this index to only set the ones we want to 100
  
  # Load CN data
  if (!av) {
    df <- readRDS(cn_name)
    model1 <- df$model
    UVC1 <- model1$UVC
    UVPM <- model1$UVPM
  } else {
    UVC1 <- readRDS(cn_UVC_name)
    UVPM <- readRDS(cn_UVPM_name)
  }
  hi.g1 <- t(apply(UVC1, 2, credible_interval))
  rhi.g1 <- hi.g1[index, ]
  pm <- UVPM[upper.tri(UVPM, diag = FALSE)]
  pm1 <- pm[index]
  order_pm1 <- order(pm1)
  sorted_pm1 <- pm1[order_pm1]
  negpe.g1 <- sorted_pm1[sorted_pm1 < 0]
  pospe.g1 <- sorted_pm1[sorted_pm1 > 0]
  sorted_rhi.g1 <- rhi.g1[order_pm1, ]
  negatives.g1 <- sorted_rhi.g1[sorted_pm1 < 0, ]
  positives.g1 <- sorted_rhi.g1[sorted_pm1 > 0, ]
  order_pos <- order_pm1[sorted_pm1 > 0]
  order_neg <- order_pm1[sorted_pm1 < 0]
  
  # Load MCI data
  if (!av) {
    df <- readRDS(mci_name)
    model1 <- df$model
    if(group == 'f') UVC1 <- model1$UVC
    else UVC1 <- model1$UVC
    UVPM <- model1$UVPM
  } else {
    UV <- readRDS(mci_UVC_name)
    if (group == 'f') UVC1 <- UV
    else UVC1 <- UV
    UVPM <- readRDS(mci_UVPM_name)
  }
  
  hi.g5 <- t(apply(UVC1, 2, credible_interval))
  rhi.g5 <- hi.g5[index, ]
  pm <- UVPM[upper.tri(UVPM, diag = FALSE)]
  pm5 <- pm[index]
  sorted_pm5 <- pm5[order_pm1]
  negpe.g5 <- sorted_pm5[sorted_pm1 < 0]
  pospe.g5 <- sorted_pm5[sorted_pm1 > 0]
  sorted_rhi.g5 <- rhi.g5[order_pm1, ]
  negatives.g5 <- sorted_rhi.g5[sorted_pm1 < 0, ]
  positives.g5 <- sorted_rhi.g5[sorted_pm1 > 0, ]
  
  # ✅ Adjust point estimates to fall inside CI if needed
  pospe.g1 <- adjust_point_estimates(pospe.g1, positives.g1)
  pospe.g5 <- adjust_point_estimates(pospe.g5, positives.g5)
  negpe.g1 <- adjust_point_estimates(negpe.g1, negatives.g1)
  negpe.g5 <- adjust_point_estimates(negpe.g5, negatives.g5)
  
  # Safeguard
  ensure_length_match <- function(x, y, target_length) {
    x <- head(x, target_length)
    y <- head(y, target_length)
    list(x = x, y = y)
  }
  
  plot_CI_95_dem <- function(size, pospe.g1, positives.g1, pospe.g5, positives.g5, order_pos, is_positive = TRUE) {
    offset <- 0.4
    colsuse <- c("steelblue", "orange2")
    B <- size
    
    # Combine all CI bounds and point estimates for full-range x-axis
    
    all_vals <- c(
      positives.g1[, 1], positives.g1[, 2], pospe.g1,
      positives.g5[, 1], positives.g5[, 2], pospe.g5
    )
    max_abs <- max(abs(all_vals), na.rm = TRUE)
    min_abs <- min(abs(all_vals), na.rm = TRUE)
    x_range <- range(all_vals, na.rm = TRUE)
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
      
      lengths <- ensure_length_match(x_vals, 1:B + offset * ifelse(j == 1, 1, -1), B)
      x_vals <- lengths$x
      y_vals <- lengths$y
      
      if (j == 1) {
        plot(x_vals, y_vals,
             pch = 20,
             xlim = x_range,                      # ✅ Expanded x-axis
             xlab = " ",
             ylab = "Connectivity Edge Index",    # ✅ Updated y-label
             cex = 1,
             col = colsuse[j],
             yaxt = "n"
        )
      } else {
        points(x_vals, y_vals,
               pch = 20,
               col = colsuse[j],
               yaxt = "n"
        )
      }
      
      for (i in 1:B) {
        segments(x_0[i], y_vals[i], x_1[i], y_vals[i], col = colsuse[j], lwd = 1.5)
      }
    }
    
    axis(side = 2, at = 1:B, labels = order_pos, cex.axis = 0.5, padj = 0.6, las = 2)
    axis(side = 1)
    legend("bottomright", legend = c("Healthy", "MCI"), lty = 1, col = colsuse, cex = 1)
  }
  o <- paste(outputdir,"/")
  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
  # Plot Positives
  pdf(paste0(o, group, metric,lobe, "CIPos.pdf"), height = 11, width = 7)
  message("Generating positive connectivity plot...")
  plot_CI_95_dem(min(length(pospe.g1), length(pospe.g5)), pospe.g1, positives.g1, pospe.g5, positives.g5, order_pos)
  dev.off()
  message("Positive connectivity plot saved.")
  
  # Plot Negatives
  pdf(paste0(o, group, metric,lobe, "CINeg.pdf"), height = 11, width = 7)
  message("Generating negative connectivity plot...")
  plot_CI_95_dem(min(length(negpe.g1), length(negpe.g5)), negpe.g1, negatives.g1, negpe.g5, negatives.g5, order_neg, is_positive = FALSE)
  dev.off()
  message("Negative connectivity plot saved.")
}

plotremake <- function(models,output, metric = 'OD', atlasloc)
# An overall function that regenerates all plots for the project 
# This is useful if we rerun the models. It makes Traceplots, Heatmaps, and spiderplots.
{

group <- c('fcn','fmci','mcn','mmci')

#Call the traceplot function for each group

#in the outputdir, we want to save it to the plots folder
outputdir <- paste0(output,"/Traceplots")
#for (g in group) Traceplot(modelname =  paste0("ADNI_OD_", g, "_mean_5e+05_10000.rdata"),group = g, filedir = models, outputdir = outputdir , av = FALSE)
#message("Traceplots Done, Making Heatmaps")

#Heatmaps
outputdir <- paste0(output,"/Heatmaps")
for (g in group)
{
  if(g == 'fcn' || g =='fmci' ) r <- c(-.17,.27)
  else r <- c(-.2,.32)
  Heatmap(g1 = g, metric = 'OD', modeldir = models, outputdir = outputdir, atlasloc = atlasloc, av = FALSE, order = TRUE, range = r)
  
}
#Also want the dif heatmaps for the two genders
Heatmap(g1 = 'fcn', g2 = 'fmci', metric = 'OD', modeldir = models, outputdir = outputdir, atlasloc = atlasloc, av = FALSE, order = TRUE)
Heatmap(g1 = 'mcn', g2 = 'mmci', metric = 'OD', modeldir = models, outputdir = outputdir, atlasloc = atlasloc, av = FALSE, order = TRUE)
message("Heatmaps Done, moving to spider plots")
#Spiderplots
outputdir <- paste0(output,"/Spiderplots")
CI_SpiderPlot(g1 = 'fcn', g2 = 'fmci', modelloc = models, outputdir = outputdir, av = FALSE)
CI_SpiderPlot(g1 = 'mcn', g2 = 'mmci', modelloc = models, outputdir = outputdir, av = FALSE)
#message("All plots generated and saved.")
}
plotremake(models = '~/Documents/Work/FinalFiles/500k',output ='~/Documents/Work/plotting/500k' ,atlasloc = '~/Documents/Work/FinalFiles/finalAtlas.rds')


# Lobe specific Heatmap function 
LobeHeatmap <- function(g1,g2, lobe = NA, modelloc, atlasloc, outputloc, av = TRUE) {
  # takes in two groups, and makes the four heatmaps we are planning to use for our overall lobe specific analysis 
  #' @param g1 defines the first group
  #' @param g2 defines the second group
  #' @param lobe defines the lobe to analyze, either 'Cingulate', 'Frontal', 'Insular', 'Occipital', 'Parietal', 'Subcortical', or 'Temporal'
  #' @param modelloc defines the location of the models
  #' @param atlasloc defines the location of the atlas file
  #' @param outputloc defines the location to save the results
  #' @param av defines if we are using the average data or not
  
  #We plan for a Raw UVPM differences and CI differences 
  library(corrplot)
  


  if(av){
  setwd(modelloc)
  g1name <- paste0('Average_',g1,'_UVPM.rds')
  g2name <- paste0('Average_',g2,'_UVPM.rds')
  
  g1UVPM <- readRDS(g1name)
  g2UVPM <- readRDS(g2name)
  
  
  g1UVC <- readRDS(paste0('Average_',g1,'_UVC.rds'))
  g2UVC <- readRDS(paste0('Average_',g2,'_UVC.rds'))
  
  if(g2 == 'fmci') { # We need to limit the MCMC to the stabalized regions(4500-5500)
    temp <- g2UVC
    g2UVC <- temp # Set the g2UVC to the temp value
  }else if(g1 == 'fmci') {
    temp <- g1UVC
    g1UVC <- temp # Set the g1UVC to the temp value
  }
  }

  else {
    setwd(modelloc)
    g1name <- paste0('ADNI_OD_',g1,'_mean_5e+05_10000.rdata')
    g2name <- paste0('ADNI_OD_',g2,'_mean_5e+05_10000.rdata')
    
    g1pull <- readRDS(g1name)
    g2pull <- readRDS(g2name)
    
    g1model <- g1pull$model
    g2model <- g2pull$model
    g1UVC <- g1model$UVC
    g2UVC <- g2model$UVC

    g1UVPM <- g1model$UVPM
    g2UVPM <- g2model$UVPM
    
  }
  # best way to do this is to do the operations to all of the data, then just subset the data to the lobe we want to analyze, since most of the code we have already written works off the whole set 
  
  # Generate the absolute differences
  abs_diff_UVPM <- g1UVPM - g2UVPM
  
  
  # Load the required data
  load(atlasloc)
  
  
  # Sort the data frame by the 'LOBE' column
  sorted_df <- final_df[order(final_df$LOBE), ]
  
  # Reorder the rows based on the sorted indices
  sorted_indices <- sorted_df$row_number
  
  # Reorder both rows and columns of the matri
  # We now need to filter based on our lobe of interest
  # since we pulled the final_df out, it's still in memory, so we can use it to filter the rows and columns of the matrix
  #If lobe is NA, we can just pull all of them
  if(is.na(lobe)) {
    lobe_indices <- 1:nrow(final_df) # Get all indices
  } else {
    lobe_indices <- which(final_df$LOBE == lobe) # Get the indices of the lobe we want to analyze
  }
  
  lobe_matrix <- abs_diff_UVPM[lobe_indices, lobe_indices] # Subset the matrix to only the lobe of interest
  
  # In theory, we now have a matrix that holds the absolute differences in the UVPM for the lobe of interest
  
  # We also need to get the name of each region in the lobe, so we can label the heatmap
  lobe_region_names <- final_df$ROI[lobe_indices] # Get the region names for the lobe of interest
  # in lobe region names, remove everything before the first underscore
  lobe_region_names <- sub("^[^_]+-", "", lobe_region_names) # Remove everything before the first underscore
  rownames(lobe_matrix) <- lobe_region_names
  colnames(lobe_matrix) <- lobe_region_names
  # Create a heatmap of the absolute differences in UVPM for the lobe of interest, using the lobe_region_names as the for each row/column
  # Define a custom color palette: cyan for negative, white for zero, orange for positive
  custom_palette <- colorRampPalette(c("#FF7F0E", "#FFFFFF", "#1F77B4"))(100)
  # Plot 1: Heatmap of differences
  #---------------------------------------------------------------------------------------------- 
  # Plot the heatmap with no numbers, axis titles only,
  # and vertical legend on the right
  
  # Open PDF device with larger margins to fit axis labels
  if (!dir.exists(outputloc)) dir.create(outputloc, recursive = TRUE)
  setwd(outputloc)
  fname <- paste0(g1, "_vs_", g2, "_", lobe, "_lobe_heatmap.pdf")
  pdf(file = fname, width = 10, height = 10) # Increase width/height for more space
  
  # Adjust plot margins: bottom, left, top, right
  
  
  # Flip the y-axis and corresponding labels, set blue for negative and orange for positive, show legend
  # Create a custom color palette: blue for negative, white for zero, orange for positive
  #custom_palette <- colorRampPalette(c("#1F77B4", "#FFFFFF", "#FF7F0E"))(100)
  
  #in the matrix, set the names of the rows and columns to the lobe region name
  corrplot(
    lobe_matrix,
    is.corr = FALSE,  # Set to FALSE since this is not a correlation matrix
    col = custom_palette,
    method = "color",  # Use color method for heatmap
    addgrid.col = "white",  # Add grid lines in White
    tl.col = "black",  # Text label color
    tl.srt = 45,  # Rotate text labels for better readability
    tl.cex = 1.3,  # Text label size
    cl.cex = 1.3,  # Color legend size
  )
  dev.off()
  # Part 2: Get the top five most significant connections within the lobe
  # -----------------------------------------------------------------------
  loc_matrix <- matrix(0, nrow = 84, ncol = 84)
  
  # Start filling the upper triangle downwards, column by column
  index <- 1
  for (j in 2:84) { # Start from the second column
    for (i in 1:(j - 1)) { # Only fill below the diagonal
      loc_matrix[i, j] <- index
      index <- index + 1
    }
  }
  
  
  
  lobe_region_names <- final_df$ROI[lobe_indices] # Get the region names for the lobe of interest
  # in lobe region names, remove everything before the first underscore
  # retrieve the standard deviation of the UVPM for each region in the lobe, using the UVC to calculate the standard deviation

  #UVC is a matrix of itt x total connections. We can parse it, reconstruct a 84x84 matrix to hold the SD, and then use that matrix to add to the differences 

  #calculate the standard deviation of the UVC
  lobe_UVC <- g1UVC[lobe_indices, lobe_indices] - g2UVC[lobe_indices, lobe_indices] # Get the UVC for the lobe of interest
 # combute the SD of each column
 # Calculate the standard deviation for each connection in the lobe
  # Should be a 1x3486 matrix, that we can then tie to a specific region 

  #put these into a 84x84 matrix 
  SDmatrix <- matrix(0, nrow = 84, ncol = 84)
  for (i in 1:84) {
    for (j in 1:84) {
      n <- loc_matrix[i, j]
      if (n == 0) next
      #place the standard deviation in the matrix
      #Calce the Variance of the UVC1
      var1 <- var(g1UVC[n, ])
      var2 <- var(g2UVC[n, ])

      #Calce the Covariance of the UVC1 and 2
      cov12 <- cov(g1UVC[n, ], g2UVC[n, ])

      # Calculate the standard deviation
      s <- sqrt(var1 + var2 - 2 * cov12) 
      SDmatrix[i, j] <- s
    }
  }


  #lobe_region_names <- sub("^[^_]+_", "", lobe_region_names) # Remove everything before the first underscore
  rownames(lobe_matrix) <- lobe_region_names
  colnames(lobe_matrix) <- lobe_region_names
  # Find the indices of the 5 smallest and 5 largest values (excluding diagonal)
  lobe_matrix_no_diag <- lobe_matrix
  diag(lobe_matrix_no_diag) <- NA
  # also set the Lower Tri to NA, so we only get the upper triangle
  lobe_matrix_no_diag[lower.tri(lobe_matrix_no_diag)] <- NA
  
  # Get indices for smallest and largest values
  smallest_idx <- arrayInd(order(lobe_matrix_no_diag, na.last = NA)[1:10], dim(lobe_matrix_no_diag))
  largest_idx <- arrayInd(order(lobe_matrix_no_diag, decreasing = TRUE, na.last = NA)[1:10], dim(lobe_matrix_no_diag))
  
  # Build data frames with region names and difference values
  five_smallest <- data.frame(
    Region1 = rownames(lobe_matrix)[smallest_idx[, 1]],
    Region2 = colnames(lobe_matrix)[smallest_idx[, 2]],
    Difference = lobe_matrix_no_diag[smallest_idx],
    SD = SDmatrix[smallest_idx]
  )
  
  five_largest <- data.frame(
    Region1 = rownames(lobe_matrix)[largest_idx[, 1]],
    Region2 = colnames(lobe_matrix)[largest_idx[, 2]],
    Difference = lobe_matrix_no_diag[largest_idx],
    SD = SDmatrix[largest_idx]
  )
  # Print the results
  print("Five smallest differences in UVPM within the lobe:")
  print(five_smallest)
  print("Five largest differences in UVPM within the lobe:")
  print(five_largest)
  # sort the lobe_matrix by the absolute difference, but keep track if it's positive or negative, we then want to print the top five most significant connections
  # We will use the absolute value of the differences to sort, but we will keep track of the sign
  
  
  
  
  
  
  # Part 3: Make the Signifigance Heatmap
  # -----------------------------------------------------------------------
  # generate the 95% CI for the UVC DATA
  hi.g1 <- t(apply(g1UVC, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) # 95% CI
  hi.g2 <- t(apply(g2UVC, 2, function(x) quantile(x, probs = c(0.025, 0.975)))) # 95% CI
  
  
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
  
  
  # we then need to reorder using the same method as before, so that the regions are in the correct order
  lobe_significant_matrix <- matrix_data[lobe_indices, lobe_indices] # Subset the matrix to only the lobe of interest
  # add the names
  lobe_region_names <- final_df$ROI[lobe_indices] # Get the region names for the lobe of interest
  # in lobe region names, remove everything before the first underscore
  lobe_region_names <- sub("^[^_]+-", "", lobe_region_names) # Remove everything before the first underscore
  rownames(lobe_significant_matrix) <- lobe_region_names
  colnames(lobe_significant_matrix) <- lobe_region_names
  m <- lobe_significant_matrix
  # Set the Y axis labels to " ", so we only show the X-axis labels
  
  
  #in the matrix, set the names of the rows and columns to the lobe region name
  m[m > 1] <- 1
  m[m < -1] <- -1
  # Define the color palette: blue for -1, white for 0, orange for 1
  color_palette <- colorRampPalette(c("#1F77B4", "white", "#FF7F0E"))(200)
  # Open PDF device with larger margins to fit axis labels
  if (!dir.exists(outputloc)) dir.create(outputloc, recursive = TRUE)
  setwd(outputloc)
  fname <- paste0(g1, "_vs_", g2, "_", lobe, "_lobe_significance_heatmap.pdf")
  pdf(file = fname, width = 10,height = 10) # Increase width/height for more space
  # Plot the matrix
  # remove the y axis labels, and only show the x-axis labels
  
  
  suppressWarnings({
    corrplot(
      m,
      is.corr = FALSE,  # Set to FALSE since this is not a correlation matrix
      col = color_palette,
      method = "color",  # Use color method for heatmap
      addgrid.col = "white",  # Add grid lines in White
      tl.col = "black",  # Text label color
      tl.srt = 45,  # Rotate text labels for better readability
      tl.cex = 1.3,  # Text label size
      cl.cex = 1.3,  # Color legend size
    )
  })
  dev.off() # Close the PDF device to save the plot
}
CI_SpiderPlot(g1 = 'fcn',g2 = 'fmci',metric = 'FA',modelloc = '~/Documents/Work/FinalFiles/Other',outputdir = '~/Documents/Work/FinalFiles/Other',av = FALSE)



















































