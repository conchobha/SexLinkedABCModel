# This file is planned to be a collection of final versions of functions that are used to gather analysis
# Working versions of functions will be in the MeasurePerformance.R file
# This file should be sourced and used by other analysis scripts


#' @title Measure Performance
#' @description This function measures the performance of a model using various metrics.



Heatmap <- function(g1, g2 = NA, av = TRUE,
                    modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                    outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps_Final/",
                    order = TRUE) {
  #' @param g1: The first group to be used in the heatmap
  #' @param g2: The second group to be used in the heatmap, if NA, it will only use g1
  #' @param av: If true, it will use the average data, if false, it will use the raw data from a model output
  #' @param modeldir: The directory where the model outputs are stored. If using av, point it to where the average data is stored
  #' @param outputdir: The directory where the heatmap will be saved\
  #' @param order: If true, it will order the heatmap by the regions of interest. If false, it will not order the heatmap
  #' @description A universal function for making heatmaps. Given a group, and if to order, it will generate a heatmap for the model given. If a second group is given, it will give the difference between the two groups
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
    print(paste("Old mean is ", mat_mean))
    return(centered_mat)
  }
  
  
  
  reorder <- function(matrix_to_reorder, wd = '~/Documents/Work/ModelFiles' ) { #function to reorder the matrix based on a given atlas
    # Load the required data
    setwd(wd)
    load("~/Documents/Work/ModelFiles/finalAtlas.rds")


    # Sort the data frame by the 'LOBE' column
    sorted_df <- final_df[order(final_df$LOBE), ]

    # Reorder the rows based on the sorted indices
    sorted_indices <- sorted_df$row_number

    # Reorder both rows and columns of the matrix
    reordered_mat <- matrix_to_reorder[sorted_indices, sorted_indices]

    return(reordered_mat)
  }

  add_RSN_borders <- function() { #function to add the borders to the heatmap if we are ordering
    load("~/Documents/Work/ModelFiles/finalAtlas.rds")
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
  if (av) { 
    # load the data for g1
    g1name <- paste0("Average_", g1, "_UVPM.rds")
    g1data <- readRDS(g1name)
  } else {
    # load the data for g1
    g1name <- paste0("ADNI_OD_", g1, "_mean_1e+05_1000.rdata")
    data <- readRDS(g1name)
    model <- data$model
    g1data <- model$UVPM
    
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
    model1name <- paste0("Average_", g1, "_UVC.rds")
    model2name <- paste0("Average_", g2, "_UVC.rds")
    UVC1 <- readRDS(model1name)
    UVC2_1 <- readRDS(model2name)
    if(g2 == 'fmci') UVC2 <- UVC2_1[4500:5500,]
    else UVC2 <- UVC2_1

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
        if (!max(tmp1[1], tmp2[1]) <= min(tmp1[2], tmp2[2])) {
          if (tmp1[2] <= tmp2[1]) {
            matrix_data[i, j] <- -1
          } else if (tmp1[1] >= tmp2[2]) {
            matrix_data[i, j] <- 1
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
  

  if (order) m <- reorder(matrix_data)
  else m <- matrix_data
  
   

  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
  setwd(o)
   #12 if we are doing full size, 8.5 if side by side 
  if(Flag){
    library(corrplot)
    pdf(file = fname, width = 12)
    
    # Define the color palette: blue for -1, white for 0, orange for 1
    color_palette <- colorRampPalette(c("#1F77B4", "white", "#FF7F0E"))(200)
    
    # Plot the matrix
    corrplot(m,
             method = "color",
             col = color_palette,
             tl.pos = "n",   # Remove axis labels
             cl.pos = "n",   # Remove color legend
             is.corr = FALSE,
             col.lim = c(-1, 1)   # This is important! Map -1 to +1
    )
    
    # Add a custom legend
    legend("right", legend = c("Smaller", "Near no difference", "Larger"), 
           fill = c("#1F77B4", "#FFFFFF", "#FF7F0E"), 
           bty = "n", cex = 1.2)
    
    
  }else {
    pdf(file = fname, width =8.5)
    corrplot(m,
             method = "color",
             col = colorRampPalette(c("#1F77B4","white", "#FF7F0E"))(200),
             tl.pos = "n", # Remove axis labels (numbers)
             # addgrid.col = "white",
             cl.pos = "b",
             is.corr = FALSE,
             col.lim = c(-0.2,.27),
             cl.cex = 1.25
    )
    }
  
  
  if (order) add_RSN_borders()
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

CI_SpiderPlot <- function(g1 = "fcn", g2 = "fmci",
                        modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                        outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/SpiderPlots", 
                        av = TRUE) 
{
  #' @param outputdir defines the output directory for the plots
  #' @param modelloc defines the location of the models
  #' @param g1 defines the first group
  #' @param g2 defines the second group
  #' This measures if we see where g1 is significantly less than g2, more than, or both

  
  reorder <- function(matrix_to_reorder) {
    # Load the required data
    load("~/Documents/Work/ModelFiles/finalAtlas.rds")
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
    if(g2 == 'fmci') UVC2 <- UVC2_1[4500:5500,]
    else UVC2 <- UVC2_1
  } else {
    setwd(modelloc)
    model1name <- paste0("ADNI_OD_", g1, "_mean_2e+05_1000.rdata")
    model1 <- readRDS(model1name)
    model1 <- model1$model
    model2name <- paste0("ADNI_OD_", g2, "_mean_2e+05_1000.rdata")
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
      if (!max(tmp1[1], tmp2[1]) <= min(tmp1[2], tmp2[2])) {
        if (tmp1[2] <= tmp2[1]) {
          matrix_data[i, j] <- -1
        } else if (tmp1[1] >= tmp2[2]) {
          matrix_data[i, j] <- 1
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

disgraph <- function(loc = '~/Downloads/APM') { # violin plot code 
  setwd(loc)
  grouplist <- c("fcn", "fscd", "fmci", "mcn", "mscd", "mmci")  # Change this to a vector
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
    fscd = "Female SCD",
    fmci = "Female MCI",
    mcn = "Healthy Male",
    mscd = "Male SCD",
    mmci = "Male MCI"
  )
  
  # Count number of datapoints per group
  counts <- datalist_df %>%
    group_by(group) %>%
    summarise(n = n(), .groups = 'drop')
  
  # Create group labels with counts
  counts$label <- paste0(group_labels[counts$group], " (n = ", counts$n, ")")
  
  # Join updated labels into datalist_df
  datalist_df <- datalist_df %>%
    left_join(counts, by = "group")
  
  # Set the order of the group factor to preserve the grouplist order
  datalist_df$label <- factor(datalist_df$label, levels = paste0(group_labels[grouplist], " (n = ", counts$n, ")"))
  
  # Violin plot with group-specific n values in x-axis labels
  ggplot(datalist_df, aes(x = label, y = value)) +
    geom_violin(fill = "gray85", color = "black", trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
    theme_minimal() +
    labs(y = "Average APM", x = "") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 18, angle = 45, hjust = 1),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold")
    )
}

CI_SpiderPlot(modelloc ='~/Documents/Work/ModelFiles' , outputdir = '~/Downloads/WorkPlots')
CI(group = 'f', modeldir = '~/Documents/Work/ModelFiles', outputdir = '~/Downloads/WorkPlots')


