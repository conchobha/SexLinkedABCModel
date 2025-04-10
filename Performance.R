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
  reorder <- function(matrix_to_reorder) { #function to reorder the matrix based on a given atlas
    # Load the required data
    load("/N/u/conlcorn/BigRed200/SexLinkedProject/data/finalAtlas.rds")


    # Sort the data frame by the 'LOBE' column
    sorted_df <- final_df[order(final_df$LOBE), ]

    # Reorder the rows based on the sorted indices
    sorted_indices <- sorted_df$row_number

    # Reorder both rows and columns of the matrix
    reordered_mat <- matrix_to_reorder[sorted_indices, sorted_indices]

    return(reordered_mat)
  }

  add_RSN_borders <- function() { #function to add the borders to the heatmap if we are ordering
    load("/N/u/conlcorn/BigRed200/SexLinkedProject/data/finalAtlas.rds")
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
    # load data for g2, if we are doing a comparison map
    if (av) {
      g2name <- paste0("Average_", g2, "_UVPM.rds")
      g2data <- readRDS(g2name)
    } else {
      g2name <- paste0("ADNI_OD_", g2, "_mean_1e+05_1000.rdata")
      data <- readRDS(g2name)
      model <- data$model
      g2data <- model$UVPM
    }
    Flag <- TRUE
    fname <- paste0(g1, "_vs_", g2, "_heatmap.pdf")
    matrix_data <- g1data - g2data
  } else matrix_data <- g1data # if we are only looking at the heatmap of one group
  

  if (order) m <- reorder(matrix_data)
  else m <- matrix_data
  

  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
  setwd(o)
   #12 if we are doing full size, 8.5 if side by side 
  if(Flag){
    library(corrplot)
    pdf(file = fname, width =12)
    
    # Define a matrix of colors based on sign, making 0 white
    col_matrix <- ifelse(abs(m) < 0.05, "white", 
                         ifelse(m < 0, "#1F77B4", "#FF7F0E"))
    
    # Plot the full correlation matrix (no triangle mask)
    corrplot(m,
             method = "color",
             col = col_matrix,
             tl.pos = "n",   # Remove axis labels (numbers)
             cl.pos = "n",   # Remove default color legend
             is.corr = FALSE
    )
    
    # Add a custom legend on the right
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
             col.lim = c(-0.2,.305),
             cl.cex = 1.25
    )
    
    }
  
  
  if (order) add_RSN_borders()
  dev.off()
}


AverageCorr <- function(location = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/DimTesting/", group = "fcn",metric = 'OD') {
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
  getCorr <- function(filename)
  {
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
        return(value)
  }
  
  setwd(location)
  file_name <- paste0("ADNI_",metric,"_", group, "_mean_50000_1000.rdata")

  # Get a matrix to store the correlations, 1-8 for Dim, 1-10 for each rep
  corr_matrix <- matrix(0, nrow = 8, ncol = 15)


  # for each Dim, we need to find the correlation for each rep
  # Get a list of all the folders in the directory
  Dimfolders <- list.dirs(location, full.names = FALSE, recursive = FALSE)
  # Loop through each folder
  for (folder in Dimfolders)
  {
    Repfolders <- list.dirs(paste0(location, folder), full.names = FALSE, recursive = FALSE)
    temp <- 1 #Used to store the index for rep, as the rep folders are not necessarily purely numerical
    for (rep in Repfolders)
    {
      # Get the file name
      filename <- paste0(folder,"/",rep, "/", file_name)
      # Check if the file exists in the folder
      if (file.exists(filename)) {
        # Load the file
        value <- getCorr(filename)
        # Store the value in the matrix
        corr_matrix[as.numeric(folder), temp] <- value
      }else warning(paste("File does not exist: ", filename))
      }
        temp <- temp + 1
  }
  # Average the matrix for each Dim
  avg_corr <- rowMeans(corr_matrix, na.rm = TRUE)

  # We also want the SD of the matrix
  sd_corr <- apply(corr_matrix, 1, sd, na.rm = TRUE)
  
  
  # Save both of these as a rdata file, to be read later
  saveRDS(sd_corr, paste0(location, "/SD_", group, "_",metric,"_Corr.rds"))
  saveRDS(avg_corr, paste0(location, "/Average_", group, "_",metric,"_Corr.rds"))
}

