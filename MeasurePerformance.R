library(coda)
library(MASS)
library(boot)
library(latentnet)
library(lvm4net)
library(amen)
library(pROC)

Heatmap <- function(g1, g2 = NA, av = TRUE,
                    modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                    outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps_Final/",
                    order = TRUE) {
  #@param g1: The first group to be used in the heatmap
  #@param g2: The second group to be used in the heatmap, if NA, it will only use g1
  #@param av: If true, it will use the average data, if false, it will use the raw data from a model output
  #@param modeldir: The directory where the model outputs are stored. If using av, point it to where the average data is stored
  #@param outputdir: The directory where the heatmap will be saved\
  #@param order: If true, it will order the heatmap by the regions of interest. If false, it will not order the heatmap


  # A universal function for making heatmaps. Given a group, and if to order, it will generate a heatmap for the model given
  # If a second group is given, it will give the difference between the two groups
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

subjectHeatmap <- function(group) {
  library(corrplot)
  # Goes to each participant, makes a heatmap of the data, and saves it in a folder
  o <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/ParticMapsor/", group)
  if (!dir.exists(o)) dir.create(o, recursive = TRUE)

  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features")
  features <- paste0(group, "_OD_mean.rds")
  labeldata <- readRDS(features)[[1]]
  setwd(o)
  for (i in 1:length(labeldata)) {
    outname <- paste0(i, "_Heatmap.pdf")
    pdf(file = outname)
    data <- labeldata[[i]]
    # Clip the data to be within the range 0 to 0.6
    data_clipped <- pmin(pmax(data, 0), 0.6)
    # Now, create the corrplot with manually set breaks and color limits
    corrplot(reorder_mat(data_clipped),
      method = "color", # Use color method for heatmap
      col = colorRampPalette(c("blue", "white", "red"))(200),
      col.lim = c(0, .6),
      tl.pos = "lt", # Place axis labels on the left and top
      tl.col = "black", # Set axis numbers color to black
      tl.cex = 0.3, # Adjust the font size of the axis numbers
      addgrid.col = "white", # Set grid lines to white
      number.cex = 0.1,
      cl.pos = "b", # Place color legend at the bottom
      is.corr = FALSE, # Force the color scale limits to always be between 0 and 0.6
    )
    dev.off()
  }
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
  UVC1 <- model1$UVC[4500:7500, 90:100]
  }else{
    UVC_name <- paste0("Average_",group,"_UVC.rds")
    UVPM_name <- paste0("Average_",group,"_UVPM.rds")
    UV <- readRDS(UVC_name)
    if(group == 'fmci') UVC1 <- UV[4500:5500, 90:100]
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
Correlation <- function(modelname = "ADNI_Dr_combined_mean_10000_1000.rdata",
                        modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output",
                        rep = 1) {
  setwd(modelloc)
  data <- readRDS(modelname)
  model1 <- data$model
  x_test <- data$testX
  l <- length(x_test)
  est <- model1$EFlPM
  est_split <- est[1:l]
  vec1 <- unlist(x_test)
  vec2 <- unlist(est_split)
  correlations <- cor(vec1, vec2)
  # View the correlation
  print(paste0("Correlation for Rep ", rep, ": ", correlations))
  return(correlations)
}
histo <- function(modelname, metric,
                  modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output",
                  outputloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots") { # makes a histogram for the specified model
  setwd(modelloc)
  data <- readRDS(modelname)
  model1 <- data$model
  dt <- unlist(model1$EFlPM)
  setwd(outputloc)
  pdf(file = paste0(metric, "_histogram.pdf"))
  hist(dt,
    main = paste0("Frequency for values for the EFLPM: ", metric),
    col = "darkmagenta"
  )

  dev.off()
}
dimtest <- function(dm) {
  start <- Sys.time()
  # Saves a PDF file, that can display the Correclation for each group, accross different
  # Dimensionality of a model.
  # Can Average accross multiple replications as well
  output <- matrix(data = NA, nrow = 8, ncol = 5)
  # Defining our 3d Output Matrix
  # These will define our DM, Row, and Column
  toutput <- array(dim = c(8, 5, 10)) # First index is Dim, Next is group, last is rep

  # We need to load each model, get the correlation, and the store it in the right spot
  dim_list <- list(1, 2, 3, 4, 5, 6, 7, 8)
  group_list <- list("fcn", "fmci", "mcn", "mmci", "combined")
  for (i in 1:10) { # loops over each Replication
    for (d in dim_list) {
      location <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/model_outputs/Dimensonality_testing/DM-", d, "/Rep_", i)
      setwd(location)
      count <- 1
      for (g in group_list) {
        # load the file
        filename <- paste0("ADNI_", dm, "_", g, "_mean_50000_1000.rdata")

        data <- tryCatch(
          {
            readRDS(filename)
          },
          error = function(e) {
            print(paste0("Error reading file: ", filename, " for DIM", d, " - Setting to NULL"))
            return(NULL)
          }
        )
        if (!is.null(data)) {
          model1 <- data$model
          x_test <- data$testX
          l <- length(x_test)
          est <- model1$EFlPM
          est_split <- est[1:l]
          vec1 <- unlist(x_test)
          vec2 <- unlist(est_split)
          correlations <- cor(vec1, vec2)
          # get the correlation

          # store it in the matrix
          toutput[d, count, i] <- correlations
        } else {
          output[d, count] <- NA
        }
        count <- count + 1
      }
    }
  }
  # Average together all of the last dim together into the output matrix
  # Average together all of the last dim together into the output matrix
  for (g in 1:5) { # For each Group set
    for (d in 1:8) { # For each Dim Set
      temp <- numeric(10) # Initialize a numeric vector with length 10
      for (i in 1:10) {
        temp[i] <- toutput[d, g, i] # Store values directly as numeric
      }
      # Get the Average of these, excluding NA values
      tav <- mean(temp, na.rm = TRUE)
      output[d, g] <- tav
    }
  }


  # Call the function to build the table with this data
  library(grid)
  library(gridExtra)
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/Dim")
  rownames(output) <- c("1", "2", "3", "4", "5", "6", "7", "8") # Dimensions
  colnames(output) <- c(
    "Female Healthy", "Female MCI", "Combined",
    "Male Healthy", "Male MCI"
  )
  name <- paste0("Dim_Testing_", dm, "_rep10.pdf")
  pdf(name, width = 12, height = 6)
  grid.newpage()
  grid.text(paste0("Dim Testing for ", dm, ". Averaged Across 10 Replications"), y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))

  # Convert the matrix to a table and display it
  grid.table(output)

  # Close the PDF device
  dev.off()
  end <- Sys.time()
  s <- as.numeric(difftime(end, start, units = "secs"))
  m <- as.numeric(difftime(end, start, units = "mins"))
  message(paste0("Execution Completed in ", s, " Seconds/ ", m, " Minutes"))
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
      "For Itterations of ", sn, ", It Took ",
      round(Comptime, 0), " Seconds/ ",
      round(Comptime / 60, 2), " Minutes/ ",
      round(Comptime / 3600, 2), " Hours/ ",
      round(Comptime / 3600 / 24, 2), " Days"
    ))
  }
}
# Confidence intervals
# Prints a confidence interval for the given group, comparing the UVPM between the healthy and unhealthy participants in those groups
CI <- function(group = "f", metric = "OD", modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/", av = TRUE, outputdir) {
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
  vec1.tmp[c(20, 36, 39, 43, 48, 69, 74, 80), ] <- 100
  vv <- vec1.tmp[upper.tri(vec1.tmp, diag = FALSE)]
  index <- which(vv != 0)
  
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
    if(group == 'f') UVC1 <- model1$UVC[4500:7500, ]
    else UVC1 <- model1$UVC
    UVPM <- model1$UVPM
  } else {
    UV <- readRDS(mci_UVC_name)
    if (group == 'f') UVC1 <- UV[4500:5500, ]
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
  pdf(paste0(o, group, metric, "CIPos.pdf"), height = 11, width = 7)
  message("Generating positive connectivity plot...")
  plot_CI_95_dem(min(length(pospe.g1), length(pospe.g5)), pospe.g1, positives.g1, pospe.g5, positives.g5, order_pos)
  dev.off()
  message("Positive connectivity plot saved.")
  
  # Plot Negatives
  pdf(paste0(o, group, metric, "CINeg.pdf"), height = 11, width = 7)
  message("Generating negative connectivity plot...")
  plot_CI_95_dem(min(length(negpe.g1), length(negpe.g5)), negpe.g1, negatives.g1, negpe.g5, negatives.g5, order_neg, is_positive = FALSE)
  dev.off()
  message("Negative connectivity plot saved.")
}



# Finds and Saves the ten regions in the Estimated Connectivity in which has the highest difference between two given models
find_difference <- function(res = 0, model1 = "ADNI_OD_fmci_mean_2e+05_1000.rdata", model2 = "ADNI_OD_mmci_mean_2e+05_1000.rdata", modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/Rep-1") {
  library(stringr)
  # Determine the Two groups we are comparing and Metric
  matches <- str_match(model1, "^ADNI_([^_]+)_([^_]+)_")[, 2:3]
  metric <- matches[1]
  group1 <- matches[2]
  matches <- str_match(model2, "^ADNI_([^_]+)_([^_]+)_")[, 3]
  group2 <- matches
  # Load Model1 and Store the UVPM data
  setwd(modeldir)

  data <- readRDS(model1)
  model <- data$model
  matrix_data1 <- model$UVPM

  data <- readRDS(model2)
  model <- data$model
  matrix_data2 <- model$UVPM


  # We need to find both the Positve And Negative Values
  pos_result <- matrix_data1 - matrix_data2
  neg_result <- matrix_data2 - matrix_data1
  abs_result <- abs(matrix_data1 - matrix_data2)
  # find the top 10 combinations
  if (res == 1) {
    result_matrix <- pos_result
  } else if (res == -1) {
    result_matrix <- neg_result
  } else {
    result_matrix <- abs_result
  }
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
    Value = round(upper_values[sorted_indices], 4),
    Row = row_col_indices[, 1],
    Column = row_col_indices[, 2]
  )
  print(top_10_with_positions)
  # Print the result
  setwd(paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/HeatMaps/", metric))
  pdf(file = paste0(group1, group2, "_heatmap.pdf"))
  corrplot(result_matrix,
    method = "color", # Use color method for heatmap
    col = colorRampPalette(c("white", "orange", "red", "black"))(200), # Color gradient
    tl.pos = "lt", # Place axis labels on the left and top
    tl.col = "black", # Set axis numbers color to black
    tl.cex = 0.3, # Adjust the font size of the axis numbers
    addgrid.col = "white", # Set grid lines to white
    number.cex = 0.1,
    cl.pos = "b", # Place color legend at the bottom
    is.corr = FALSE
  ) # Indicate this is not a correlation matrix
  dev.off()
}

pull_connect <- function(group1, group2, edge1, edge2) {
  # Pulls the connectivity data for a given edge, from each group, and is returned as a dataframe

  # Make the filenames for the two groups
  group1name <- paste0(group1, "_OD_mean.rds")
  group2name <- paste0(group2, "_OD_mean.rds")
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/data/processedData/Features/")

  # pull the data from both files

  data1 <- readRDS(group1name)
  data2 <- readRDS(group2name)
  list1 <- list()
  list2 <- list()
  # Add the data to the list
  for (person in data1[[1]]) list1 <- append(list1, person[edge1, edge2])
  for (person in data2[[1]]) list2 <- append(list2, person[edge1, edge2])
  listcomb <- list()
  listcomb[[1]] <- list1
  listcomb[[2]] <- list2
  # Create the dataframe

  return(listcomb)
}

find_smallest <- function(g1 = "fcn", g2 = "fmci", test = "t") {
  smallest_p_values <- data.frame(p_value = numeric(0), i = integer(0), j = integer(0))
  for (i in 1:84) {
    for (j in 1:84) {
      d <- pull_connect(g1, g2, i, j)
      sublist1 <- d[[1]] # List of length 13
      sublist2 <- d[[2]] # List of length 32

      # Access individual vectors within the sublists
      sublist1 <- unlist(sublist1)
      sublist2 <- unlist(sublist2)

      # Perform a t-test
      if (test == "t") {
        t_test_result <- t.test(sublist1, sublist2)
      } else {
        t_test_result <- var.test(sublist1, sublist2)
      }
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
  print(paste0("For ", g1, " and ", g2, "\n"))
  print(smallest_p_values)
}


findEdge <- function(number, filename = "iADRC_Struture_Diffusion_Tau_Abeta_84ROIcombo.xlsx", dataloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/", order = FALSE) {
  # Function to search for a given edge name, given its number
  library(readxl)
  if (!order) {
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
    # make a 2d matrix to hold the results

    # Go to the Index
    result <- df$ICV[df$`...5` == number]

    # Return the String
    rm(df)
    return(result[9]) # Files have a structure R does not like, this works fine however
  } else {



  }
}

# Function to sort through a given group, and find the average standard deviation for that group
# in their raw data.
findavSD <- function(group) {
  # create the 84x84 matrix
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
        s <- sd(list1, na.rm = TRUE) # Use `na.rm = TRUE` to ignore NAs
        SDs[i, j] <- s
      } else {
        SDs[i, j] <- NA # Assign NA if list1 is empty (no valid data)
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
  if (length(flat_SDs) > 0) {
    avg_SD <- mean(flat_SDs)
    print(paste0("group: ", group))
    # Find and print the highest and lowest SD
    max_SD <- max(flat_SDs)
    min_SD <- min(flat_SDs)
    print(paste0("Average SD: ", avg_SD))
    print(paste0("Highest SD: ", max_SD))
    print(paste0("Lowest SD: ", min_SD))
  } else {
    print("No valid data to calculate the average SD.")
  }
}








CI_analysis <- function(g1 = "fcn", g2 = "fmci",
                        modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD",
                        type = "c",
                        outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/SpiderPlots", av = TRUE) 
  {
  #' @param outputdir defines the output directory for the plots
  #' @param modelloc defines the location of the models
  #' @param g1 defines the first group
  #' @param g2 defines the second group
  #' @param type defines the type of analysis. Can be l,m,a,c
  #' This measures if we see where g1 is significantly less than g2, more than, or both
  #' c denotes a combined plot

  reorder <- function(matrix_to_reorder) {
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


  o <- paste0(outputdir, "/", type)
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
    UVCg <- readRDS(model2name)
    UVC2 <- UVCg[4500:5500,]
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

  r <- reorder(loc_matrix)
  loc_matrix <- r

  # Now, we can use this matrix to parse our data
  Greater_Important <- list(
    list(name = "Cingulate", regions = list()),
    list(name = "Frontal", regions = list()),
    list(name = "Insular", regions = list()),
    list(name = "Occipital", regions = list()),
    list(name = "Parietal", regions = list()),
    list(name = "Subcortical", regions = list()),
    list(name = "Temporal", regions = list())
  )

  Lower_Important <- list(
    list(name = "Cingulate", regions = list()),
    list(name = "Frontal", regions = list()),
    list(name = "Insular", regions = list()),
    list(name = "Occipital", regions = list()),
    list(name = "Parietal", regions = list()),
    list(name = "Subcortical", regions = list()),
    list(name = "Temporal", regions = list())
  )


  # if we are making a combined plot
  for (lobe in lobe_info) {
    # Get the range of values
    range <- lobe$range
    # Get the name of the lobe
    name <- lobe$name

    # create empty list
    Lower_regions <- list()
    Greater_regions <- list()
    # We need to parse the reordered loc matrix, and use the correct index to check the CI's
    for (i in range[1]:range[2]) { # Parse for each column
      for (j in range[1]:range[2]) { # Parse for each Row
        # if it is zero, we can pass this one
        n <- loc_matrix[i, j]
        if (n == 0) next
        # if it isnt, we need to grab the CI stored at that point
        tmp1 <- hi.g1[n, ]
        tmp2 <- hi.g2[n, ]
        # check if they overlap
        if (!max(tmp1[1], tmp2[1]) <= min(tmp1[2], tmp2[2])) {
          if (tmp1[2] <= tmp2[1]) {
            Lower_regions <- append(Lower_regions, list(c(i, j)))
          } else if (tmp1[1] >= tmp2[2]) {
            Greater_regions <- append(Greater_regions, list(c(i, j)))
          }
        }
      }
    }
    # add the lower and greater to the total list
    Greater_Important[[which(sapply(Greater_Important, function(x) x$name) == name)]]$regions <- Greater_regions
    Lower_Important[[which(sapply(Lower_Important, function(x) x$name) == name)]]$regions <- Lower_regions
  }

  Greater_DF <- data.frame(
    Name = sapply(Greater_Important, function(x) x$name),
    TotalConnections = sapply(Greater_Important, function(x) length(x$regions))
  )

  Lower_DF <- data.frame(
    Name = sapply(Lower_Important, function(x) x$name),
    TotalConnections = sapply(Lower_Important, function(x) length(x$regions))
  )

  max_val <- max(Greater_DF$TotalConnections, Lower_DF$TotalConnections)
  min_val <- 0

  # make a radar chart of the data using the fmsb library
  library(fmsb)

  radar_data <- as.data.frame(rbind(
    rep(max_val, nrow(Greater_DF)), # Max values
    rep(min_val, nrow(Greater_DF)), # Min values
    Greater_DF$TotalConnections,
    Lower_DF$TotalConnections # Actual values
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
  
  # Update legend with new colors
  #legend("topright",
   #      legend = c("Larger", "Smaller"),
    #     bty = "n", pch = 20, col = c("cyan", "orange"),
     #    text.col = "grey25", pt.cex = 2
  #)
  


  dev.off()
}








# Function to search through every folder in a DIR, find a given file, then return the average UVC or UVPM for that group

MakeAverage <- function(location = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/", est = "UVC", group = "fcn", metric = 'OD') {
  #' @param est defines the type of estimate to use, either UVC,UVPM, or TAC
  # Load the files
  setwd(location)
  # Get a list of all the folders in the directory
  folders <- list.dirs(location, full.names = FALSE, recursive = FALSE)
  # Initialize a list to store the data
  data_list <- list()
  # Loop through each folder
 dim <- 2
  for (folder in folders)
  {
    # Get the file name
    filename <- paste0("ADNI_",metric,"_", group, "_mean_2e+05_1000.rdata")
    # Check if the file exists in the folder
    if (file.exists(paste0(folder, "/", filename))) {
      # Load the file
      data <- readRDS(paste0(folder, "/", filename))
      # Extract the model
      model <- data$model
      input <- model$input
      if(input$K != 2) stop("You have the wrong dim in the file structure")
      # Extract the UVC or UVPM
      if (est == "UVC") {
        data_list[[folder]] <- model$UVC
      } else if (est == "UVPM") {
        data_list[[folder]] <- model$UVPM
      } else if (est == "THETAPM") {
        data_list[[folder]] <- model$THETAPM
      } else if(est == "TAC") {
        data_list[[folder]] <- model$TAC
      } else {
        stop("Invalid estimate type. Use 'UVC', 'UVPM', or 'TAC'.")
      }
      
    }
  }

  # make the average matrix, and save it as a file in location
  if (length(data_list) > 0) {
    # Get the dimensions of the matrix
    dim1 <- dim(data_list[[1]])
    # Initialize a matrix to store the average
    average_matrix <- matrix(0, nrow = dim1[1], ncol = dim1[2])
    # Loop through each matrix and add it to the average
    for (matrix in data_list)
    {
      average_matrix <- average_matrix + matrix
    }
    # Divide by the number of matrices to get the average
    average_matrix <- average_matrix / length(data_list)
    # Save the average matrix as an RDS file
    print("Saving file")
    saveRDS(average_matrix, paste0(location, "/Average_", group, "_", est,".rds"))
  }
}

MakeAverage(location = '/N/slate/conlcorn/SexLinkedProject/FinalModels/OD_DM1',est = 'UVPM',group = 'fcn')
MakeAverage(location = '/N/slate/conlcorn/SexLinkedProject/FinalModels/OD_DM1',est = 'UVC',group = 'fcn')
MakeAverage(location = '/N/slate/conlcorn/SexLinkedProject/FinalModels/OD_DM1',est = 'UVPM',group = 'fscd')

grouplist <- list("fcn", "fmci", "mcn", "mmci")

for(g in grouplist){
  MakeAverage(location = '/N/slate/conlcorn/SexLinkedProject/OldModels/Dim2Models',est = 'UVC',group = g)
  MakeAverage(location = '/N/slate/conlcorn/SexLinkedProject/OldModels/Dim2Models',est = 'UVPM',group = g)
}

 AverageCorr <- function(location = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/DimTesting/", group = "fcn",metric = 'OD') {
  # Plan:
  # Build a matrix of the Correlation in each model
  # Find the average of each model
  # Save the average as a file
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
    # for each rep in that folder
    for (rep in 1:15)
    {
      # Get the file name
      filename <- paste0(folder,"/Rep-",rep, "/", file_name)
      # Check if the file exists in the folder
      if (file.exists(filename)) {
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
        corr_matrix[as.numeric(folder), rep] <- value
      }else warning(paste("File does not exist: ", filename))
      }
    # average all results, and save as a file
  }
  # Average the matrix for each Dim
  avg_corr <- rowMeans(corr_matrix, na.rm = TRUE)

  # We also want the SD of the matrix
  sd_corr <- apply(corr_matrix, 1, sd, na.rm = TRUE)
  
  
  # Save both of these as a rdata file, to be read later
  saveRDS(sd_corr, paste0(location, "/SD_", group, "_",metric,"_Corr.rds"))
  saveRDS(avg_corr, paste0(location, "/Average_", group, "_",metric,"_Corr.rds"))
}



grouplist <- list("fcn", "fmci", 
                  "fscd", "mcn", "mmci", "mscd"
                  )

for(g in grouplist){
  setwd("/N/slate/conlcorn/SexLinkedProject/DimTesting/OD")
  #AverageCorr(location ="/N/slate/conlcorn/SexLinkedProject/DimTesting/OD" , group = g,metric = 'OD')
  cordata <- readRDS(paste0("Average_", g, "_OD_Corr.rds"))
  sddata <- readRDS(paste0("SD_", g, "_OD_Corr.rds"))
  print(g)
  cat(paste(round(cordata, 4), " \\pm ", round(sddata, 3)))
  
  print("\n")
}




replist <- c(1:10)
group <- "fcn"

grouplist <- list("fcn", "fmci", 
                  "fscd", "mcn", "mmci", "mscd")

for (g in grouplist) {
  print("---------------------")
  print(g)
  print("---------------------")
  for (r in replist) {
    # set wd to the head folder
    setwd("/N/slate/conlcorn/SexLinkedProject/FinalModels/OD")
    # build the folder and file name
    modelloc <- paste0("/N/slate/conlcorn/SexLinkedProject/FinalModels/OD/Rep-", r)
    filename <- paste0("ADNI_OD_", g, "_mean_2e+05_1000.rdata")
    # call the correlation function
    Correlation(modelname = filename, modelloc = modelloc, rep = r)
  }
}



# Function to Print the top five regions in terms of differences between groups
difference <- function(g1, g2, modelloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/", av = TRUE,atlasloc) {
  #We are going to find the difference between two groups, and print the top five regions in terms of difference
  #' @param g1 defines the first group
  #' @param g2 defines the second group
  #' @param modelloc defines the location of the models
  #' @param av defines if we are using the average data or not
# Load the averaged Data

#Determine the Abs Difference


#Add the difference to a dataframe that also holds the Region name

# Print the top 10 regions in terms of difference
GetRegion <- function(indices, atlasloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/finalAtlas.rds") {
  # This function will take a vector of indices and return a 2-column matrix: index and region name
  load(atlasloc)
  # Sort final df by the lobe
  region_names <-final_df$ROI[indices]
  result <- cbind(Index = indices, Region = region_names)
  return(result)
}
  
  if (av) {
    # Load the average data
    setwd(modelloc)
    model1name <- paste0("Average_", g1, "_UVPM.rds")
    model2name <- paste0("Average_", g2, "_UVPM.rds")
    UVC1 <- readRDS(model1name)
    UVC2 <- readRDS(model2name)
  } else {
    # Load the individual data
    setwd(modelloc)
    model1name <- paste0("ADNI_OD_", g1, "_mean_2e+05_1000.rdata")
    model2name <- paste0("ADNI_OD_", g2, "_mean_2e+05_1000.rdata")
    model1 <- readRDS(model1name)
    model2 <- readRDS(model2name)
    UVC1 <- model1$model$UVPM
    UVC2 <- model2$model$UVPM
  }
  
  # Calculate the absolute difference
  abs_diff <- abs(UVC1 - UVC2)

  # set the lower tri to -1
  abs_diff[lower.tri(abs_diff)] <- -1 # Set the lower triangle to -1 to ignore it in the results
  # load the atlas data




  if(av){
    model1name <- paste0("Average_", g1, "_UVC.rds")
    model2name <- paste0("Average_", g2, "_UVC.rds")
    UVC1_1 <- readRDS(model1name)
    UVC2_1 <- readRDS(model2name)
  }else{
    g1name <- paste0("ADNI_",metric,"_", g1, "_mean_2e+05_1000.rdata")
    data <- readRDS(g1name)
    model <- data$model
    UVC1_1 <- model$UVC
    
    g2name <- paste0("ADNI_",metric,"_", g2, "_mean_2e+05_1000.rdata")
    data <- readRDS(g2name)
    model <- data$model
    UVC2_1<- model$UVC
  }

    UVC  <- UVC1_1 - UVC2_1 # This is the difference between the two groups, we will use this to find the regions that are significantly different

    CI <- t(apply(UVC, 2, function(x) quantile(x, probs = c(0.025, 0.975))))

    loc_matrix <- matrix(0, nrow = 84, ncol = 84)
  
  # Start filling the upper triangle downwards, column by column
    index <- 1
    for (j in 2:84) { # Start from the second column
      for (i in 1:(j - 1)) { # Only fill below the diagonal
        loc_matrix[i, j] <- index
        index <- index + 1
     }
    }

    # loc_matrix is now a matrix that holds the index that can be used to match the CI's to the regions
    CIMatrix <- matrix(vector("list", 84 * 84), nrow = 84, ncol = 84)
    
    for (i in 1:84) {
      for (j in 1:84) {
        n <- loc_matrix[i, j]
        if(n == 0) next
        CIMatrix[[i, j]] <- unname(c(CI[n, 1], CI[n, 2]))  # Store vector in a cell
      }
    }
    
  # Create a data frame to hold the differences and region name 1 and 2
  # turn the abs_diff into a into a list that holds 3 values, the index(x and y) in which it is stored, and the difference value 
  # Create a data frame to hold Region 1, Region 2, and the difference value
    diff_df <- data.frame()
    
    for (i in 1:nrow(abs_diff)) {
      for (j in 1:ncol(abs_diff)) {
        if (abs_diff[i, j] == -1) next
        
        ci_vals <- CIMatrix[[i, j]]
        
        # Skip if no CI was stored
        if (is.null(ci_vals)) next
        
        diff_df <- rbind(diff_df, data.frame(
          Region1 = i,
          Region2 = j,
          Difference = abs_diff[i, j],
          CI_low = ci_vals[1],
          CI_high = ci_vals[2]
        ))
      }
    }

# We also want to compute the CI's for these regions, so we can see if they are significantly different as well as numerically different
# load the UVC data for both groups

 

  regionlookup <- GetRegion(1:84, atlasloc = atlasloc) # we can refrence this value to get the region name
  # Convert the list to a data frame
  regionlookup <- as.data.frame(regionlookup)



  #take the number stored in the first two values, and replace them with the region name from the lookup table
  diff_df$Region1 <- regionlookup$Region[match(diff_df$Region1, regionlookup$Index)]
  diff_df$Region2 <- regionlookup$Region[match(diff_df$Region2, regionlookup$Index)]

  #print the top 5 regions in terms of difference
  # sort by the difference value
  diff_df <- diff_df[order(-diff_df$Difference), ]
  top_regions <- head(diff_df, 5)
  print(paste("Top 5 regions in terms of difference between", g1, "and", g2))
  print(top_regions)

}


avAPM <- function(g='fcn',metric='OD',modelloc = "/N/slate/conlcorn/SexLinkedProject/FinalModels/OD")
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

for(g in grouplist) avAPM(g = g,modelloc = '/N/slate/conlcorn/SexLinkedProject/OldModels/Dim2Models')


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


grouplist <- list("fcn", "fmci", 
                  "fscd", "mcn", "mmci", "mscd")

for (g1_idx in 1:(length(grouplist) - 1)) {
  for (g2_idx in (g1_idx + 1):length(grouplist)) {
    g1 <- grouplist[g1_idx]
    g2 <- grouplist[g2_idx]
    APMTesting(g1 = g1, g2 = g2)
  }
}


AttributeCICompare <- function(g1,g2,modelloc,av = TRUE) # CI's seem to be too wide to gleam any insights
{
  # Function that takes in two regions, and checks the THETAPM for both attributes, and sees which regions are display the largest difference between groups, 
  # It then returns both the number, the difference in CI, and the name of the region 
  
  if(!av) {
    warning("Not Supported Yet")
    return()
    }
  else{
    g1name <- paste0('Average_',g1,'_TAC.rds')
    g2name <- paste0('Average_',g2,'_TAC.rds')
    
    g1data <- readRDS(g1name)
    g2data <- readRDS(g2name)
  }
  
  # Generate the CI's 
  CI <- function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }
  
  hi.g1 <- t(apply(g1data, 2, function(x) quantile(x, probs = c(0.05, 0.95)))) # 90% CI
  hi.g2 <- t(apply(g2data, 2, function(x) quantile(x, probs = c(0.05, 0.95)))) # 90% CI
  
  # Determine if they overlap
  Tau <- matrix(0,nrow=84)
  Amy <- matrix(0,nrow=84)
  for(i in 1:84)#84 total regions
  {
    tauindex <- i*2-1
    AmyIndex <- i*2
    #Tau
    tmp1 <- hi.g1[tauindex, ]
    tmp2 <- hi.g2[tauindex, ]
    if (max(tmp1[1], tmp2[1]) > min(tmp1[2], tmp2[2])) { #checks if they don't overlap
      if (tmp1[2] <= tmp2[1]) { #checks if g2 is larger or smaller than g1
        Tau[i,] <- 1
      } else if (tmp1[1] >= tmp2[2]) {
        Tau[i,] <- -1
      } 
    }
    
    #Amy
    tmp1 <- hi.g1[AmyIndex, ]
    tmp2 <- hi.g2[AmyIndex, ]
    if (max(tmp1[1], tmp2[1]) > min(tmp1[2], tmp2[2])) { #checks if they don't overlap
      if (tmp1[2] <= tmp2[1]) { #checks if g2 is larger or smaller than g1
        Amy[i,] <- 1
      } else if (tmp1[1] >= tmp2[2]) {
        Amy[i,] <- -1
      } 
    }
    
  }
  
  #current matrices now have the values of 1, -1, or 0 for each 84 regions

  # get the count of -1's and 1's for each matrix
  TauCount <- sum(Tau == 1)
  AmyCount <- sum(Amy == 1)
  TauCountNeg <- sum(Tau == -1)
  AmyCountNeg <- sum(Amy == -1)
  # Print the results
  print(paste("For group ", g1, " vs ", g2))
  print(paste("Tau: ", TauCount, " regions are larger in ", g2, " and ", TauCountNeg, " regions are larger in ", g1))
  print(paste("Amy: ", AmyCount, " regions are larger in ", g2, " and ", AmyCountNeg, " regions are larger in ", g1))
  
}


AttributeAbsCompare <- function(g1,g2,modelloc, av = TRUE) # Measures the Abs differences of attributes, and returns the 5 most signifigant regions
{
  if(!av) {
    warning("Not Supported Yet")
    return()
    }
  else{
    g1name <- paste0('Average_',g1,'_THETAPM.rds')
    g2name <- paste0('Average_',g2,'_THETAPM.rds')
    
    g1data <- readRDS(g1name)
    g2data <- readRDS(g2name)
  }
  
  matrix <- abs(g1data - g2data) # Get the absolute difference between the two groups
  # Get the top 5 regions with the largest absolute difference in each column

  Tau <- matrix[,1] # Get the Tau columns
  Amy <- matrix[,2] # Get the Amy columns
  # Add a column for the region index number
  Tau <- cbind(Region = 1:84, Tau)
  Amy <- cbind(Region = 1:84, Amy)
  
  # sort the Tau and Amy matrices by the absolute difference in descending order
  Tau <- Tau[order(abs(Tau[,2]), decreasing = TRUE), ]
  Amy <- Amy[order(abs(Amy[,2]), decreasing = TRUE), ]
  # Get the top 5 regions for each column
  top_Tau <- Tau[1:5, ]
  top_Amy <- Amy[1:5, ]
  # replace the region index with the region name using the GetRegion function
  top_Tau <- GetRegion(top_Tau[,1])
  top_Amy <- GetRegion(top_Amy[,1])
  # Print the results
  print(paste("For group ", g1, " vs ", g2))
  print("Top 5 regions for Tau:")
  print(top_Tau)
  print("Top 5 regions for Amy:")
  print(top_Amy)
}

GetRegion <- function(indices, atlasloc = "/N/u/conlcorn/BigRed200/SexLinkedProject/data/finalAtlas.rds") {
  # This function will take a vector of indices and return a 2-column matrix: index and region name
  load(atlasloc)
  # Sort final df by the lobe
  region_names <-final_df$ROI[indices]
  result <- cbind(Index = indices, Region = region_names)
  return(result)
}


LobebarPlot <- function(g1 = 'fcn',g2 = 'fmci', modelloc, outputloc,atlas, av = TRUE){
# Finds the aboslute differences in the estimated connectivty between each group, broken
# up by lobe, and then plots the differences in a bar plot
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
setwd(modelloc)
  if(!av) {
    warning("Not Supported Yet")
    return()
  }
  else{
    g1name <- paste0('Average_',g1,'_UVPM.rds')
    g2name <- paste0('Average_',g2,'_UVPM.rds')
    
    g1data <- readRDS(g1name)
    g2data <- readRDS(g2name)
  }
  
  # Generate the absolute differences
  abs_diff <- abs(g1data - g2data)
  # Reorder the matrix to match the final atlas order
  matrixreordered <- reorder(abs_diff,atlas)
  abs_diff <- matrixreordered 
  # Define lobe ranges
  lobe_info <- list(
    list(name = "Cingulate", range = c(1, 8)),
    list(name = "Frontal", range = c(9, 30)),
    list(name = "Insular", range = c(31, 32)),
    list(name = "Occipital", range = c(33, 40)),
    list(name = "Parietal", range = c(41, 50)),
    list(name = "Subcortical", range = c(51, 64)),
    list(name = "Temporal", range = c(65, 82))
  )
  
  # Initialize a vector to store the differences for each lobe
  lobe_diffs <- numeric(length(lobe_info))
  
  # Calculate the absolute differences for each lobe
  for (i in seq_along(lobe_info)) {
    lobe <- lobe_info[[i]]
    range <- lobe$range
    lobe_diffs[i] <- sum(matrixreordered[range[1]:range[2], range[1]:range[2]])
  }
  
  # Create a bar plot of the differences
  barplot(lobe_diffs,
          names.arg = sapply(lobe_info, function(x) x$name),
          main = paste("Absolute Differences in Connectivity\nbetween", g1, "and", g2),
          ylab = "Absolute Difference",
          col = "lightblue",
          las = 2)
  
  # Save the plot
  if (!dir.exists(outputloc)) dir.create(outputloc, recursive = TRUE)
  setwd(outputloc)
  
  fname <- paste0(g1, "_vs_", g2, "_lobe_differences.pdf")
  
  pdf(file = fname, width = 10, height = 6)
  dev.off() # Close the PDF device to save the plot
}


brainplotregiondifferences <- function (g1 = 'fcn',g2 = 'fmci', modelloc, outputloc,atlas, av = TRUE)
{
  # take a sum of the absolute differences in the estimated connectivty between each group, broken
  # up by brain region, and then plots the differences in a brain plot using ggseg

  # trim the matrix to only the upper triangle
trim_upper_triangle <- function(matrix) {
  # Create a copy of the matrix
  trimmed_matrix <- matrix
  
  # Set the lower triangle to NA
  trimmed_matrix[lower.tri(trimmed_matrix)] <- NA
  
  return(trimmed_matrix)
}

setwd(modelloc)
  if(!av) {
    warning("Not Supported Yet")
    return()
  }
  else{
    g1name <- paste0('Average_',g1,'_UVPM.rds')
    g2name <- paste0('Average_',g2,'_UVPM.rds')
    
    g1data <- readRDS(g1name)
    g2data <- readRDS(g2name)
  }
  
  # Generate the absolute differences
  abs_diff <- abs(g1data - g2data)
  m <- trim_upper_triangle(abs_diff) # Trim the matrix to only the upper triangle
  # for each row, find the sum of the absolute differences, save them in a df that holds
    # the name of the region, and the sum of the absolute differences, and its index number
  region_sums <- data.frame(Index = integer(0), Region = character(0), Sum = numeric(0))
  # Load the atlas data
  #remember, region 36-49 are the subcortical regions, so we will ignore them when we plot
  load(atlas)
  # fill the index and region columns
  for (i in 1:84) {
    region_name <- final_df$ROI[i] # Get the region name
    region_sum <- sum(m[i, ], na.rm = TRUE) # Get the sum of the absolute differences for that region
    # Append to the data frame
    region_sums <- rbind(region_sums, data.frame(Index = i, Region = region_name, Sum = region_sum))
  }

  # Filter out the subcortical regions (36-49)
  region_sums <- region_sums[!region_sums$Index %in% 36:49, ]
  # also remove Index 35 and 84
  region_sums <- region_sums[!region_sums$Index %in% c(35, 84), ]
  # we need to take this df to then plot the value in ggseg3d for the dk atlas
  library(ggseg3d)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  # Create a ggseg3d plot
  # change all _ in regions to -
  region_sums$Region <- gsub("-", "_", region_sums$Region)
  # Merge region_sums with dk_3d region names to ensure correct mapping
  someData <- dk_3d %>%
    #filter(surf == "inflated") %>%
    unnest(ggseg_3d) %>%
    ungroup() %>%
    select(label) %>%
    na.omit() %>%
    left_join(region_sums, by = c("label" = "Region")) %>%
    mutate(Sum = ifelse(is.na(Sum), 0, Sum))
     # Set missing regions to 0
  # Now we can plot the data

  ggseg3d(.data = someData,
          atlas = dk_3d,
          colour = "Sum", text = "Sum") 
  
}
