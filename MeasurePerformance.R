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
  # A universal function for making heatmaps. Given a group, and if to order, it will generate a heatmap for the model given
  # If a second group is given, it will give the difference between the two groups
  library(corrplot)

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

  add_RSN_borders <- function() {
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

      text_width <- strwidth(group$name, cex = 1.5) * 1.1  # Slightly larger for padding
      text_height <- strheight(group$name, cex = 1.5) * 1.1 
      
      # Text position
      text_x <- range[2] + 10
      text_y <- .5+(flipped_ybottom + flipped_ytop) / 2
      
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
  if (av){
    # load the data for g1
    g1name <- paste0("Average_", g1, "_UVPM.rds")
    g1data <- readRDS(g1name)
  }else{
  # load the data for g1
  g1name <- paste0("ADNI_OD_", g1, "_mean_1e+05_1000.rdata")
  data <- readRDS(g1name)
  model <- data$model
  g1data <- model$UVPM
  }
  fname <- paste0(g1, "_heatmap.pdf")
  if (!is.na(g2)) {
    # load data for g2, if we are doing a comparison map
    if(av){
      g2name <- paste0("Average_", g2, "_UVPM.rds")
      g2data <- readRDS(g2name)
      matrix_data <- g1data - g2data
    }else{
    g2name <- paste0("ADNI_OD_", g2, "_mean_1e+05_1000.rdata")
    data <- readRDS(g2name)
    model <- data$model
    g2data <- model$UVPM
    matrix_data <- g1data - g2data
    fname <- paste0(g1, "_vs_", g2, "_heatmap.pdf")
    }
  }else {
    matrix_data <- g1data
  }

  if (order) {
    m <- reorder(matrix_data)
  } else {
    m <- matrix_data
  }

  if (!dir.exists(o)) dir.create(o, recursive = TRUE)
  setwd(o)
  pdf(file = fname, width = 9)

  # Needs to be update eventually to be tracked via the regions. I have the .py code for that, need to integrate.
  corrplot(m,
    method = "color",
    col = colorRampPalette(c("blue", "white", "red"))(200),
    tl.pos = "n", # Remove axis labels (numbers)
    #addgrid.col = "white",
    cl.pos = "b",
    is.corr = FALSE,
    cl.cex = 1.25
  )
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



Traceplot <- function(modelname = "ADNI_Da_combined_mean_10000_1000.rdata", group, filedir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/Rep-1/") {
  setwd(filedir)
  data <- readRDS(modelname)
  model1 <- data$model
  UVPM1 <- model1$UVPM
  UVC1 <- model1$UVC[4500:7500, 90:100]
  # model1$TAC
  l <- data$sn / 10
  mc1 <- mcmc(data = UVC1, start = 1, end = nrow(UVC1), thin = 1)
  setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/Traceplots")
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
CI <- function(group = "f", metric = "OD", modeldir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/Rep-1/") {
  library(dplyr)
  setwd(modeldir) # Set the Directory to where we store the models
  cn_name <- paste0("ADNI_", metric, "_", group, "cn_mean_1e+05_1000.rdata") # file name for the cn group
  mci_name <- paste0("ADNI_", metric, "_", group, "mci_mean_1e+05_1000.rdata") # file name for mci group

  credible_interval <- function(x) {
    quantile(x, probs = c(0.025, 0.975))
  }

  vec1.tmp <- array(rep(0, 84 * 84), dim = c(84, 84))
  vec1.tmp[c(20, 36, 39, 43, 48, 69, 74, 80), ] <- 100
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
  UVC1 <- model1$UVC[4500:7500, ]
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
    colsuse <- c("steelblue", "orange2")
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
        plot(x_vals, y_vals,
          pch = 20,
          # xlim = c(-0.4, 0.4),
          # ylim = c(0,B),
          xlab = "Connectivity",
          ylab = "Different positions",
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
        if (between(0, x_0[i], x_1[i])) {
          segments(x_0[i], y_vals[i], x_1[i], y_vals[i], col = colsuse[j], lwd = 1.5)
        } else {
          segments(x_0[i], y_vals[i], x_1[i], y_vals[i], col = colsuse[j], lwd = 1.5)
        }
      }
    }

    axis(side = 2, at = 1:B, labels = order_pos, cex.axis = 0.5, padj = 0.6, las = 2)
    axis(side = 1)
    legend("bottomright", legend = c("Healthy", "Dementia"), lty = 1, col = colsuse, cex = 1)
  }

  pdf(paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/CIs/", group, metric, "CIPos.pdf"), height = 11, width = 7)
  message("Generating positive connectivity plot...")
  plot_CI_95_dem(min(length(pospe.g1), length(pospe.g5)), pospe.g1, positives.g1, pospe.g5, positives.g5, order_pos)
  dev.off()
  message("Positive connectivity plot saved.")

  ### Negative group
  pdf(paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/CIs/", group, metric, "CINeg.pdf"), height = 11, width = 7)
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
                        outputdir = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/plots/SpiderPlots", av = TRUE) {
  #' @param outputdir defines the output directory for the plots
  #' @param modelloc defines the location of the models
  #' @param g1 defines the first group
  #' @param g2 defines the second group
  #' @param type defines the type of analysis. Can be l,m,a,c
  #' This measures if we see where g1 is significantly less than g2, more than, or both
  #' c denotes a combined plot
  #' @TODO: Add support for combined plots
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
  if(av){ # if we are using average data 
    setwd(modelloc)
    model1name <- paste0("Average_", g1, "_UVC.rds")
    model2name <- paste0("Average_", g2, "_UVC.rds")
    UVC1 <- readRDS(model1name)
    UVC2 <- readRDS(model2name)
  }else{
 
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

  # Initialize lists to store results
  important_regions <- list(
    list(name = "Cingulate", regions = list()),
    list(name = "Frontal", regions = list()),
    list(name = "Insular", regions = list()),
    list(name = "Occipital", regions = list()),
    list(name = "Parietal", regions = list()),
    list(name = "Subcortical", regions = list()),
    list(name = "Temporal", regions = list())
  )
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

  # Now, we can use this matrix to parse our data
  if (type != "c") {
    for (lobe in lobe_info) {
      # Get the range of values
      range <- lobe$range
      # Get the name of the lobe
      name <- lobe$name

      # create empty list
      regions <- list()
      # if we aren't combining the two plots

      for (i in range[1]:range[2]) { # Parse for each column
        for (j in range[1]:range[2]) { # Parse for each Row
          # check if we are in an out of range spot using the loc_matrix
          n <- loc_matrix[i, j]
          if (n == 0) next
          # Get the CI's for each group using n
          tmp1 <- hi.g1[n, ]
          tmp2 <- hi.g2[n, ]

          # check if they overlap
          if (!max(tmp1[1], tmp2[1]) <= min(tmp1[2], tmp2[2])) {
            if (type == "l" && tmp1[2] <= tmp2[1]) {
              regions <- append(regions, list(c(i, j)))
            } else if (type == "g" && tmp1[1] >= tmp2[2]) {
              regions <- append(regions, list(c(i, j)))
            } else if (type == "a") regions <- append(regions, list(c(i, j)))
          }
        }
      }
      # Add the regions to the important regions list
      important_regions[[which(sapply(important_regions, function(x) x$name) == name)]]$regions <- regions

      df_connections <- data.frame(
        Name = sapply(important_regions, function(x) x$name),
        TotalConnections = sapply(important_regions, function(x) length(x$regions))
      )

      max_val <- max(df_connections$TotalConnections)
      min_val <- 0

      # make a radar chart of the data using the fmsb library
      library(fmsb)

      radar_data <- as.data.frame(rbind(
        rep(max_val, nrow(df_connections)), # Max values
        rep(min_val, nrow(df_connections)), # Min values
        df_connections$TotalConnections # Actual values
      ))

      colnames(radar_data) <- df_connections$Name # Set column names
      if (type == "a") {
        chart_title <- paste0("Signifigant regions for ", g1, " and ", g2)
      } else if (type == "l") {
        chart_title <- paste0(g1, " is less than ", g2)
      } else if (type == "g") chart_title <- paste0(g1, " is greater than ", g2)
      # Create Spider Plot
      fname <- paste0(g1, "_", g2, "_", type, "_spiderplot.pdf")
      if (!dir.exists(o)) dir.create(o, recursive = TRUE)
      setwd(o)
      pdf(file = fname)
      step <- ceiling(max_val / 4)  # Ensure step is an integer
      
      areas <- c(rgb(1, 0, 0, 0.25),  # Red fill with transparency
                 rgb(0, 1, 0, 0.25))  # Green fill with transparency
      
      radarchart(radar_data,
                 axistype = 1,
                 pcol = 2:3,
                 pfcol = areas,  # Fill the inner area with colors
                 plwd = 3,
                 title = chart_title, 
                 vlcex = 2,
                 cglcol = "gray", 
                 cglty = 2, 
                 axislabcol = "black",
                 caxislabels = seq(0, max_val, by = step)  # Correctly spaced integer labels
      )
      

      dev.off()
    }
  } else {
    Greater_Important <- list(
      list(name = "Cingulate", regions = list()),
      list(name = "Frontal", regions = list()),
      list(name = "Insular", regions = list()),
      list(name = "Occipital", regions = list()),
      list(name = "Parietal", regions = list()),
      list(name = "Subcortical", regions = list()),
      list(name = "Temporal", regions = list()))
      
    Lower_Important <- list(
        list(name = "Cingulate", regions = list()),
        list(name = "Frontal", regions = list()),
        list(name = "Insular", regions = list()),
        list(name = "Occipital", regions = list()),
        list(name = "Parietal", regions = list()),
        list(name = "Subcortical", regions = list()),
        list(name = "Temporal", regions = list()))
    
    
    # if we are making a combined plot
        for (lobe in lobe_info) {
      # Get the range of values
      range <- lobe$range
      # Get the name of the lobe
      name <- lobe$name

      # create empty list
      Lower_regions <- list()
      Greater_regions <- list()
      # if we aren't combining the two plots

      for (i in range[1]:range[2]) { # Parse for each column
        for (j in range[1]:range[2]) { # Parse for each Row
          # check if we are in an out of range spot using the loc_matrix
          n <- loc_matrix[i, j]
          if (n == 0) next
          # Get the CI's for each group using n
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
      
      # Add the regions to the important regions list
      Greater_Important[[which(sapply(Greater_Important, function(x) x$name) == name)]]$regions <- Greater_regions
      Lower_Important[[which(sapply(Lower_Important, function(x) x$name) == name)]]$regions <- Lower_regions
      
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
      chart_title <- paste0("Signifigant regions for ", g1, " and ", g2)
      # Create Spider Plot
      fname <- paste0(g1, "_", g2,"_spiderplot.pdf")
      if (!dir.exists(o)) dir.create(o, recursive = TRUE)
      setwd(o)
      pdf(file = fname,width = 9)
      step <- ceiling(max_val / 4)  # Ensure step is an integer
      
      areas <- c(rgb(1, 0, 0, 0.25),  # Red fill with transparency
                 rgb(0, 1, 0, 0.25))  # Green fill with transparency
      
      # Generate axis labels as character strings to prevent formatting issues
      axis_labels <- as.character(seq(0, max_val, by = step))
      
      radarchart(radar_data,
                 axistype = 1,
                 pcol = 2:3,
                 pfcol = areas,  # Fill the inner area with colors
                 plwd = 3,
                 title = chart_title, 
                 vlcex = 2,
                 cglcol = "gray", 
                 cglty = 2, 
                 axislabcol = "black",
                 caxislabels = axis_labels  # Ensure max value appears correctly
      )
      
      
      
      
      legend("topright",
             legend = c("Larger","Smaller"),
             bty = "n", pch = 20, col = areas,
             text.col = "grey25", pt.cex = 2)
      
      
      dev.off()
    }
  
  
  
  }
}

CI_analysis(g1 = 'fscd',g2 = "mmci")
grouplist <- list("fcn", "fmci", "fscd", "mcn", "mmci", "mscd")

for (g1_idx in 1:(length(grouplist) - 1)) {
  for (g2_idx in (g1_idx + 1):length(grouplist)) {
    g1 <- grouplist[g1_idx]
    g2 <- grouplist[g2_idx]
    CI_analysis(g1 = g1,g2 = g2)
    
  }
}

# Function to search through every folder in a DIR, find a given file, then return the average UVC or UVPM for that group

MakeAverage <- function(location = "/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/", est = "UVC", group = 'fcn')
{
  # Load the files
  setwd(location)
  # Get a list of all the folders in the directory
  folders <- list.dirs(location, full.names = FALSE, recursive = FALSE)
  # Initialize a list to store the data
  data_list <- list()
  # Loop through each folder
  for(folder in folders)
  {
    # Get the file name
    filename <- paste0("ADNI_OD_", group, "_mean_2e+05_1000.rdata")
    # Check if the file exists in the folder
    if(file.exists(paste0(folder, "/", filename)))
    {
      # Load the file
      data <- readRDS(paste0(folder, "/", filename))
      # Extract the model
      model <- data$model
      # Extract the UVC or UVPM
      if(est == "UVC")
      {
        data_list[[folder]] <- model$UVC
      }
      else if(est == "UVPM")
      {
        data_list[[folder]] <- model$UVPM
      }
    }
  }

  #make the average matrix, and save it as a file in location
  if(length(data_list) > 0)
  {
    # Get the dimensions of the matrix
    dim1 <- dim(data_list[[1]])
    # Initialize a matrix to store the average
    average_matrix <- matrix(0, nrow = dim1[1], ncol = dim1[2])
    # Loop through each matrix and add it to the average
    for(matrix in data_list)
    {
      average_matrix <- average_matrix + matrix
    }
    # Divide by the number of matrices to get the average
    average_matrix <- average_matrix / length(data_list)
    # Save the average matrix as an RDS file
    saveRDS(average_matrix, paste0(location, "Average_", group, "_", est, ".rds"))
  }
}

CI_analysis(av = TRUE)





replist <- c(1:10, 42)
group <- 'fcn'

for (g in grouplist){
  print("---------------------")
  print(g)
  print("---------------------")
  for (r in replist){
  #set wd to the head folder
    setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD")
  #build the folder and file name 
    modelloc <- paste0("/N/u/conlcorn/BigRed200/SexLinkedProject/output/FinalFiles/OD/Rep-",r)
    filename <- paste0("ADNI_OD_",g,"_mean_2e+05_1000.rdata")
  #call the correlation function 
    Correlation(modelname = filename,modelloc = modelloc,rep = r)
  }
}

