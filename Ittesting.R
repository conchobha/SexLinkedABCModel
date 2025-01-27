# Source the Functions needed to run the model
library(MASS)
library(foreach)
library(doParallel)
library(parallel)
library(readr) # Used for reading CSV file faster than the built in r function
library(dplyr) # used for selecting collums in the read functions
library(caTools) #Used for splitting data 
library(doParallel)




setwd("/N/u/conlcorn/BigRed200/SexLinkedProject/methods")
source("dMRIABC.R") # Adds the Functions we need in order to run the model
#Set up The variables
dn <- 1:8 
sn <- 10000
fn <- 100
burns <- 1000
rep <- 1:10
grouplist <- list('fcn', 'fmci', 'mcn', 'mmci', 'combined') # Ensure correct names for groups
datametriclist <- list(
  'FA', 'OD', 'meanlength', 'numoffibers' # Simplified for testing
) # For the testing, we will only use combined and FA, since this is only for itteration testing 
#Set up the List of Iterations
# We know 500k is too much, so we will do those under it 
itts <-seq(100000,400000,by=50000) #builds a list from 100k - 400k in increments of 50000
# This creates seven runs to do, which we can parallel 
#Build the Parallel 
# We don't care about saving the model, but we want to report the final time to compute for each model
# I edited the saving of the model to save how long it takes. This should help

#registerDoParallel(7) # We only need 7 for this
#results <-foreach(i = itts) %dopar%{
args = commandArgs(trailingOnly = TRUE)
k <-as.numeric(args[1])


RunModel(sn=itts[k])
  #}
