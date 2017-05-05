library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(xlsx)
library(plater)
library(Hmisc)

#Raw format function for output of SoftMax plate-reader files.
plate_input_96 <- function(filename, ...){
  data1 <- read.table(file = filename, 
                      skip = 2,
                      sep = "\t", 
                      quote = "\"", 
                      dec = ".", 
                      fill = TRUE, 
                      comment.char = "",
                      nrows = 8,
                      colClasses = c(NULL, rep("numeric", 12), NULL))
  data2 <- data1[, c(2:13)]
  data2 <- data.frame(LETTERS[1:8], data2)#making a column of letters 
  #absorbance is used to input for the plate_read() function, used later
  names(data2) <- c("Absorbance", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  return(data2)
}

tidy_plate <- function(plate_input){
  data3 <- melt(plate_input)#making dataframe tidy
  Wells <- paste(data3$Absorbance, data3$variable, sep = "")#putting the letters and numbers together
  finale <- data.frame(Wells, Absorbance = data3$value)#combining them well names and absorbance vals
  finale2 <- finale[order(finale$Wells), ]#reordering so the columns are compatible with later dfs
  return(finale2)
}

#loading the transformed SoftMax file, and the sample map file (which contains plate info)
plate_load <- function(sample_map, softmax_file){
  reads_raw <- plate_input_96(filename = softmax_file)
  reads_tidy <- tidy_plate(reads_raw)
  
  sample_tidy <- read_plate(sample_map)#reading file using plater function
  subset_df <- reads_tidy[reads_tidy$Wells %in% sample_tidy$Wells, ]#only the absorbances that have treatment labels
  tidied_df <- merge(subset_df, sample_tidy)
  return(tidied_df)
}
