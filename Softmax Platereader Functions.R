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
  #absorbance is used to input for tOhe plate_read() function, used later
  names(data2) <- c("Absorbance", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  return(data2)
}

tidy_plate <- function(plate_input){
  data3 <- melt(plate_input)# making dataframe tidy
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


setwd("C:\\Users\\Scott\\Google Drive\\Projects\\Helping Data Files\\Elden")

prot_conc <- plate_load(sample_map = "bca-sample-map.csv", 
                        softmax_file = "ER180502_BCA1.txt")

# blank subtraction
blank_sub_val <- mean(prot_conc[prot_conc$standards == 0,]$Absorbance, na.rm = TRUE)

prot_conc$zero_abs <- prot_conc$Absorbance - blank_sub_val


# fit model and predict proteins

lm_out1 <- lm(formula = standards ~ zero_abs + I(zero_abs^2), data = prot_conc)

new_df <- data.frame(zero_abs = prot_conc[prot_conc$is_standard == FALSE,]$zero_abs)

predicted_concs <- predict(object = lm_out1, 
                           newdata = new_df)

prot_conc$predicted_prot <- NULL

prot_conc[prot_conc$is_standard == FALSE,"predicted_prot"] <- as.vector(predicted_concs)

# outside detection limit

prot_conc2 <- prot_conc %>% 
  filter(is_standard == FALSE) %>% 
  dplyr::group_by(sample_id) %>%   
  summarise(sample_mean = mean(predicted_prot, na.rm = ),
            sample_median = median(predicted_prot),
            sample_sd = sd(predicted_prot),
            sample_sdmean = sample_sd/sample_mean)

# warning column
prot_conc2$warning <- ifelse(test = prot_conc2$sample_sdmean > 0.2, yes = 'warning', no = 'good-to-go')

#detection limit warning
prot_conc2$detection_limit <- ifelse(test = prot_conc2$sample_median < 1, yes = 'below-ld', no = 'gtg')

dilution_factor <- 37.5

prot_conc2$transformed_prot <- prot_conc2$sample_median*dilution_factor

prot_conc2[prot_conc2$detection_limit == 'below-ld',]$transformed_prot <- NULL

# FIGURE
# WRITE A CSV FILE


#this function takes only one plate at a time. it assumes that your Standards column is called "Standards"
blanked_df <- function(plate_load_output) {
  
  blank_abs <- plate_load_output[plate_load_output$Standards == 0, ]$Absorbance %>% mean(na.rm = TRUE)
  treatment_df <- plate_load_output %>% filter(Treatment != "standard")
  blank_subtracted_df <- treatment_df %>% mutate(abs_corrected = Absorbance - blank_abs)
  return(blank_subtracted_df)
  
  }

#blanked_df(plate_load_output = standards_treat)

#standard_curve_creation





