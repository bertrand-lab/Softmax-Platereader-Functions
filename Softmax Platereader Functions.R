library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(xlsx)
library(plater)
library(Hmisc)
library(gridExtra)
library(ggforce)

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


# setwd("C:\\Users\\Scott\\Google Drive\\Projects\\Helping Data Files\\Elden")

# prot_conc <- plate_load(sample_map = "bca-sample-map.csv", 
                        # softmax_file = "ER180502_BCA1.txt")

# blank subtraction
blank_subtraction <- function(tidy_prot_file){
  
  blank_sub_val <- mean(tidy_prot_file[tidy_prot_file$standards == 0,]$Absorbance, na.rm = TRUE)
  tidy_prot_file$zero_abs <- tidy_prot_file$Absorbance - blank_sub_val
  return(tidy_prot_file)
  
}
# blank_prot_conc <- blank_subtraction(prot_conc)


# fit model and predict proteins
fit_and_predict <- function(blanked_prot_conc){
  
  lm_out1 <- lm(formula = standards ~ zero_abs + I(zero_abs^2), data = blanked_prot_conc)
  
  new_df <- data.frame(zero_abs = blanked_prot_conc[blanked_prot_conc$is_standard == FALSE,]$zero_abs)
  
  predicted_concs <- predict(object = lm_out1, 
                             newdata = new_df)
  
  blanked_prot_conc$predicted_prot <- NULL
  
  blanked_prot_conc[blanked_prot_conc$is_standard == FALSE,"predicted_prot"] <- as.vector(predicted_concs)
  
  return(blanked_prot_conc)
}
# fitted_conc <- fit_and_predict(blanked_prot_conc = blank_prot_conc)


# outside detection limit
detection_limit <- function(fitted_conc, dilution_factor){
  prot_conc2 <- fitted_conc %>% 
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
  
  # dilution_factor <- 37.5
  
  prot_conc2$transformed_prot <- prot_conc2$sample_median*dilution_factor
  
  prot_conc2[prot_conc2$detection_limit == 'below-ld',]$transformed_prot <- NULL
  
  return(prot_conc2)
}


# detect_test <- detection_limit(fitted_conc = fitted_conc, 20)

finale_wrapper <- function(sample_map, softmax_file, dilution_factor, 
                           gen_plot = TRUE, save_plot = TRUE, make_csv = TRUE){
  
  loaded_plate <- plate_load(sample_map = sample_map, softmax_file = softmax_file)
  blank_prot_conc <- blank_subtraction(loaded_plate)
  fitted_conc <- fit_and_predict(blanked_prot_conc = blank_prot_conc)
  outside_detect_conc <- detection_limit(fitted_conc = fitted_conc, dilution_factor = dilution_factor)
  
  if(gen_plot){
    if(save_plot){
      if(!is.null(dev.list())){
        dev.off()
      }
      png(paste(softmax_file, "_dilution", dilution_factor, ".png", sep = ""), 
          width=23*0.75, height=27.94*0.75, units="cm", res=700)
    }
    st_curve <- blank_prot_conc %>% 
      ggplot(aes(x = standards, y = zero_abs)) + 
      labs(y = 'Absorbance (zeroed)', x = 'Protein Concentration (Standards)') +
      geom_point(shape = 18, size = 3, alpha = 0.8) + 
      theme_bw()  + facet_zoom(x = standards < 25, y = zero_abs < 0.25)
    
    conc_plot <- fitted_conc %>% 
      filter(is_standard == FALSE) %>% 
      ggplot() + geom_point(aes(x = sample_id %>% as.factor(), 
                            y = predicted_prot, 
                            colour = sample_id %>% as.factor()), 
                          size = 3, alpha = 0.9) +
      labs(x = 'Sample ID', y = 'Protein Concentration (ug/ml)') +
      geom_point(data = outside_detect_conc, aes(x = sample_id, y = sample_mean), size = 2, alpha = 0.4) +
      geom_point(data = outside_detect_conc, aes(x = sample_id, y = sample_median), size = 3, 
                 shape = 18, alpha = 0.8) + 
      scale_colour_discrete(guide = FALSE) +
      theme_bw()
    
    calculate_prot <- outside_detect_conc %>% 
      ggplot(aes(x = sample_id %>% as.factor(), y = transformed_prot)) + 
      geom_point() + scale_colour_discrete(guide = FALSE) + 
      labs(x = 'Sample ID', y = 'Dilution Factor Transformed \nProtein Concentration (ug/ml)') +
      theme_bw()
    
    grid.arrange(st_curve, conc_plot, calculate_prot)
    if(save_plot){
      dev.off()
      cur_dir <- getwd()
      print(paste0('Plot saved in ', cur_dir))
    }
  }
  if(make_csv){
    write.csv(outside_detect_conc, file = paste(softmax_file, "_dilution", dilution_factor, ".csv", sep = ""))
  }
  
  return(outside_detect_conc)
}



# finale_wrapper(sample_map = "bca-sample-map.csv", "ER180502_BCA1.txt", dilution_factor = 20)






