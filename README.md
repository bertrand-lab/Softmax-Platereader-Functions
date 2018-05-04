# Softmax-Platereader-Functions
Functions for reading in data from the SoftMax software for a 96 well plate.

How to:

1) Open R script and copy from `protein-bca-standcurve.R`:

`library(RCurl)`

# Upload the latest version of softmax reader functions.
`script <- getURL("https://raw.githubusercontent.com/bertrand-lab/Softmax-Platereader-Functions/master/Softmax%20Platereader%20Functions.R", ssl.verifypeer = FALSE)`

# Run the lastest version
`eval(parse(text = script))`

# open project in working directory (or set working directory)

# edit out top and bottom lines from txt file

# create sample map

# Run

`finale_wrapper(sample_map = "", `
               `softmax_file = "", `
               `dilution_factor = 20)`

To Do:

- Figure out how to *not* have to edit the top and bottom lines from softmax pro output.
- Improve functions defensibility. 
- Put fitted model into plot output.

Function Descriptions:

plate_input_96()
 - Reads in data from a .txt file exported from SoftMax
 - Converts the columns and row names to the appropriate headings (A, B, C, etc.) and (1-8).
 
tidy_plate()
  - Converts the above object into a tidy dataframe, where observations correspond to a well ID and an absorbance value
  - Specifically written for working with the package {plater}
  
plate_load()
  - Wrapper function including the two above.
  - Also requires another file, the "sample map". This conforms to the format fromt {plater}, where the treatments applied to each well are put in a .csv file representing a 96 well plate. 
  - The read-in absorbance data is then merged with the sample map, all in a tidy data frame.
