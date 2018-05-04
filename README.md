# Softmax-Platereader-Functions
Functions for reading in data from the SoftMax software for a 96 well plate.

How to:

1) Set working directory
2) Edit out top and bottom lines from txt file from softmax pro output (see example file)
3) Create sample map (see example)
4) Load all functions.
5) Run:

`finale_wrapper(sample_map = "sample-map.csv", "softmax-output.txt", dilution_factor = 20)`

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
