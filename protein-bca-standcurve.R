library(RCurl)

script <- getURL("https://raw.githubusercontent.com/bertrand-lab/Softmax-Platereader-Functions/master/Softmax%20Platereader%20Functions.R", 
                 ssl.verifypeer = FALSE)
eval(parse(text = script))

# open project in working directory (or set working directory)

# edit out top and bottom lines from txt file

# create sample map

# Run

finale_wrapper(sample_map = "", 
               softmax_file = "", 
               dilution_factor = 20)