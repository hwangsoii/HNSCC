# Load required libraries
library(tools)

# Define the path
msisensor_pro_path <- '/Users/hwangso_macmini/Google Drive/My Drive/Work/project/HNSCC_2/TRIUMPH_result/msisensor'

# Initialize file list
file_list <- list.files(path = msisensor_pro_path, full.names = TRUE)

# Sort the file list
MSIsensor_list <- sort(file_list)

# Commented out part of Python code for writing to a new sample
# new_file <- paste0(msisensor_pro_path, "/msi_count.txt")
# write.table(data.frame(Sample = character(), Total_msi_locus = character()),
#             new_file, sep = '\t', row.names = FALSE, col.names = TRUE)


# Define the function to extract value from 2nd line and 3rd column
read_value_from_file <- function(file_path) {
  # Read the 2nd line of the file
  second_line <- readLines(file_path, n = 2)[2]
  
  # Split the line on the tab character to get the columns
  columns <- unlist(strsplit(second_line, "\t"))
  
  # Return the third column's value
  return(columns[3])
}

# Initialize dataframe to store results
result_df <- data.frame(Sample = character(), MSI = character(), stringsAsFactors = FALSE)
# Loop over sorted list of MSIsensor files
for (file in MSIsensor_list) {
  if (grepl('_dis$', file)) next
  if (grepl('_germline$', file)) next
  if (grepl('_somatic$', file)) next
  
  # print(file)
  
  # Extract the sample name from the file path
  sample <- file_path_sans_ext(basename(file))
  msi <- read_value_from_file(file)
  # print(c(sample, msi))  
  # Append to the dataframe
  result_df <- rbind(result_df, data.frame(Sample2 = sample, MSI = msi, stringsAsFactors = FALSE))
}


# ... [The rest of your code remains unchanged]

# Create a new column 'Status' based on MSI values
# result_df$MS <- ifelse(as.numeric(result_df$MSI) >= 10, "MSI-H", ifelse(as.numeric(result_df$MSI) >= 3.5), "MSI-L", "MSS")
result_df$MS <- ifelse(as.numeric(result_df$MSI) >= 3.5, "MSI", "MSS")

# Print the result dataframe
# print(result_df)

