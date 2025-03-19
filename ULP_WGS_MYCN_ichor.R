# Load required package
library(dplyr)

# Set the parent directory path
parent_dir <- "./ichor_cna_seg_files"

# Initialize an empty dataframe to store results
results <- data.frame(file_name = character(), stringsAsFactors = FALSE)

# Get a list of all subdirectories
subdirs <- list.dirs(parent_dir, recursive = TRUE, full.names = TRUE)

# Loop over each subdirectory
for (subdir in subdirs) {
  # Get a list of files in the subdirectory
  files <- list.files(subdir, pattern = "\\.cna\\.seg$", full.names = TRUE)
  
  # Loop over each file in the subdirectory
  for (file in files) {
    # Read the file into a dataframe
    df <- tryCatch(read.table(file, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    
    # If file is successfully read
    if (!is.null(df)) {
      # Filter the relevant row
      df <- df %>%
        mutate(V1 = as.numeric(V1), V2 = as.numeric(V2), V3 = as.numeric(V3))
      filtered_row <- df %>%
        filter(V1 == 2, V2 == 15000001, V3 == 16000000)
      
      # If the row exists, extract the required data
      if (nrow(filtered_row) > 0) {
        # Extract the file name and modify it
        file_name <- basename(file) %>%
          gsub("RP-2413_", "", .) %>%
          gsub("_v1_WGS_OnPrem\\.cna\\.seg", "", .)
        
        # Combine the file name and the entire filtered row
        filtered_row <- cbind(file_name = file_name, filtered_row)
        
        # Append the data to the results dataframe
        results <- rbind(results, filtered_row)
      }
    }
  }
}

# Assign proper column names to the results dataframe
colnames(results)[1] <- "file_name"  # Ensure the first column has the correct name

# Write the results to a file or print
write.table(results, file = "output_file.txt", sep = "\t", row.names = FALSE, quote = FALSE)
