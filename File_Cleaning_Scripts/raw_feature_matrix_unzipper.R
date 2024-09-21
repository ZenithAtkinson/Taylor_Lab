setwd("~/Lab_Data/Cloud_Analysis_Output_Files/")
# Current directory: ~/Lab_Data/Cloud_Analysis_Output_Files
tar_files <- list.files(pattern = "\\.tar\\.gz$")
print(tar_files)
# Location of output: ~/Lab_Data/(given folder)

# Function to create a folder name from a tar.gz file, excluding "raw_feature_bc_matrix"
create_folder_name <- function(tar_file) {
  # Remove the .tar.gz suffix
  folder_name <- gsub("\\.tar\\.gz$", "", tar_file)
  
  # Remove "raw_feature_bc_matrix" from the folder name
  folder_name <- gsub("_raw_feature_bc_matrix", "", folder_name)
  
  return(folder_name)
}

# Iterate over the tar.gz files
for (tar_file in tar_files) {
  # Create the folder name without "raw_feature_bc_matrix"
  folder_name <- create_folder_name(tar_file)
  
  # Define the path for the new folder (one level up)
  folder_path <- file.path("..", folder_name)
  
  # Create the main folder if it doesn't exist
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }
  
  # Define the subfolder path (raw_feature_bc_matrix)
  subfolder_path <- file.path(folder_path, "raw_feature_bc_matrix")
  
  # Create the subfolder inside the main folder
  if (!dir.exists(subfolder_path)) {
    dir.create(subfolder_path)
  }
  
  # Extract the tar.gz file into the subfolder
  untar(tar_file, exdir = subfolder_path)
  
  # Print a message to confirm
  cat("Extracted", tar_file, "into", folder_path, "\n")
}
