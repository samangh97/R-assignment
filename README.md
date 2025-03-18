# README
# Read the files and store them in variables
fang_data <- read.table("fang_et_al_genotypes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
snp_data <- read.table("snp_position.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Inspect dimensions of the data
dim(fang_data)  # Dimensions of the fang_et_al_genotypes data
dim(snp_data)   # Dimensions of the snp_position data

# Transpose the fang_et_al_genotypes data
transposed_genotypes <- t(fang_data)

transposed_genotypes <- as.data.frame(transposed_genotypes)
write.table(transposed_genotypes, file = "transposed_genotypes.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

# Read the file
snp_lines <- readLines("snp_position.txt")

# Extract the header and add two blank lines
header <- snp_lines[1]
modified_lines <- c(header, "", "", snp_lines[-1])

# Write the modified content to a new file
writeLines(modified_lines, "snp_position_mod.txt")
# Read the input file
transposed_genotypes <- read.table("transposed_genotypes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# Identify columns to keep based on the 3rd row
cols_to_keep <- which(transposed_genotypes[3, ] %in% c("ZMMLR", "ZMMIL", "ZMMMR"))
# Subset the data to keep only the marked columns
maize_data <- transposed_genotypes[, cols_to_keep]

# Save the subsetted data to a new file
write.table(maize_data, file = "maize.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the modified file
maize_check <- read.table("maize.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# Read snp_position_mod.txt
snp_position <- read.table("snp_position_mod.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Read maize.txt
maize_data <- read.table("maize.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# Extract columns 1, 3, and 4
snp_subset <- snp_position[, c(1, 3, 4)]
# Remove the first two rows from maize_data
maize_data <- maize_data[-(1:2), ]
# Combine the subsetted SNP data with maize data
maize_joined <- cbind(snp_subset, maize_data)
# Save the combined data to a new file
write.table(maize_joined, file = "maize_joined.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Read the modified file
maize_joined_check <- read.table("maize_joined.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# Load required package
library(dplyr)

# Read the data
maize_data <- read.table("maize_joined.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Function to check if a value is numeric
is_numeric_str <- function(x) {
  !is.na(suppressWarnings(as.numeric(x)))
}

# Function to process and save each chromosome
process_chromosome <- function(chrom_num, data) {
  subset_data <- data[data[[2]] == chrom_num, ]  # Filter by chromosome number
  
  # Identify numeric and non-numeric SNP positions
  numeric_positions <- sapply(subset_data[[3]], is_numeric_str)
  
  # Convert only numeric values to numbers for sorting
  subset_data$SNP_numeric <- as.numeric(subset_data[[3]])
  
  # Sort by SNP position (numeric first, "?" last)
  subset_data <- subset_data[order(subset_data$SNP_numeric, na.last = TRUE), ]
  
  # Restore original SNP column with "?" for missing values
  subset_data[[3]][!numeric_positions] <- "?"
  
  # Remove the extra SNP_numeric column after sorting
  subset_data$SNP_numeric <- NULL
  
  # Define file name
  output_file <- paste0("maize_sortedv1_", chrom_num, ".txt")
  
  # Write to file
  write.table(subset_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("File successfully written: ", output_file)
}

# Apply function to chromosomes 1 to 10
lapply(1:10, process_chromosome, data = maize_data)


# Read the data
maize_data <- read.table("maize_joined.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Function to check if a value is numeric
is_numeric_str <- function(x) {
  !is.na(suppressWarnings(as.numeric(x)))
}

# Function to process and save each chromosome (sorted in decreasing order)
process_chromosome_desc <- function(chrom_num, data) {
  subset_data <- data[data[[2]] == chrom_num, ]  # Filter by chromosome number
  
  # Identify numeric and non-numeric SNP positions
  numeric_positions <- sapply(subset_data[[3]], is_numeric_str)
  
  # Convert only numeric values to numbers for sorting
  subset_data$SNP_numeric <- as.numeric(subset_data[[3]])
  
  # Sort by SNP position in decreasing order (numeric first, "-" last)
  subset_data <- subset_data[order(-subset_data$SNP_numeric, na.last = TRUE), ]
  
  # Restore original SNP column with "-" for missing values
  subset_data[[3]][!numeric_positions] <- "-"

  # Remove the extra SNP_numeric column after sorting
  subset_data$SNP_numeric <- NULL
  
  # Define file name
  output_file <- paste0("maize_sortedv2_", chrom_num, ".txt")
  
  # Write to file
  write.table(subset_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("File successfully written: ", output_file)
}

# Apply function to chromosomes 1 to 10
lapply(1:10, process_chromosome_desc, data = maize_data)
#repeating the process for teosinte data
# Read the files and store them in variables
fang_data <- read.table("fang_et_al_genotypes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
snp_data <- read.table("snp_position.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Transpose the fang_et_al_genotypes data
print("Transposing fang data...")
transposed_genotypes <- t(fang_data)
transposed_genotypes <- as.data.frame(transposed_genotypes)
write.table(transposed_genotypes, file = "transposed_genotypes.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)


# Read and modify snp_position file
print("Modifying snp_position file...")
snp_lines <- readLines("snp_position.txt")
header <- snp_lines[1]
modified_lines <- c(header, "", "", snp_lines[-1])
writeLines(modified_lines, "snp_position_mod.txt")



# Read the transposed data
print("Reading transposed data...")
transposed_genotypes <- read.table("transposed_genotypes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Identify columns to keep (teosinte categories: ZMPBA, ZMPIL, ZMPJA)
print("Identifying columns to keep...")
cols_to_keep <- which(transposed_genotypes[3, ] %in% c("ZMPBA", "ZMPIL", "ZMPJA"))
teosinte_data <- transposed_genotypes[, cols_to_keep]


# Save the filtered teosinte data
print("Saving teosinte data...")
write.table(teosinte_data, file = "teosinte.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Read files
print("Reading teosinte and SNP data...")
teosinte_check <- read.table("teosinte.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
snp_position <- read.table("snp_position_mod.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)


# Extract relevant SNP columns
print("Extracting SNP columns...")
snp_subset <- snp_position[, c(1, 3, 4)]
teosinte_data <- read.table("teosinte.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
teosinte_data <- teosinte_data[-(1:2), ]  # Remove first two rows


# Merge SNP data with teosinte genotypes
print("Merging SNP and teosinte data...")
teosinte_joined <- cbind(snp_subset, teosinte_data)
write.table(teosinte_joined, file = "teosinte_joined.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Load required package
library(dplyr)


# Read the joined data
print("Reading joined data...")
teosinte_data <- read.table("teosinte_joined.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)



# Function to check if a value is numeric
is_numeric_str <- function(x) {
  !is.na(suppressWarnings(as.numeric(x)))
}

# Function to process and save each chromosome (sorted in increasing order)
process_chromosome <- function(chrom_num, data) {
  print(paste("Processing chromosome", chrom_num, "..."))
  subset_data <- data[data[[2]] == chrom_num, ]  
  numeric_positions <- sapply(subset_data[[3]], is_numeric_str)
  subset_data$SNP_numeric <- as.numeric(subset_data[[3]])
  subset_data <- subset_data[order(subset_data$SNP_numeric, na.last = TRUE), ]
  subset_data[[3]][!numeric_positions] <- "?"
  subset_data$SNP_numeric <- NULL  
  output_file <- paste0("teosinte_sortedv1_", chrom_num, ".txt")
  write.table(subset_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("File successfully written:", output_file))
}

# Apply function to chromosomes 1 to 10
print("Processing chromosomes 1 to 10 (increasing order)...")
lapply(1:10, process_chromosome, data = teosinte_data)

# Function to process and save each chromosome (sorted in decreasing order)
process_chromosome_desc <- function(chrom_num, data) {
  print(paste("Processing chromosome", chrom_num, "(decreasing order)..."))
  subset_data <- data[data[[2]] == chrom_num, ]  
  numeric_positions <- sapply(subset_data[[3]], is_numeric_str)
  subset_data$SNP_numeric <- as.numeric(subset_data[[3]])
  subset_data <- subset_data[order(-subset_data$SNP_numeric, na.last = TRUE), ]
  subset_data[[3]][!numeric_positions] <- "-"
  subset_data$SNP_numeric <- NULL  
  output_file <- paste0("teosinte_sortedv2_", chrom_num, ".txt")
  write.table(subset_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("File successfully written:", output_file))
}

# Apply function to chromosomes 1 to 10
print("Processing chromosomes 1 to 10 (decreasing order)...")
lapply(1:10, process_chromosome_desc, data = teosinte_data)

#visualization 
# Load libraries
library(ggplot2)
library(dplyr)

# Step 1: Read the data
maize_data <- read.table("maize_joined.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
maize_data <- maize_data[, 2:3]
colnames(maize_data) <- c("Chromosome", "SNP_Position")

teosinte_data <- read.table("teosinte_joined.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
teosinte_data <- teosinte_data[, 2:3]
colnames(teosinte_data) <- c("Chromosome", "SNP_Position")

# Step 2: Filter out invalid SNP positions
filter_data <- function(data) {
    data %>%
        filter(!is.na(SNP_Position) & grepl("^[0-9]+$", SNP_Position)) %>%
        mutate(SNP_Position = as.numeric(SNP_Position))
}

maize_data <- filter_data(maize_data)
teosinte_data <- filter_data(teosinte_data)

# Step 3: Combine the data
maize_data$Species <- "Maize"
teosinte_data$Species <- "Teosinte"
combined_data <- bind_rows(maize_data, teosinte_data)

# Ensure Chromosome is a factor with the correct order
chromosome_order <- c(1:10, "unknown", "multiple")
combined_data$Chromosome <- factor(combined_data$Chromosome, levels = chromosome_order)

# Step 4: Bar plot with error bars
snp_counts <- combined_data %>%
    group_by(Chromosome, Species) %>%
    summarise(SNP_Count = n(), .groups = "drop") %>%
    mutate(SE = sqrt(SNP_Count))  # Standard error for error bars

bar_plot <- ggplot(snp_counts, aes(x = Chromosome, y = SNP_Count, fill = Species)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +  # Bar plot with dodge
    geom_errorbar(aes(ymin = SNP_Count - SE, ymax = SNP_Count + SE),  # Error bars
                position = position_dodge(width = 0.9), width = 0.25) +
    scale_fill_manual(values = c("Maize" = "blue", "Teosinte" = "pink")) +  # Custom colors
    labs(title = "SNP Counts per Chromosome",
         x = "Chromosome",
         y = "SNP Count",
         fill = "Species") +
    theme_minimal()

print(bar_plot)

# Step 5: Dot plot for SNP positions
dot_plot <- ggplot(combined_data, aes(x = Chromosome, y = SNP_Position, color = Species)) +
    geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.6) +  # Dots with jitter
    scale_color_manual(values = c("Maize" = "blue", "Teosinte" = "pink")) +  # Custom colors
    labs(title = "SNP Positions by Chromosome",
         x = "Chromosome",
         y = "SNP Position",
         color = "Species") +
    theme_minimal()

print(dot_plot)
#heterozygosity and homozygosity
# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Read the data
maize_data <- read.table("maize_joined.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
teosinte_data <- read.table("teosinte_joined.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Step 2: Define a function to classify sites
classify_sites <- function(data) {
    # Select columns with nucleotide data (from column 4 to the last column)
    nucleotide_data <- data[, 4:ncol(data)]
    
    # Classify each site as homozygous, heterozygous, or missing
    result <- apply(nucleotide_data, 2, function(sample) {
        # Count homozygous, heterozygous, and missing sites
        homozygous <- sum(grepl("^(A/A|C/C|G/G|T/T)$", sample))  # Count homozygous sites
        heterozygous <- sum(grepl("^(A/C|A/G|A/T|C/A|C/G|C/T|G/A|G/C|G/T|T/A|T/C|T/G)$", sample))  # Count heterozygous sites
        missing <- sum(grepl("^(\\?/\\?|\\./\\.|./.)$", sample))  # Count missing sites
        
        # Calculate total nucleotides (excluding missing data)
        total_nucleotides <- sum(homozygous, heterozygous)
        
        # Calculate proportions (normalized to total nucleotides)
        data.frame(
            Homozygous = homozygous / total_nucleotides,
            Heterozygous = heterozygous / total_nucleotides,
            Missing = missing / total_nucleotides  # Missing data normalized to total sites
        )
    })
    
    # Combine results into a data frame
    result <- do.call(rbind, result)
    result$Sample <- colnames(nucleotide_data)  # Add sample names
    return(result)
}

# Step 3: Classify sites for maize and teosinte
maize_classified <- classify_sites(maize_data)
teosinte_classified <- classify_sites(teosinte_data)

# Add a "Species" column
maize_classified$Species <- "Maize"
teosinte_classified$Species <- "Teosinte"

# Combine the data
combined_classified <- bind_rows(maize_classified, teosinte_classified)

# Step 4: Calculate average percentages for maize and teosinte
average_percentages <- combined_classified %>%
    group_by(Species) %>%
    summarise(
        Avg_Homozygous = mean(Homozygous) * 100,
        Avg_Heterozygous = mean(Heterozygous) * 100,
        Avg_Missing = mean(Missing) * 100,
        .groups = "drop"
    )

# Print the average percentages
print("Average percentages of homozygous, heterozygous, and missing sites:")
print(average_percentages)

# Step 5: Reshape the data for plotting
combined_long <- combined_classified %>%
    pivot_longer(cols = c(Homozygous, Heterozygous, Missing), names_to = "Type", values_to = "Proportion")

# Step 6: Visualize the results
plot <- ggplot(combined_long, aes(x = Sample, y = Proportion, fill = Type)) +
    geom_bar(stat = "identity", position = "fill") +  # Stacked bar plot with normalized height
    facet_wrap(~ Species, scales = "free_x") +  # Separate panels for maize and teosinte
    scale_fill_manual(values = c("Homozygous" = "blue", "Heterozygous" = "pink", "Missing" = "gray")) +  # Custom colors
    labs(title = "Proportion of Homozygous, Heterozygous, and Missing Sites",
         x = "Sample",
         y = "Proportion",
         fill = "Site Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels

print(plot)  

#homozygous and heterozygous
# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)

# Read the genotype file
genotype_data <- read.table("fang_et_al_genotypes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Split into maize and teosinte based on column 3
maize_data <- genotype_data[genotype_data[[3]] %in% c("ZMMIL", "ZMMLR", "ZMMMR"), ]
teosinte_data <- genotype_data[genotype_data[[3]] %in% c("ZMPBA", "ZMPIL", "ZMPJA"), ]

# Function to calculate missing, homozygous, and heterozygous proportions
calculate_counts <- function(df, group_name) {
  # Select genotype columns (4th column onward)
  genotype_matrix <- df[, 4:ncol(df)]
  
  # Initialize result dataframe
  result <- data.frame(Missing = integer(nrow(df)),
                       Homozygous = integer(nrow(df)),
                       Heterozygous = integer(nrow(df)),
                       Group = group_name)  # Add a group identifier
  
  # Loop through each row
  for (i in 1:nrow(genotype_matrix)) {
    row_values <- genotype_matrix[i, ]
    
    # Count missing, homozygous, and heterozygous values
    missing_count <- sum(row_values == "?/?")
    homozygous_count <- sum(row_values %in% c("A/A", "T/T", "G/G", "C/C"))
    heterozygous_count <- sum(!(row_values %in% c("?/?", "A/A", "T/T", "G/G", "C/C")))
    
    # Normalize counts
    total <- missing_count + homozygous_count + heterozygous_count
    result[i, 1:3] <- c(missing_count / total, homozygous_count / total, heterozygous_count / total)
  }
  
  return(result)
}

# Apply function to maize and teosinte data
maize_counts <- calculate_counts(maize_data, "Maize")
teosinte_counts <- calculate_counts(teosinte_data, "Teosinte")

# Combine data for plotting
combined_counts <- rbind(maize_counts, teosinte_counts)

# Convert data to long format for ggplot
long_data <- pivot_longer(combined_counts, cols = c("Missing", "Homozygous", "Heterozygous"),
                          names_to = "Category", values_to = "Proportion")

# Generate normalized stacked bar plot
ggplot(long_data, aes(x = Group, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +  # Normalized height
  labs(title = "Normalized Genotype Proportions in Maize and Teosinte",
       x = "Group", y = "Proportion") +
  scale_fill_manual(values = c("Missing" = "gray", "Homozygous" = "blue", "Heterozygous" = "red")) +
  theme_minimal()


