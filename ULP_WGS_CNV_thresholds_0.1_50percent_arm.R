
# Read chromosome arms information
x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
            col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
x <- x[ , .(length = sum(chromEnd - chromStart)), by = .(chrom, arm = substring(name, 1, 1)) ]
chr_arms <- as.data.frame(x)
chr_arms$chrom <- sub('chr','',chr_arms$chrom)
#chr_arms$chrom <- as.numeric(chr_arms$chrom)
colnames(chr_arms) <- c('chromosome','arm','length')

# Define centromere locations
centromere <- chr_arms[chr_arms$arm == 'p',]
centromere <- centromere[,c(1,3)]
colnames(centromere) <- c('chromosome','end')


# Define arms of interest
gain_arms <- c("1q", "2p", "17q")
loss_arms <- c("1p", "3p", "4p", "11q")
all_arms <- c(gain_arms, loss_arms)



# Function to process a single seg file
process_seg_file <- function(seg_file) {
  # Read seg file
  seg_data <- fread(seg_file)
  setDT(seg_data)
  
  seg_data <- seg_data[,c(1,2,3,4,10,5)]
  colnames(seg_data) <- c("Sample", "Chromosome", "Start", "End", "CN", "Num_Markers")

  # Convert CN to log2(copy ratio)
  seg_data[, log2_CR := log2(CN / 2)]

  # Get sample name
  sample_name <- unique(seg_data$Sample)
  sample_name <- sub('RP-2413_','',sample_name)
  sample_name <- sub('_v1_WGS_OnPrem','',sample_name)

  # Initialize results
  sample_result <- setNames(rep("no", length(all_arms)), all_arms)

  # Loop over only the arms of interest
  for (chr_arm in all_arms) {
    chr <- as.integer(gsub("[pq]", "", chr_arm))
    arm <- gsub("[0-9]", "", chr_arm)

    # Get centromere position
    centromere_pos <- centromere[centromere$chromosome == chr, 2]

    # Get arm length and 50% threshold
    arm_length <- chr_arms[chr_arms$chromosome == chr & chr_arms$arm == arm, 3]
    half_arm_length <- arm_length * 0.5

    # Subset segments for this chromosome
    chr_segments <- seg_data[Chromosome == chr]

    # Initialize list for split segments
    split_segments <- list()

    # Loop through each segment and adjust it
    for (i in 1:nrow(chr_segments)) {
      row <- chr_segments[i]
      start <- row$Start
      end <- row$End
      log2_CR <- row$log2_CR

      if (end <= centromere_pos) {
        # Segment fully in p-arm
        split_segments[[length(split_segments) + 1]] <- data.table(
          Chromosome = row$Chromosome, Start_Adjusted = start, End_Adjusted = end,
          Arm = "p", log2_CR = log2_CR, Segment_Length = end - start
        )
      } else if (start > centromere_pos) {
        # Segment fully in q-arm
        split_segments[[length(split_segments) + 1]] <- data.table(
          Chromosome = row$Chromosome, Start_Adjusted = start, End_Adjusted = end,
          Arm = "q", log2_CR = log2_CR, Segment_Length = end - start
        )
      } else {
        # Segment crosses centromere, split into p and q parts
        split_segments[[length(split_segments) + 1]] <- data.table(
          Chromosome = row$Chromosome, Start_Adjusted = start, End_Adjusted = centromere_pos,
          Arm = "p", log2_CR = log2_CR, Segment_Length = centromere_pos - start
        )
        split_segments[[length(split_segments) + 1]] <- data.table(
          Chromosome = row$Chromosome, Start_Adjusted = centromere_pos + 1, End_Adjusted = end,
          Arm = "q", log2_CR = log2_CR, Segment_Length = end - (centromere_pos + 1)
        )
      }
    }

    # Combine the split segments
    arm_segments <- rbindlist(split_segments)

    # Filter only the arm of interest
    arm_segments <- arm_segments[Arm == arm]

    # Check if at least 50% of the arm is affected
    affected_length <- arm_segments[(log2_CR >= 0.25 & chr_arm %in% gain_arms) | (log2_CR <= -0.25 & chr_arm %in% loss_arms), sum(Segment_Length, na.rm = TRUE)]

    if (affected_length >= half_arm_length) {
      sample_result[chr_arm] <- "yes"
    }
  }

  # Return results as data.table
  return(data.table(Arm = names(sample_result), Sample = sample_name, Status = sample_result))
}

# Process all seg files
seg_files <- list.files("./ULP_WGS_ichor_sln_other", pattern = "\\.seg\\.txt$", full.names = TRUE, recursive = TRUE)
all_results <- rbindlist(lapply(seg_files, process_seg_file), fill = TRUE)

# Reshape to have arms as rows and samples as columns
final_results <- dcast(all_results, Arm ~ Sample, value.var = "Status", fill = "no")
write_clip(final_results)





