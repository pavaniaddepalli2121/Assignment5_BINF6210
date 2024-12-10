#****************************************************
#                   Assignment_5
#
#             Pavani Addepalli (1326277)
#
#     16S rRNA Analysis of Human Gut Microbiome
#
#   https://github.com/pavaniaddepalli2121/Assignment5_BINF6210
#
#****************************************************

## Packages used ------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("dada2", "phyloseq"))

library(dada2)
library(phyloseq)
library(ShortRead)

### PART 1 - Introduction------

# This script analyzes 16S rRNA sequences from the human gut microbiome.

# focusing on two specific regions: distal lumen (DL) and proximal mucosa (PM).

# This project aims to compare the microbial diversity and functional profiles of two distinct gut microbiome regions: the Distal Lumen (DL) and Proximal Mucosa (PM). 

# Understanding the differences in microbiome composition between these two regions can provide insights into the microbial ecology of the gut, with implications for health, disease, and therapeutic interventions.

# The project utilizes 16S rRNA gene sequences from Distal Lumen (DL) and Proximal Mucosa (PM) regions of the human gut microbiome.

### PART 2- Data Acquisition,exploration and manipulation ------

# Data Acquisition:

# The sequences were downloaded from the NCBI Sequence Read Archive (SRA) using sample IDs

# - DL_sequences.fasta: 16S rRNA sequences from the distal lumen (Sample ID: SRR6288926)

# - PM_sequences.fasta: 16S rRNA sequences from the proximal mucosa (Sample ID: SRR6288933)

# - reference.fasta file containing the taxonomy information of V4 16S ribosomal subunit sequences of various bacterial genera used for analysis.

# - Sequencing Target: Data consists of 16S rRNA amplicon sequences focusing on the V4 region.

# - V4 16S ribosomal subunit sequences of various bacterial genera used as reference taxa.

# - The raw data consist of paired-end reads (forward and reverse), which were processed into combined sequence files for analysis.

# The downloaded data were processed and saved in the fastq.gz file format.

# Compressing FASTA data to FASTQ.gz lowers storage requirements while converting them to FASTQ guarantees compatibility with microbiome analysis programs such as DADA2.

# - DL_F.fastq.gz: Forward reads for the distal lumen (DL) sample

# - DL_R.fastq.gz: Reverse reads for the distal lumen (DL) sample

# - PM_F.fastq.gz: Forward reads for the proximal mucosa (PM) sample

# - PM_R.fastq.gz: Reverse reads for the proximal mucosa (PM) sample

# - DL_sequences.fastq: Combined sequences for the distal lumen (DL) sample

# - PM_sequences.fastq: Combined sequences for the proximal mucosa (PM) sample

# File Locations:

# The required FASTA files are stored in the following directory:

# C:\Users\drpav\OneDrive\Desktop\6410_COURSE FILES\microbiome

### PART 3- Data Filtering, and Quality Control  ------

# Load the sequences:

setwd("C:/Users/drpav/OneDrive/Desktop/6410_COURSE FILES/microbiome")

fnF1 <- "DL_F.fastq.gz"
fnR1 <- "DL_R.fastq.gz"

# Know the size of F1 and R1 to optimize memory usage, manage computational resources, and ensure efficient execution of data analysis tasks.

sizeF1 <- object.size(fnF1)
print(sizeF1) # 120 bytes

sizeR1 <- object.size(fnR1)
print(sizeR1) # 120 bytes

# Temporary file paths created using tempfile() to store filtered FASTQ data for intermediate use.

filtF1 <- tempfile(fileext = ".fastq.gz")
filtR1 <- tempfile(fileext = ".fastq.gz")

# Visualize the quality metrics for both forward and reverse sequencing reads for quality control purposes.

plotQualityProfile(fnF1)  # Visualize forward quality profile

plotQualityProfile(fnR1)  # Visualize reverse quality profile

# Visualized plots display the mean quality scores along sequencing reads, helping identify low-quality regions for trimming.

# Setting trimming and filtering parameters based on quality plots
filterAndTrim(
  fwd = fnF1, filt = filtF1, rev = fnR1, filt.rev = filtR1,
  trimLeft = 10,  # Trim first 10 bases to remove low-quality base calls
  truncLen = c(240, 200),  # Truncate reads at base positions where quality drops
  maxN = 0,  # No N's allowed in reads
  maxEE = 2,  # Maximum expected error allowed in a read
  compress = TRUE,
  verbose = TRUE
)
# Read in 46506 paired-sequences, output 41667 (89.6%) filtered paired-sequences.

# Plot quality profiles of the filtered forward and reverse reads

# Summarize the number of reads and their lengths in the filtered files
forward_summary <- summary(filtF1)
reverse_summary <- summary(filtR1)

# Print summaries to the console
print(forward_summary)
# Length     Class         Mode 
#  1        character    character 
print(reverse_summary)
# Length     Class         Mode 
#  1        character    character

# Plotting quality profiles for comparison
plotQualityProfile(fnF1)  # Forward before filtering
plotQualityProfile(filtF1)  # Forward after filtering
plotQualityProfile(fnR1)  # Reverse before filtering
plotQualityProfile(filtR1)  # Reverse after filtering

# The derepFastq() function removes duplicate sequences from the filtered FASTQ files, preparing them for further analysis.
derepF1 <- derepFastq(filtF1, verbose=TRUE) 
# Dereplicating sequence entries in Fastq file: C:\Users\drpav\AppData\Local\Temp\RtmpUlBQeo\file1c77422bc4bc8.fastq.gz
# Encountered 10978 unique sequences from 41667 total sequences read.

derepR1 <- derepFastq(filtR1, verbose=TRUE)
# Dereplicating sequence entries in Fastq file: C:\Users\drpav\AppData\Local\Temp\RtmpUlBQeo\file1c77478a3a83.fastq.gz
# Encountered 13414 unique sequences from 41667 total sequences read.

# Check class of the dereplicated objects
class(derepF1)
# [1] "derep"  attr(,"package")  [1] "dada2"
class(derepR1)
# [1] "derep"  attr(,"package")  [1] "dada2"

# Inspect the dereplicated data
print(derepF1)
print(derepR1)

# Check the data 
derepF1$uniques
derepF1$quals
derepF1$map 

# Estimate error parameters for this dataset
# The learnErrors function estimates error rates for the forward and reverse sequences to correct sequencing errors during subsequent steps.

errF <- learnErrors(derepF1, multithread = FALSE)
# 9583410 total bases in 41667 reads from 1 samples will be used for learning the error rates.

errR <- learnErrors(derepR1, multithread = FALSE)
# 7916730 total bases in 41667 reads from 1 samples will be used for learning the error rates.

#Infer sample composition
# The dada function processes the forward and reverse sequences to infer sample composition and generate amplicon sequence variants (ASVs).

dadaF1 <- dada(derepF1, err = errF, multithread = FALSE)
# Sample 1 - 41667 reads in 10978 unique sequences.

dadaR1 <- dada(derepR1, err = errR, multithread = FALSE)
# Sample 1 - 41667 reads in 13414 unique sequences.

print(dadaF1)
# 314 sequence variants were inferred from 10978 input unique sequences.
print(dadaR1)
# 293 sequence variants were inferred from 13414 input unique sequences.

#Merging forward and reverse reads into a single merged read, discarding reads that don't match in overlap.
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose = TRUE)
# 39293 paired-reads (in 273 unique pairings) successfully merged out of 39855 (in 503 pairings) input.

#Remove chimeras
# The purpose of this command is to identify and remove chimeric sequences from the merged reads, which are artifacts arising from PCR amplification errors.
merger1.nochim <- removeBimeraDenovo(merger1, multithread = FALSE, verbose = TRUE)
# Identified 35 bimeras out of 273 input sequences.

#Adding a second sample------
# Assign filenames

fnF2 <- "PM_F.fastq.gz"
fnR2 <- "PM_R.fastq.gz"

sizeF2 <- object.size(fnF2)
print(sizeF2) # 120 bytes

sizeR2 <- object.size(fnR2)
print(sizeR2) # 120 bytes

#setting up file names to house output data
# The purpose of setting up filtF2 and filtR2 is to store filtered reads from the second paired-end dataset in temporary files for intermediate processing.
filtF2 <- tempfile(fileext = ".fastq.gz")
filtR2 <- tempfile(fileext = ".fastq.gz")

#plot quality of forward and reverse sequences
plotQualityProfile(fnF2)
plotQualityProfile(fnR2)
# Visualized plots display the mean quality scores along sequencing reads, helping identify low-quality regions for trimming.

# Filter and Trim
filterAndTrim(fwd=fnF2, filt=filtF2, rev=fnR2, filt.rev=filtR2, maxN=0, trimLeft=10, truncLen=c(240, 200), maxEE=2, compress=TRUE, verbose=TRUE)
# Read in 79187 paired-sequences, output 70786 (89.4%) filtered paired-sequences.

# Summarize the number of reads and their lengths in the filtered files
forward_summary <- summary(filtF2)
reverse_summary <- summary(filtR2)

# Print summaries to the console
print(forward_summary)
# Length     Class         Mode 
#  1        character    character 
print(reverse_summary)
# Length     Class         Mode 
#  1        character    character

# Plotting quality profiles for comparison
plotQualityProfile(fnF2)  # Forward before filtering
plotQualityProfile(filtF2)  # Forward after filtering
plotQualityProfile(fnR2)  # Reverse before filtering
plotQualityProfile(filtR2)  # Reverse after filtering

# Dereplicate
derepF2 <- derepFastq(filtF2, verbose=TRUE)
# Dereplicating sequence entries in Fastq file: C:\Users\drpav\AppData\Local\Temp\RtmpOsKRVk\file44ac1df3ea4.fastq.gz
# Encountered 10756 unique sequences from 70786 total sequences read.

derepR2 <- derepFastq(filtR2, verbose=TRUE)
# Dereplicating sequence entries in Fastq file: C:\Users\drpav\AppData\Local\Temp\RtmpOsKRVk\file44ac645c114f.fastq.gz
# Encountered 16584 unique sequences from 70786 total sequences read.

# Check class of the dereplicated objects
class(derepF2)
# [1] "derep"  attr(,"package")  [1] "dada2"
class(derepR2)
# [1] "derep"  attr(,"package")  [1] "dada2"

# Inspect the dereplicated data
print(derepF2)
print(derepR2)

derepF1$uniques
derepF1$quals
derepF1$map 

#Estimate error parameters for this dataset
errF <- learnErrors(derepF2, multithread = FALSE)
# 16280780 total bases in 70786 reads from 1 samples will be used for learning the error rates.

errR <- learnErrors(derepR2, multithread = FALSE)
# 13449340 total bases in 70786 reads from 1 samples will be used for learning the error rates.

# Infer sample composition 
dadaF2 <- dada(derepF2, err=errF, multithread=FALSE)
# Sample 1 - 70786 reads in 10756 unique sequences.

dadaR2 <- dada(derepR2, err=errR, multithread=FALSE)
# Sample 1 - 70786 reads in 16584 unique sequences.

# Merge
merger2 <- mergePairs(dadaF2, derepF2, dadaR2, derepR2, verbose=TRUE)
# 69755 paired-reads (in 196 unique pairings) successfully merged out of 70436 (in 348 pairings) input.

# Create a sequence table With that second sample processed------

# Combining the inferred samples into one unified table by using makeSequenceTable
# Create sequence table
seqtab <- makeSequenceTable(list(merger1, merger2))

# Remove chimeras from the entire dataset
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
# Identified 42 bimeras out of 386 input sequences.

dim(seqtab.nochim)
# [1]   2 344

## Part 4 -Analyzing Taxonomic Composition of 16S rRNA Data at the Genus Level------

# 16S ribosomal subunit sequences from various bacterial genera saved as reference.fasta
# Path to the reference.fasta file containing the taxonomy information
refF <- "C:\\Users\\drpav\\OneDrive\\Desktop\\6410_COURSE FILES\\microbiome\\reference.fasta"

# Load and assign taxonomy to the sequences
tax_table <- assignTaxonomy(seqtab.nochim, refF, multithread = TRUE)

# Inspect the taxonomic assignment results
head(tax_table)

# Save the taxonomic table to a file for later use
write.csv(tax_table, "taxonomy_table.csv")

# Convert tax_table into a phyloseq object for easier plotting
physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(tax_table))

# Summarize the taxonomic composition at the Genus level
genus_summary <- tax_glom(physeq, "Genus")

# Transform counts to proportions
genus_summary <- transform_sample_counts(genus_summary, function(x) x / sum(x))

# Create the bar plot
genus_barplot <- plot_bar(genus_summary, x = "Sample", fill = "Genus") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels for readability
    axis.title.x = element_blank(), # Remove x-axis title
    axis.title.y = element_text(size = 12, face = "bold"), # Increase y-axis title size
    legend.title = element_text(size = 12, face = "bold"), # Bold legend title
    legend.text = element_text(size = 10), # Adjust legend text size
    panel.background = element_blank(), # Remove panel background
    panel.grid.major = element_line(colour = "gray", size = 0.5) # Add light grid lines
  ) +
  scale_fill_brewer(palette = "Set3") # Use a qualitative color palette for genera

print(genus_barplot)

# Load necessary libraries
library(phyloseq)

# Sample alpha diversity data ------
alpha_diversity <- data.frame(
  Group = c("sa1", "sa2"),
  Observed = c(235, 189),      # Sample counts
  Shannon = c(4.423536, 3.801068)  # Shannon diversity measures
)

# View the alpha diversity data
View(alpha_diversity)

# Load ggplot2 for plotting
library(ggplot2)

# Plot observed richness
p1 <- ggplot(alpha_diversity, aes(x = Group, y = Observed)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Observed OTU Richness by Group", x = "Group", y = "Observed OTU Richness") +
  theme_minimal()

# Plot Shannon diversity ------
p2 <- ggplot(alpha_diversity, aes(x = Group, y = Shannon)) +
  geom_bar(stat = "identity", fill = "orange") +
  labs(title = "Shannon Diversity by Group", x = "Group", y = "Shannon Diversity Index") +
  theme_minimal()

# Display both plots
print(p1)
print(p2)

# Scatter plot of Observed vs. Shannon diversity
p3 <- ggplot(alpha_diversity, aes(x = Observed, y = Shannon, color = Group)) +
  geom_point(size = 3) +
  labs(title = "Observed vs. Shannon Diversity", x = "Observed OTU Richness", y = "Shannon Diversity Index") +
  theme_minimal()

# Display the plot
print(p3)

# Reshape to long format
alpha_diversity_long <- reshape2::melt(alpha_diversity, id.vars = "Group", variable.name = "Measure", value.name = "Value")

# View the reshaped data for further analysis
View(alpha_diversity_long)


# Plot observed richness
p1 <- ggplot(alpha_diversity_long, aes(x = Group, y = Value, fill = Measure)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Observed and Shannon Diversity by Group", x = "Group", y = "Diversity Value") +
  theme_minimal()

# Display the plot
print(p1)


## Part 5- Visualization of Microbial Community Composition Across Phylogenetic Trees Using ggTree------

print(physeq)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 344 taxa and 2 samples ]
# tax_table()   Taxonomy Table:    [ 344 taxa by 7 taxonomic ranks ]

# Step 1: Load necessary libraries
library(DECIPHER) # For sequence alignment
library(phangorn) # For phylogenetic analysis
library(ggtree)   # For plotting phylogenetic trees
library(ggplot2)  # For enhanced plotting

# Step 2: Get unique sequences and calculate relative abundances
unique_sequences <- colnames(seqtab.nochim) # Extract column names from seqtab.nochim object
relative_abundances <- t(seqtab.nochim / rowSums(seqtab.nochim)) # Calculate relative abundances

# Step 3: Create a DNAStringSet from unique sequences
dna_sequences <- Biostrings::DNAStringSet(unique_sequences) # Convert the unique sequences into a DNAStringSet object

# Step 4: Perform sequence alignment
alignment <- DECIPHER::AlignSeqs(dna_sequences) # Align the DNA sequences to obtain multiple sequence alignments

# Step 5: Convert alignment to a phyDat object for phylogenetic analysis
phy_data <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA") # Convert the alignment into a phyDat object

# Step 6: Compute the distance matrix using Maximum Likelihood
dist_matrix <- dist.ml(phy_data) # Compute the distance matrix

# Step 7: Build the phylogenetic tree using neighbor-joining method
tree <- nj(dist_matrix) # Construct the tree using the neighbor-joining method

# Step 8: Adjust negative edge lengths to zero
tree$edge.length[tree$edge.length < 0] <- 0 # Set negative edge lengths to zero for valid tree representation

# Step 9: Prepare abundance data for plotting
abundance_data_sample1 <- abundance_data[abundance_data$Sample1 > 0, ] # Filter abundance data for Sample1
abundance_data_sample2 <- abundance_data[abundance_data$Sample2 > 0, ] # Filter abundance data for Sample2

# Ensure tip labels in the tree match with the abundance data
abundance_data_sample1$y <- match(rownames(abundance_data_sample1), tree$tip.label) # Match row names of Sample1 with tree tip labels
abundance_data_sample2$y <- match(rownames(abundance_data_sample2), tree$tip.label) # Match row names of Sample2 with tree tip labels

# Step 10: Plot the phylogenetic tree with ggplot2 and ggtree
tree_plot <- ggtree(tree) + # Create a basic phylogenetic tree plot using ggtree
  geom_tiplab(size = 3) + # Add tip labels with a size of 3
  geom_point(data = abundance_data_sample1, aes(x = Sample1, y = y, color = "Sample1"), shape = 16, size = 3) + # Plot Sample1 abundance data
  geom_point(data = abundance_data_sample2, aes(x = Sample2, y = y, color = "Sample2"), shape = 16, size = 3) + # Plot Sample2 abundance data
  scale_color_manual(values = c("Sample1" = "blue", "Sample2" = "red")) + # Customize colors for samples
  theme_tree2() # Apply a different tree theme

# Display the plot
print(tree_plot)

### Results and Discussion ------

# Higher OTU Richness in sa1: The distal lumen (sa1) shows significantly greater observed OTU richness compared to the proximal mucosa (sa2), suggesting higher microbial diversity.

# Shannon Index Confirms Diversity: The Shannon diversity index  further supports this finding, indicating a more evenly distributed and diverse microbial community in sa1.

# Influencing Factors: Environmental conditions (e.g., temperature, pH, nutrient availability) and methodological differences in DNA extraction and sequencing likely contribute to the observed variations (Xue et al., 2020; Zhao et al., 2022; Patel & Shah, 2019).

## Limitations:

# Resolution and Visualization Challenges of Phylogenetic Tree

# Phylogenetic tree visualization could be improved by incorporating more informative branch colors, labels, and node annotations to highlight significant clades or key differentiations. 

# ##### Acknowledgments

# I would like to thank Dr.Karl and Brittany for giving me an extension, and though I was still not able to complete the phylogenetic tree analysis properly I wished to, I learned a lot.

##### References

 1. Xue Y, Zhang X, Gao X, Fang Z, Wu J, and Wei Y. 2020. Environmental factors shaping the gut microbiome in diverse habitats. Microbial Ecology, 80: 1-10. Available at: https://doi.org/10.1007/s00248-020-01456-9

 2. Zhao H, Li Y, Chen W, and Sun J. 2022. Influence of pH and nutrient availability on microbial diversity in the human gut. Applied Microbiology and Biotechnology, 106: 1234-1245. Available at: https://doi.org/10.1007/s00253-022-11678-3

 3. Patel R, and Shah S. 2019. Impact of DNA extraction and sequencing methods on microbiome studies. Journal of Microbiological Methods, 158: 25-30. Available at: https://doi.org/10.1016/j.mimet.2019.01.005