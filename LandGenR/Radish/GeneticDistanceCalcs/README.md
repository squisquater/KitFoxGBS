There are several different ways to calculate genetic distance. 
* You can do it between individuals or groups of individuals (subpopulations).
* There are also several different metrics and the one you choose to calculate will depend on the downstream analyses you plan to run.


## Calculate Pairwise Fst and Nei's D
```
# Install and load required packages
install.packages("adegenet")
install.packages("vcfR")
install.packages("hierfstat")
install.packages("StAMPP")
library(adegenet)
library(vcfR)
library(hierfstat)
library(StAMPP)

# Set working directory to the directory containing your files
# setwd("/path/to/your/directory")

# List all files in the directory to verify the presence of populations_info.txt
print(list.files())

# Read the population information file
pop_data <- read.table("population_info.txt", header = TRUE)
print(pop_data)

# Read the VCF file
vcf <- read.vcfR("populations.snps.vcf")

# Convert VCF to genlight object
genlight_obj <- vcfR2genlight(vcf)

# Debugging step: check the sample IDs in the genlight object
print(indNames(genlight_obj))

# Ensure the sample IDs in pop_data match those in genlight_obj
if (!all(pop_data$IndividualID %in% indNames(genlight_obj))) {
  stop("Some sample IDs in populations_info.txt do not match those in the VCF file.")
} else {
  message("All sample IDs in populations_info.txt match those in the VCF file.")
}

# Assign population information to genlight object
pop(genlight_obj) <- pop_data$Population[match(indNames(genlight_obj), pop_data$IndividualID)]

# Verify the populations have been correctly assigned
print(pop(genlight_obj))

# Convert genlight object to StAMPP format
stampp_data <- stamppConvert(genlight_obj, "genlight")

# Calculate pairwise Fst
stampp_fst <- stamppFst(stampp_data, nboots = 100, percent = 95, nclusters = 1)

# Extract the Fst matrix
fst_matrix <- stampp_fst$Fsts

# Print the Fst matrix
print(fst_matrix)

# Save the Fst matrix to a CSV file
write.csv(fst_matrix, "pairwise_fst_matrix.csv", row.names = TRUE)

# Calculate Nei's genetic distance (D)
neis_d_matrix <- stamppNeisD(stampp_data)

# Print Nei's D matrix
print(neis_d_matrix)

# Save Nei's D matrix to a CSV file
write.csv(as.matrix(neis_d_matrix), "neis_d_matrix.csv", row.names = TRUE)
```
