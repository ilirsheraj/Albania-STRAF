# Prepared by Ilir Sheraj (Dr. Gjilpera)
# This script prepares a text file for input in STRUCTURE software
library(poppr)
library(dplyr)
library(tibble)

# Save output in STRUCTURE folder
structure_dir <- file.path(getwd(), "STRUCTURE_analysis")
dir.create(structure_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------
# Part 1: Create Region MetaData
# ------------------------------------------------
# This is an excel file with anonymized IDs for each person, City of origin 
# And based on some historical info, they are converted into 3 regions,
# South, North and Center
region_meta <- readxl::read_xlsx("QARQET-te dhena 2000 profile.xlsx", sheet = 1)
# head(region_meta)
region_meta <- region_meta %>%
  dplyr::select(`Sample URN`, `Item No`)
colnames(region_meta) <- c("City", "ind")

# Based on city, determine the regions
region_meta <- region_meta %>%
  mutate(
    City = trimws(toupper(City)),
    pop = case_when(
      City %in% c("SHKODER", "KUKES", "DIBER", "LEZHE") ~ "NORTH",
      City %in% c("TIRANE", "DURRES", "ELBASAN") ~ "CENTER",
      City %in% c("FIER", "KORCE", "BERAT", "VLORE", "GJIROKASTER") ~ "SOUTH"
    )
  ) %>%
  as.data.frame()

# Check the distribution
table(region_meta$City, useNA = "ifany")
table(region_meta$pop, useNA = "ifany")

# --------------------------------------------------------
# Allele Frequency Table
# --------------------------------------------------------
dat <- read.table("STRAF 2000 profile Full name Loci.txt",
                  header = TRUE, sep = "\t", 
                  stringsAsFactors = FALSE, check.names = FALSE)

regions <- region_meta[match(dat$ind, region_meta$ind), ]
all(regions$ind == dat$ind)

ind <- dat$ind
pop <- regions$pop

# Convert the table into a Diploid Frequency Table
geno2 <- dat[, -(1:2)]

# One Column per Locus
loci <- colnames(geno2)[seq(1, ncol(geno2), by = 2)]

# Numerical regional codes for STRUCTURE If you wanna use later
region_code <- c(
  "NORTH"  = 1L,
  "CENTER" = 2L,
  "SOUTH"  = 3L
)

structure_pop <- unname(region_code[as.character(pop)])

stopifnot(
  length(structure_pop) == nrow(dat),
  !anyNA(structure_pop)
)

# Allele columns from the original STRAF table
structure_genotypes <- dat[, -(1:2), drop = FALSE]

# This was the shitty part: STRUCTURE requires integer allele labels, so all those 
# with decimals have to be converted into full integers, so ultiplication by 10 
# preserves microvariant alleles such as 9 -> 90; 9.3 -> 93; 31.2 -> 312
structure_genotypes[] <- lapply(
  structure_genotypes,
  function(x) {
    x <- suppressWarnings(as.numeric(x))
    ifelse(is.na(x), -9L, as.integer(round(x * 10)))
  }
)

# One row per individual: ID, sampling-region code, then two allele columns per locus
structure_data <- data.frame(
  Label = ind,
  PopData = structure_pop,
  structure_genotypes,
  check.names = FALSE
)

# Another critical step here ;)
# STRUCTURE requires one marker name per locus in the header, although there 
# are two allele columns per locus in each individual row.
marker_header <- c("", "", loci)

structure_file <- file.path(structure_dir, "Albania_2000_STRUCTURE.txt")

# Write marker-name row
write.table(
  t(marker_header),
  file = structure_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE)

# Append the individual genotype data
write.table(
  structure_data,
  file = structure_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE,
  na = "-9")

# Save the regional code in case you wanna see if there is heterogeneity among regions
write.csv(
  data.frame(
    Region = names(region_code),
    Structure_code = unname(region_code)),
  file.path(structure_dir, "STRUCTURE_region_codes.csv"),
  row.names = FALSE)

# EOF