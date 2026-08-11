# Prepared by Ilir Sheraj (Dr. Gjilpera)
# This script provides an additional analysis of regional variation.
# Here, I use AMOVA to assess how much of the variance is explained by
# the available metadata, including city, region, and, of course, individual.
library(poppr)
library(dplyr)
library(tibble)
library(ggplot2)

# Save output in regional folder
out_dir   <- file.path(getwd(), "Regions_outputs")
plot_dir  <- file.path(getwd(), "Regions_plots")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

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

# Convert them to diploid frequencies
geno1 <- as.data.frame(sapply(seq_along(loci), function(i) {
  paste(geno2[[2*i - 1]], geno2[[2*i]], sep = "/")}),
  stringsAsFactors = FALSE)

colnames(geno1) <- loci
rownames(geno1) <- ind

# Building the genind object
gen <- df2genind(geno1, 
                 sep = "/",
                 pop = pop,
                 ind.names = ind,
                 ploidy = 2,
                 type = "codom")

gen

# Make sure city_info is aligned with the genind individuals
city_info <- region_meta$City[match(indNames(gen), region_meta$ind)]

# Sanity check
stopifnot(all(!is.na(city_info)))
stopifnot(length(city_info) == nInd(gen))

# Make a Hierarchy dataframe (2 level actually, but fuck it)
hier_df <- data.frame(
  Region = as.character(pop(gen)),
  City   = city_info,
  stringsAsFactors = FALSE)

head(hier_df)

# Assign the hierarchy to genind object
strata(gen) <- hier_df

# Verify it worked
head(strata(gen))

# Run hierarchical AMOVA: poppr package
amova_result <- poppr.amova(
  gen,
  hier = ~Region/City,
  method = "ade4")

print(amova_result)

names(amova_result)

# Save the ANOVA Table
amova_anova <- as.data.frame(
  amova_result$results,
  check.names = FALSE) %>%
  rownames_to_column("Source")

write.csv(amova_anova, file.path(out_dir, "AMOVA_anova_table.csv"),
          row.names = FALSE)

# Percentage of Variance Explained
amova_variance <- as.data.frame(
  amova_result$componentsofcovariance,
  check.names = FALSE) %>%
  rownames_to_column("Source")

write.csv(amova_variance, file.path(out_dir, "AMOVA_variance_components.csv"),
          row.names = FALSE)

amova_phi <- as.data.frame(
  amova_result$statphi,
  check.names = FALSE) %>%
  rownames_to_column("Statistic")

# Phi Statistics
amova_phi

write.csv(amova_phi, file.path(out_dir, "AMOVA_phi_statistics.csv"),
          row.names = FALSE)

##############################################
# A bit of Polishing
amova_components <- amova_variance %>%
  mutate(Source = trimws(gsub("^Variations\\s+", "", Source)))

names(amova_components)[names(amova_components) == "%"] <- "Percent"

amova_components

var_plot_df <- data.frame(
  Source = amova_components$Source,
  Percent = amova_components$Percent,
  stringsAsFactors = FALSE)

# Clean and order components: Im bored, so im doing it long way :)
var_plot_df <- var_plot_df %>%
  filter(!grepl("total", Source, ignore.case = TRUE)) %>%
  mutate(
    Source = case_when(
      grepl("between region", Source, ignore.case = TRUE) &
        !grepl("city", Source, ignore.case = TRUE) ~
        "Among regions",
      
      grepl("city", Source, ignore.case = TRUE) &
        grepl("region", Source, ignore.case = TRUE) ~
        "Among cities within regions",
      
      grepl("samples", Source, ignore.case = TRUE) &
        grepl("city", Source, ignore.case = TRUE) ~
        "Among individuals within cities",
      
      grepl("within samples", Source, ignore.case = TRUE) ~
        "Within individuals",
      
      TRUE ~ Source),
    
    Source = factor(
      Source,
      levels = c(
        "Within individuals",
        "Among individuals within cities",
        "Among cities within regions",
        "Among regions"
      )),
    
    # Factor decimals
    Label = ifelse(
      Percent >= 1,
      sprintf("%.2f%%", Percent),
      ifelse(
        Percent >= 0.01,
        sprintf("%.3f%%", Percent),
        sprintf("%.4f%%", Percent)))
  )

# Use a color map
amova_cols <- c(
  "Among regions" = "#D55E00",
  "Among cities within regions" = "#E69F00",
  "Among individuals within cities" = "#56B4E9",
  "Within individuals" = "#999999"
)

# Create a bar plot showing percentage of variance explained by all components
p_full <- ggplot(var_plot_df,aes(x = Percent, y = Source, fill = Source)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = Label),
    hjust = -0.15,
    size = 3.2) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.30))) +
  scale_fill_manual(values = amova_cols) +
  labs(
    subtitle = "",
    x = "Percentage of total variance",
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text = element_text(size = 10))

ggsave("AMOVA_All_Components.jpeg", 
       plot = p_full, 
       width = 5, 
       height = 2, 
       units = "in",
       path = plot_dir)

ggsave("AMOVA_All_Components.pdf", 
       plot = p_full, 
       width = 5, 
       height = 2, 
       units = "in",
       path = plot_dir)

# Another plot zooming in in three marginal components
small_df <- var_plot_df %>%
  filter(Source != "Within individuals") %>%
  mutate(
    Source = factor(
      Source,
      levels = c(
        "Among individuals within cities",
        "Among cities within regions",
        "Among regions")
      )
    )

p_small <- ggplot(small_df, aes(x = Percent, y = Source, fill = Source)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = Label),
    hjust = -0.15,
    size = 3.2) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.35))) +
  scale_fill_manual(values = amova_cols) +
  labs(
    subtitle = "",
    x = "Percentage of total variance",
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text = element_text(size = 10))


ggsave("AMOVA_selected_Components.jpeg", 
       plot = p_small, 
       width = 5, 
       height = 2, 
       units = "in",
       path = plot_dir)

ggsave("AMOVA_selected_Components.pdf", 
       plot = p_small, 
       width = 5, 
       height = 2, 
       units = "in",
       path = plot_dir)


# This is part of analysis to use permutation test
# Since AMOVA is very slow, unless you have 20 cores, dont try it
# We used STRUCTURE, so no need anyway but it was fun trying
# EOF
