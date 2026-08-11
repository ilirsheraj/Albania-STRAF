# Prepared by Ilir Sheraj (Dr. Gjilpera)
# This script downloades frequency data from STRidER Database,
# prepares frequency table for ALbania and a number of other countries for 
# Population level analysis. Then it generates distances 
library(tidyverse)
library(ape)
library(forensicpopdata)
library(straf)
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
library(ggrepel)
library(graph4lg)
# -------------------------------------------
# Part 1: Extract data from STRidER
# -------------------------------------------
freqs <- read_STRidER_xml()

# See what populations exist
names(freqs)

# Remove database/continental aggregate entries and other unrelated countries
exclude <- c(
  "Entire Database",
  "Africa",
  "America",
  "Asia",
  "Europe",
  "SAUDI ARABIA",
  "KENYA",
  "VIETNAM",
  "THAILAND",
  "DOMINICAN REPUBLIC",
  "SOUTH AFRICA"
)

freqs_country <- freqs[!names(freqs) %in% exclude]

# Countries retained
names(freqs_country)
length(freqs_country)

# Keep this so you dont have to download it over and over since the data
# to be used for comparison will keep changing all the time
freqs_sub <- freqs_country

# Some Diagnostics
# Check the Loci
# All locus names occurring anywhere
all_loci <- sort(unique(unlist(lapply(freqs_sub, names))))
all_loci

# Union of alleles for each locus
allele_union <- setNames(
  lapply(all_loci, function(locus) {
    
    alleles <- unlist(
      lapply(freqs_sub, function(pop) {
        if (locus %in% names(pop)) {
          names(pop[[locus]])
        } else {
          NULL
        }
      })
    )
    
    sort(unique(alleles))
  }),
  all_loci
)

allele_union

# Get the list of common loci for all
loci_per_pop <- lapply(freqs_sub, names)
common_loci <- Reduce(intersect, loci_per_pop)
common_loci
length(common_loci)

# ----------------------------------------------
# Part 2: Create the Albanian Object
# ----------------------------------------------
alb_path <- "STRAF 2000 profile Full name Loci.txt"

alb <- read.delim(alb_path, check.names = FALSE, stringsAsFactors = FALSE)

# keep only ALBANIA rows
alb <- alb[alb$pop == "ALBANIA", ]

# locus names are duplicated, so get the unqiue ones
locus_cols <- names(alb)[-(1:2)]
loci_alb <- unique(locus_cols)

locus_freq <- function(df, locus_name){
  idx <- which(names(df) == locus_name)
  a1 <- df[[idx[1]]]; a2 <- df[[idx[2]]]
  alleles <- c(a1, a2) %>% as.character()
  alleles <- alleles[!is.na(alleles) & alleles != "R" & alleles != ""]
  tab <- table(alleles)
  tibble(Allele = names(tab), Freq = as.numeric(tab) / sum(tab))
}

albania <- setNames(lapply(loci_alb, function(L) locus_freq(alb, L)), loci_alb)

# Convert it to list
albania <- lapply(
  albania,
  function(x) {
    freq <- x$Freq
    names(freq) <- as.character(x$Allele)
    freq})

str(albania)

# Check they sum to 1
sapply(albania, sum)

# # Hairy Part: Make sure names are the same, otherwise it sucks
# intersect(names(albania), common_loci)
# # 15
# setdiff(names(albania), common_loci)
# # "SE33"
# setdiff(common_loci, names(albania))
# # character(0)

# -----------------------------------------------
# Part 3: Add other countries of Interest
# -----------------------------------------------
# This is quite repetitive and boring, most of the se files have the same
# format since the tables were extracted from papers pdfs, use this function to 
# return a clean dataframe which can ge further manipulated
df_to_af <- function(df) {
  
  # First column must be named Allele
  if (!"Allele" %in% names(df)) {
    stop("Dataframe must contain a column named 'Allele'.")
  }
  
  af_list <- lapply(
    names(df)[names(df) != "Allele"],
    function(locus) {
      
      x <- df[, c("Allele", locus)]
      
      # Keep only reported frequencies
      x <- x[!is.na(x[[locus]]), ]
      
      freq <- as.numeric(x[[locus]])
      names(freq) <- as.character(x$Allele)
      
      freq
    }
  )
  
  names(af_list) <- names(df)[names(df) != "Allele"]
  
  af_list
}

# ---------------------------------------------
# Kosovo
# ---------------------------------------------
kosovo_df <- read.csv(
  "NewPopulations/Kubat_2004_Kosovo_Albanians_Table1_allele_frequencies.csv",
  check.names = FALSE
)

kosovo <- df_to_af(kosovo_df)
sapply(kosovo, sum)

# intersect(names(kosovo), common_loci)
# # 10: "D8S1179" "D21S11"  "D3S1358" "TH01" "D16S539" "D2S1338" "D19S433" "VWA" 
# # "D18S51"  "FGA"
# 
# setdiff(names(kosovo), common_loci)
# # "D7S820"  "CSF1PO"  "D13S317" "TPOX"    "D5S818"
# 
# setdiff(common_loci, names(kosovo))
# # "D1S1656"  "D2S441"   "D10S1248" "D12S391"  "D22S1045"

# ---------------------------------------------------------------
# Serbia
# ---------------------------------------------------------------
# This is fucking problematic: The table is corrupted
serbia_df <- readxl::read_xlsx(
  "NewPopulations/serbia_new_data.xlsx", sheet = 1)

# Fix the discrepancy with row 2
row6  <- which(serbia_df$Allele == 6)
rowNA <- which(is.na(serbia_df$Allele))

freq_cols <- setdiff(names(serbia_df), "Allele")

for (loc in freq_cols) {
  
  x <- serbia_df[[loc]][row6]
  y <- serbia_df[[loc]][rowNA]
  
  x <- ifelse(is.na(x), 0, x)
  y <- ifelse(is.na(y), 0, y)
  
  new_value <- x + y
  
  # Keep genuinely empty cells as NA
  if (is.na(serbia_df[[loc]][row6]) &&
      is.na(serbia_df[[loc]][rowNA])) {
    new_value <- NA_real_
  }
  
  serbia_df[[loc]][row6] <- new_value
}

serbia_df <- serbia_df[-rowNA, ]

serbia <- df_to_af(serbia_df)

sapply(serbia, sum)

names(serbia)[names(serbia) == "vWA"] <- "VWA"

# intersect(names(serbia), common_loci)
# # 15: "D10S1248" "D12S391"  "D16S539"  "D18S51"   "D19S433"  "D1S1656"  "D21S11"   
# # "D22S1045" "D2S1338"  "D2S441"  "D3S1358"  "D8S1179"  "FGA" "TH01" "VWA"
# 
# setdiff(names(serbia), common_loci)
# # "CSF1PO"  "D13S317" "D5S818"  "D7S820"  "Penta D" "Penta E" "TPOX"
# 
# setdiff(common_loci, names(serbia))
# # character(0)

# -----------------------------------------------------
# Turkish Cyprus
# ----------------------------------------------------
cyprus_df <- readxl::read_xlsx(
  "NewPopulations/Turkish_Cypriot_15_STR_allele_frequencies.xlsx", sheet = 1)

cyprus <- df_to_af(cyprus_df)

names(cyprus)[names(cyprus) == "THO1"] <- "TH01"
names(cyprus)[names(cyprus) == "vWA"] <- "VWA"

sapply(cyprus, sum)

# intersect(names(cyprus), common_loci)
# # "D8S1179" "D21S11"  "D3S1358" "TH01"  "D16S539" "D2S1338" "D19S433" "VWA"
# # "D18S51"  "FGA"
# 
# setdiff(names(cyprus), common_loci)
# # "D7S820"  "CSF1PO"  "D13S317" "TPOX" "D5S818"
# 
# setdiff(common_loci, names(cyprus))
# # "D1S1656"  "D2S441"   "D10S1248" "D12S391"  "D22S1045"

# --------------------------------------------
# Turkey and Its Regions
# --------------------------------------------
# Turkish Data is in a single excel file: Use this function to extract numbers
# related to allele frequencies and ignore the rest of the file

read_af_sheet <- function(path, sheet = 1, skip = 0, n_max = Inf) {
  
  # Read sheet
  df <- readxl::read_excel(path, sheet = sheet, skip = skip, n_max = n_max)
  
  # Always change first column to Allele
  colnames(df)[1] <- "Allele"
  
  # Remove rows without an allele
  df <- df[!is.na(df$Allele), ]
  
  # Replace "-" with NA, otherwise it sucks
  df[df == "-"] <- NA
  
  # Convert all columns to numeric
  df[] <- lapply(df,
    function(x) suppressWarnings(as.numeric(as.character(x))))
  
  # Convert to locus -> named allele-frequency vectors
  af_list <- lapply(names(df)[-1], function(locus) {
      
      x <- df[, c("Allele", locus)]
      
      # Keep only reported frequencies
      x <- x[!is.na(x[[locus]]), ]
      
      freq <- x[[locus]]
      names(freq) <- as.character(x$Allele)
      
      freq
      })
  
  # Stick allele names
  names(af_list) <- names(df)[-1]
  
  return(af_list)
}

## Region 1: Marmaris without Istanbul
marmaris <- read_af_sheet(
  path = "NewPopulations/iahb_a_1183709_sm2605.xls",
  sheet = 1,
  skip = 2,
  n_max = 74
)

sapply(marmaris, sum)
# 
# intersect(names(marmaris), common_loci)
# # "D10S1248" "D12S391"  "D16S539"  "D18S51"   "D19S433"  "D1S1656"  "D21S11"   "D22S1045"
# # "D2S1338"  "D2S441"   "D3S1358"  "D8S1179"  "FGA"      "TH01"     "VWA"
# 
# setdiff(names(marmaris), common_loci)
# # character(0)
# 
# setdiff(common_loci, names(marmaris))
# # character(0)


# Region 2: Istanbul
istanbul <- read_af_sheet(
  path = "NewPopulations/iahb_a_1183709_sm2605.xls",
  sheet = 8,
  skip = 2,
  n_max = 71
)

sapply(istanbul, sum)

# intersect(names(istanbul), common_loci)
# # "D10S1248" "D12S391"  "D16S539"  "D18S51"   "D19S433"  "D1S1656"  "D21S11"   
# # "D22S1045" "D2S1338"  "D2S441"   "D3S1358"  "D8S1179"  "FGA"   "TH01"   "VWA"
# 
# setdiff(names(istanbul), common_loci)
# # character(0)
# 
# setdiff(common_loci, names(istanbul))
# # character(0)

# Region 3: Aegian
aegian <- read_af_sheet(
  path = "NewPopulations/iahb_a_1183709_sm2605.xls",
  sheet = 7,
  skip = 2,
  n_max = 48
)

sapply(aegian, sum)

# intersect(names(aegian), common_loci)
# # "D10S1248" "D12S391"  "D16S539"  "D18S51"   "D19S433"  "D1S1656"  "D21S11"   
# # "D22S1045" "D2S1338"  "D2S441"   "D3S1358"  "D8S1179"  "FGA"   "TH01"   "VWA"
# 
# setdiff(names(aegian), common_loci)
# # character(0)
# 
# setdiff(common_loci, names(aegian))
# # character(0)

# Region 4: Mediterrenean
med_region <- read_af_sheet(
  path = "NewPopulations/iahb_a_1183709_sm2605.xls",
  sheet = 4,
  skip = 2,
  n_max = 60
)

sapply(med_region, sum)

# intersect(names(med_region), common_loci)
# # "D10S1248" "D12S391"  "D16S539"  "D18S51"   "D19S433"  "D1S1656"  "D21S11"   
# # "D22S1045" "D2S1338"  "D2S441"   "D3S1358"  "D8S1179"  "FGA"   "TH01"   "VWA"
# 
# setdiff(names(med_region), common_loci)
# # character(0)
# 
# setdiff(common_loci, names(med_region))

# ---------------------------------------------------
# Part 4: Add all countries to make comparion
# ---------------------------------------------------
freqs_sub[["ALBANIA"]] <- albania
# freqs_sub[["KOSOVO"]] <- kosovo
# freqs_sub[["SERBIA"]] <- serbia
# freqs_sub[["TURKISH CYPRUS"]] <- cyprus
freqs_sub[["TURKEY MARMARIS"]] <- marmaris
freqs_sub[["TURKEY ISTANBUL"]] <- istanbul
freqs_sub[["TURKEY AEGIAN"]] <- aegian
freqs_sub[["TURKEY MED"]] <- med_region

length(freqs_sub)

# -------------------------------------------------
# Part 5: Find Common Loci
# ------------------------------------------------

loci_per_pop <- lapply(freqs_sub, names)
common_loci <- Reduce(intersect, loci_per_pop)
length(common_loci)

# Run some diagnostics here
freqs_shared <- lapply(freqs_sub, function(pop) pop[common_loci])


check_freq_sums <- function(pop) {
  sapply(pop[common_loci], function(x) {
    z <- af_extract(x)
    sum(z$freqs, na.rm = TRUE)
  })
}

freq_sums <- lapply(freqs_shared, check_freq_sums)

range(unlist(freq_sums))

freq_sum_df <- do.call(
  rbind,
  lapply(names(freq_sums), function(pop) {
    data.frame(
      Population = pop,
      Locus = names(freq_sums[[pop]]),
      Sum = as.numeric(freq_sums[[pop]]),
      row.names = NULL
    )
  })
)

write.csv(freq_sum_df, "straf_database_loci_frequency_diagnostic.csv")

freq_sum_df[
  freq_sum_df$Sum < 0.95 | freq_sum_df$Sum > 1.05,
]

# ----------------------------------------------
# Part 6: Distance Trees
# ----------------------------------------------
# Use the following function to extract loci distances
# Run boostrap and plot the data

# Helpers: extract allele freqs + Nei distances + D matrix
af_extract <- function(x){
  if (is.numeric(x) && !is.null(names(x))) {
    alleles <- names(x); freqs <- as.numeric(x)
    ok <- is.finite(freqs) & !is.na(alleles) & alleles != ""
    return(list(alleles = alleles[ok], freqs = freqs[ok]))
  }
  if (is.data.frame(x) || is.matrix(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    allele_col <- grep("allele", names(df), ignore.case = TRUE, value = TRUE)
    alleles <- if (length(allele_col) >= 1) as.character(df[[allele_col[1]]]) else rownames(df)
    
    freq_col <- grep("freq|frequency", names(df), ignore.case = TRUE, value = TRUE)
    if (length(freq_col) >= 1) {
      freqs <- suppressWarnings(as.numeric(df[[freq_col[1]]]))
    } else {
      scores <- sapply(names(df), function(nm) sum(is.finite(suppressWarnings(as.numeric(df[[nm]])))))
      freqs <- suppressWarnings(as.numeric(df[[names(df)[which.max(scores)]]]))
    }
    
    ok <- is.finite(freqs) & !is.na(alleles) & alleles != ""
    return(list(alleles = as.character(alleles[ok]), freqs = freqs[ok]))
  }
  stop("Unsupported locus format: ", paste(class(x), collapse = "/"))
}

af_standardize <- function(x) {
  
  z <- af_extract(x)
  
  alleles <- trimws(as.character(z$alleles))
  freqs <- as.numeric(z$freqs)
  
  ok <- (
    !is.na(alleles) &
      alleles != "" &
      is.finite(freqs) &
      freqs >= 0)
  
  alleles <- alleles[ok]
  freqs <- freqs[ok]
  
  # Combine duplicate allele rows, if present
  freq_by_allele <- tapply(freqs, alleles, sum)
  
  if (length(freq_by_allele) == 0L ||
    sum(freq_by_allele) <= 0) {
    stop("No valid positive allele frequencies found.")
  }
  
  # Ensures frequencies sum to one
  freq_by_allele / sum(freq_by_allele)
}

# For Nei’s standard distance, accumulate the identity components across loci 
# before taking the logarithm:

nei_standard_distance <- function(pop_A, pop_B, loci_use) {
  
  Jxy <- 0
  Jx  <- 0
  Jy  <- 0
  loci_used <- 0L
  
  for (loc in loci_use) {
    
    p <- af_standardize(pop_A[[loc]])
    q <- af_standardize(pop_B[[loc]])
    
    alleles <- union(names(p), names(q))
    
    p_aligned <- setNames(rep(0, length(alleles)), alleles)
    q_aligned <- p_aligned
    
    p_aligned[names(p)] <- p
    q_aligned[names(q)] <- q
    
    Jxy <- Jxy + sum(p_aligned * q_aligned)
    Jx  <- Jx  + sum(p_aligned^2)
    Jy  <- Jy  + sum(q_aligned^2)
    
    loci_used <- loci_used + 1L
  }
  
  if (loci_used == 0L || Jx <= 0 || Jy <= 0) {
    return(NA_real_)
  }
  
  I <- Jxy / sqrt(Jx * Jy)
  
  # Numerical protection
  I <- min(max(I, .Machine$double.eps), 1)
  
  -log(I)
}

# Build the distance matrix
build_nei_Dmat <- function(freqs_shared, loci_use) {
  
  pops <- names(freqs_shared)
  n <- length(pops)
  
  D <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops))
  
  if (n < 2L) {
    return(D)
  }
  
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      
      d <- nei_standard_distance(
        pop_A = freqs_shared[[pops[i]]],
        pop_B = freqs_shared[[pops[j]]],
        loci_use = loci_use
      )
      
      D[i, j] <- d
      D[j, i] <- d
    }
  }
  
  D
}

# -------------------------------------------------
# This is all testing before bootstrap
# -------------------------------------------------
# Calculate Distance Matrix
Dmat_nei <- build_nei_Dmat(freqs_shared, common_loci)

round(Dmat_nei, 5)

# Nei Distance
D_nei <- as.dist(Dmat_nei)

mds_nei <- cmdscale(D_nei,
                    k = 2,
                    eig = TRUE,
                    add = TRUE)

# These plots have no statistical power
tree_nei <- ape::nj(D_nei)

plot(tree_nei,
     main = "NJ tree (Nei distance)",
     use.edge.length = FALSE,
     cex = 0.7,
     font = 2,
     no.margin = TRUE)

plot(
  tree_nei,
  main = "Neighbor-joining tree based on Nei's standard genetic distance",
  use.edge.length = TRUE,
  cex = 0.7,
  font = 2,
  no.margin = TRUE
)

# Order the countries by closest distance to Albania
sort(Dmat_nei["ALBANIA", ])

# --------------------------------------------------
# Calculations with Bootstrap
# -------------------------------------------------
bootstrap_nj_loci_nei_parallel <- function(
    freqs_shared,
    B = 500,
    seed = 42,
    outgroup = NULL,
    return_boot_trees = FALSE,
    cores = max(1, parallel::detectCores() - 1)
) {
  
  loci_per_pop <- lapply(freqs_shared, names)
  common_loci <- Reduce(intersect, loci_per_pop)
  
  if (length(common_loci) < 2L) {
    stop("Too few loci shared across populations.")
  }
  
  # Reference tree
  D0 <- build_nei_Dmat(
    freqs_shared = freqs_shared,
    loci_use = common_loci)
  
  if (
    anyNA(D0) ||
    any(!is.finite(D0)) ||
    any(D0 < 0)
  ) {
    stop("Invalid values detected in the reference distance matrix.")
  }
  
  ref <- ape::nj(as.dist(D0))
  
  root_tip <- if (!is.null(outgroup)) {
    outgroup
  } else {
    ref$tip.label[1]
  }
  
  if (!root_tip %in% ref$tip.label) {
    stop("Specified outgroup is not present among populations.")
  }
  
  ref <- ape::root(
    ref,
    outgroup = root_tip,
    resolve.root = TRUE)
  
  ref <- ape::ladderize(ref)
  
  L <- length(common_loci)
  
  # Reproducible seeds for each replicate
  set.seed(seed)
  boot_seeds <- sample.int(
    .Machine$integer.max,
    B)
  
  boot_trees <- parallel::mclapply(
    seq_len(B),
    
    function(b) {
      
      set.seed(boot_seeds[b])
      
      loci_b <- sample(
        common_loci,
        size = L,
        replace = TRUE)
      
      Db <- build_nei_Dmat(
        freqs_shared = freqs_shared,
        loci_use = loci_b)
      
      if (
        anyNA(Db) ||
        any(!is.finite(Db)) ||
        any(Db < 0)
      ) {
        stop(
          "Invalid distance matrix in bootstrap replicate ",
          b
        )
      }
      
      tb <- ape::nj(as.dist(Db))
      
      tb <- ape::root(
        tb,
        outgroup = root_tip,
        resolve.root = TRUE)
      
      ape::ladderize(tb)
    },
    
    mc.cores = cores
  )
  
  counts <- ape::prop.clades(
    ref,
    boot_trees
  )
  
  support <- 100 * counts / B
  
  list(
    tree = ref,
    support = support,
    common_loci = common_loci,
    Dmat = D0,
    boot_trees = if (return_boot_trees) boot_trees else NULL
  )
}

# -------------------------------------------
# Get the Results
# -------------------------------------------
res_nei <- bootstrap_nj_loci_nei_parallel(
  freqs_shared,
  B = 500,
  seed = 42,)

# Kosovo and TCyprus lower number of markers to 10
# save(res_nei, file = "serbia_kosovo_cyprus_bootstrap.RData")
# Serbian data is problematic, removed
# save(res_nei, file = "kosovo_cyprus_noserbia_bootstrap.RData")
# Remove Serbia, Kosovo and Cyprus
save(res_nei, file = "loci15_bootstrap.RData")

res_nei$common_loci
length(res_nei$common_loci)

res_nei$support

# Save the stuff
pdf(file = "NJ_tree_Nei_500_bootstrap_15loci.pdf",
  width = 10,
  height = 7)

plot(
  res_nei$tree,
  main = "NJ Tree Nei's Standard Genetic distance",
  use.edge.length = TRUE,
  cex = 0.7,
  font = 2,
  no.margin = TRUE
)
dev.off()

# JPEG
jpeg(
  filename = "NJ_tree_Nei_500_bootstrap_15loci.jpg",
  width = 10,
  height = 7,
  units = "in",
  res = 600,
  quality = 100
)

plot(
  res_nei$tree,
  main = "",
  use.edge.length = TRUE,
  cex = 0.7,
  font = 2,
  no.margin = TRUE
)
dev.off()

# No edgelength
jpeg(
  filename = "NJ_tree_Nei_500_bootstrap_smaller_15loci.jpg",
  width = 5,
  height = 3,
  units = "in",
  res = 600,
  quality = 100
)

plot(
  res_nei$tree,
  main = "",
  use.edge.length = FALSE,
  cex = 0.5,
  font = 1,
  no.margin = TRUE
)
dev.off()

support_labels <- ifelse(
  res_nei$support >= 50,
  round(res_nei$support),
  ""
)

# -----------------------------------------------
# Distance Matrix
# -----------------------------------------------
# Order populations by Nei distance from Albania
alb_order <- names(
  sort(res_nei$Dmat["ALBANIA", ]))

alb_order

# Reorder original distance matrix
D_keep <- res_nei$Dmat[alb_order, alb_order]

# Keep lower triangle only
D_keep[upper.tri(D_keep, diag = FALSE)] <- NA

dm_long <- as.data.frame(D_keep) %>%
  rownames_to_column("Pop_1") %>%
  pivot_longer(
    -Pop_1,
    names_to = "Pop_2",
    values_to = "Dist"
  ) %>%
  filter(!is.na(Dist)) %>%
  mutate(
    Pop_1 = factor(Pop_1, levels = rev(alb_order)),
    Pop_2 = factor(Pop_2, levels = alb_order)
  )

matrix_plot <- ggplot(dm_long, aes(x = Pop_2, y = Pop_1, fill = Dist)) +
  geom_tile(color = "white",inewidth = 0.3) +
  scale_fill_viridis_c(
    option = "magma",
    direction = 1,
    name = "Nei D") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 7,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 7),
    plot.title = element_text(
      hjust = 0.5
    ),
    legend.position = c(0.78, 0.78),
    legend.background = element_rect(
      fill = alpha("white", 0.75),
      color = "white",
      linewidth = 0.3)) +
  labs(
    title = "Nei distance",
    x = NULL,
    y = NULL) +
  coord_fixed() +
  scale_x_discrete(
    drop = FALSE,
    limits = alb_order
  ) +
  scale_y_discrete(
    drop = FALSE,
    limits = rev(alb_order))

ggsave("loci_15_countries_distance_matrix.jpeg", 
       plot = matrix_plot, 
       width = 6, 
       height = 6, 
       units = "in")

# EOF
