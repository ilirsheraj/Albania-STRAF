library(tidyverse)
library(ape)
library(forensicpopdata)
library(straf)
library(adegenet)
library(hierfstat)
library(pegas)
library(poppr)
library(ggrepel)

# Working Here: "/home/ilir/EZaimi"
# Make two subdirectories
## One to save csv tables
out_dir   <- file.path(getwd(), "Regions_outputs")
## Second, to save the figures
plot_dir  <- file.path(getwd(), "Regions_plots")
# dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
# dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load the cities data
regions <- readxl::read_xlsx("QARQET-te dhena 2000 profile.xlsx", sheet = 1)
regions <- regions %>% dplyr::select(`Sample URN`, `Item No`)
colnames(regions) <- c("City", "ind")

regions <- regions %>% mutate(pop = case_when(
  City %in% c("SHKODER", "KUKES", "DIBER", "LEZHE") ~ "NORTH",
  City %in% c("TIRANE", "DURRES", "ELBASAN") ~ "CENTER",
  City %in% c("FIER", "KORCE", "BERAT", "VLORE","GJIROKASTER") ~ "SOUTH"
)) %>% as.data.frame()
head(regions)

# Load the frequency Table
dat <- read.table(
  "STRAF 2000 profile Full name Loci.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat[1:5, 1:6]
colnames(dat)

# Create the population index based on cities
regions <- regions[match(dat$ind, regions$ind), ]
all(regions$ind == dat$ind)

ind <- dat$ind
pop <- regions$pop

table(pop)
# CENTER  NORTH  SOUTH 
# 1039    358    603 

# Convert to Diploid Frequency Table
geno2 <- dat[, -(1:2)]

# One Column per Locus
loci <- colnames(geno2)[seq(1, ncol(geno2), by = 2)]

# Convert them to diploid frequencies
geno1 <- as.data.frame(
  sapply(seq_along(loci), function(i) {
    paste(geno2[[2*i - 1]], geno2[[2*i]], sep = "/")}),
  stringsAsFactors = FALSE
)

colnames(geno1) <- loci
rownames(geno1) <- ind


# Building the genind object
gen <- df2genind(
  geno1,
  sep = "/",
  pop = pop,
  ind.names = ind,
  ploidy = 2,
  type = "codom"
)

# Explore
# number of individuals
nInd(gen)

# Numbe rof Loci
nLoc(gen)

# Summary of Population Sizes: This case single population
table(pop(gen))

# Summary of All
summary(gen)

# Name sof Loci
locNames(gen)


# Visualize
af <- makefreq(gen, missing = "mean")

af_df <- as.data.frame(af) %>%
  rownames_to_column("ID") %>% 
  pivot_longer(-ID, names_to = "AlleleID", values_to = "Freq") %>%
  filter(!is.na(Freq) & Freq > 0) %>%
  separate(AlleleID, into = c("Locus", "Allele"), sep = "\\.", extra = "merge") %>%
  mutate(Allele = gsub("_", ".", Allele))

ggplot(af_df, aes(x = as.numeric(Allele), y = Freq)) +
  geom_col(width = 0.3, color = "brown1") +
  facet_wrap(~ Locus, scales = "free_x") +
  theme_classic() +
  labs(title = "Allele frequencies per locus", 
       x = "Alleles", 
       y = "Frequency") +
  theme(axis.text.x = element_text(size = 7),
        strip.text = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

# Summary Stats
stats <- hierfstat::basic.stats(gen)
# observed heterozygosity
stats$Ho
# expected heterozygosity
stats$Hs

# Convert this all into a clean table
rownames(stats$perloc) <- colnames(geno1)
perloc_df <- as.data.frame(stats$perloc)

perloc_df$Locus <- rownames(perloc_df)
perloc_df <- perloc_df[, c("Locus", setdiff(names(perloc_df), "Locus"))]

# Ho and Hs measure genetic diversity, Fis measures deviation from Hardy–Weinberg equilibrium, 
# and Fst quantifies how much populations differ genetically.
# perloc_df <- perloc_df[, !names(perloc_df) %in% c("Dstp", "Htp", "Fstp", "Dest")]
head(perloc_df)
write.csv(perloc_df, file = paste0(out_dir, "/Summary_statistics.csv"))

# Build a PCA
X <- scaleGen(gen, NA.method = "mean")
pca <- prcomp(X)
pve_percent <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

# Use ggplot to plot PCA
pop_cols <- c("#D55E00", "#CC79A7", "#332288")

# Use ggplot to plot PCA
pdf(file = paste0(plot_dir, "/Albanian_Population_PCA.pdf"), width = 5, height = 4)
as.data.frame(pca$x[, 1:2]) %>% 
  ggplot(aes(PC1, PC2, color = as.factor(pop))) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pop_cols, name = "Region") +
  labs(title = "Albanian Population PCA",
       x = paste0("PC1 (", pve_percent[1], "%)"),
       y = paste0("PC2 (", pve_percent[2], "%)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()

# Get a Scree Plot
pca_df <- tibble(
  PC = seq_along(pca$sdev),
  Variance = pca$sdev^2 / sum(pca$sdev^2)
)

pdf(file = paste0(plot_dir, "/Albanian_Population_Scree_Plot.pdf"), width = 6, height = 5)
ggplot(pca_df, aes(PC, Variance)) +
  geom_line() +
  geom_point(size = 1) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Variance Explained"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()

# t-SNE and UMAP: There is no structure here
# cumulative variance
cum_pve <- cumsum(pca$sdev^2) / sum(pca$sdev^2)

# see how many PCs give 50–70% variance
which(cum_pve >= 0.50)[1]
which(cum_pve >= 0.70)[1]

npc <- 86
X_pca <- pca$x[, 1:npc]

library(Rtsne)
set.seed(42)

perp <- min(30, floor((nrow(X_pca) - 1) / 3))

tsne <- Rtsne(
  X_pca,
  perplexity = perp,
  pca = FALSE,
  check_duplicates = FALSE
)

tsne_df <- data.frame(
  Dim1 = tsne$Y[, 1],
  Dim2 = tsne$Y[, 2],
  Region = as.factor(pop)
)

pdf(file = paste0(plot_dir, "/Albanian_Population_tSNE.pdf"), width = 5, height = 4)
ggplot(tsne_df, aes(Dim1, Dim2, color = Region)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pop_cols, name = "Region") +
  labs(title = paste0("Albanian Population t-SNE (perplexity=", perp, ")"),
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()


# UMAP
library(uwot)
set.seed(42)

umap <- umap(
  X_pca,
  n_neighbors = 15,
  min_dist = 0.1
)

umap_df <- data.frame(
  UMAP1 = umap[, 1],
  UMAP2 = umap[, 2],
  Region = as.factor(pop)
)

pdf(file = paste0(plot_dir, "/Albanian_Population_UMAP.pdf"), width = 5, height = 4)
ggplot(umap_df, aes(UMAP1, UMAP2, color = Region)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pop_cols, name = "Region") +
  labs(title = "Albanian Population UMAP",
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()

# Try a nuclear metric: Discriminant Analysis of PCA
set.seed(42)

# Pick a K (e.g., 3 if you have CENTER/NORTH/SOUTH; or try 2–6)
grp <- find.clusters(
  gen,
  max.n.clust = 10,
  n.pca = 50,
  n.clust = 3
)

str(grp)
table(grp$grp)

dapc1 <- dapc(gen, grp$grp, n.pca = 50, n.da = 2)
pdf(file = paste0(plot_dir, "/Albanian_Population_DA_PC.pdf"), width = 5, height = 4)
scatter(dapc1)
dev.off()

# Use regions and see if there is any clustering
dapc_reg <- dapc(gen, pop(gen), n.pca = 50, n.da = 2)
pdf(file = paste0(plot_dir, "/Albanian_Population_Region_DA_PC.pdf"), width = 5, height = 4)
scatter(dapc_reg)
dev.off()

dapc_reg <- dapc(gen, pop(gen), n.pca = 50, n.da = 2)
dapc_reg$eig


# Final, try different option
set.seed(42)
Ks <- 2:6
res <- lapply(Ks, function(K){
  grpK <- find.clusters(gen, max.n.clust = 10, n.pca = 50, n.clust = K)
  dapcK <- dapc(gen, grpK$grp, n.pca = 50, n.da = min(K-1, 2))
  list(K = K, grp = grpK, dapc = dapcK)
})

# Example: plot K = 3
scatter(res[[which(Ks==3)]]$dapc)

sapply(res, function(x) {
  sum(x$dapc$eig)
})

set.seed(42)
xval <- xvalDapc(
  gen,
  grp = res[[which(Ks == 3)]]$grp$grp,
  n.pca.max = 80,
  training.set = 0.9,
  result = "groupMean",
  center = TRUE,
  scale = FALSE
)

xval

# grp is the output of find.clusters()
cluster_table_long <- bind_rows(
  lapply(res, function(x) {
    data.frame(
      Ind     = names(x$grp$grp),
      K       = x$K,
      Cluster = as.integer(x$grp$grp),
      stringsAsFactors = FALSE)
  })) %>% 
  dplyr::filter(K == 3)

head(cluster_table_long)

# How many individuals per cluster per K
table(cluster_table_long$K, cluster_table_long$Cluster)

k3 <- res[[which(Ks == 3)]]
dapc_coords <- as.data.frame(k3$dapc$ind.coord)
dapc_coords$Ind <- rownames(dapc_coords)

dapc_membership <- dapc_coords %>%
  left_join(
    cluster_table_long %>% filter(K == 3),
    by = "Ind"
  )

head(dapc_membership)
write.csv(dapc_membership, file = paste0(out_dir, "/ld_on_pca_unsupervised_clustering.csv"))

set.seed(1)
rand_grp <- sample(pop(gen))   # random labels
dapc_rand <- dapc(gen, rand_grp, n.pca = 50, n.da = 2)
dapc_rand$eig
# Albanian STR variation shows no discrete regional genetic structure; 
# observed differentiation is weak, diffuse, and consistent with high gene flow and panmixia.

##############################
# MDS: Make it Stable, otherwise it sucks
replen <- rep(4, nLoc(gen))
names(replen) <- locNames(gen)
if ("D22S1045" %in% names(replen)) replen["D22S1045"] <- 3

# Check 
length(replen)
nLoc(gen)
all(names(replen) == locNames(gen))

# CMD By Regions
d_ind <- bruvo.dist(gen, replen = replen)
mds <- cmdscale(d_ind, k = 2, eig = TRUE)
var_exp <- round(100 * mds$eig / sum(mds$eig[mds$eig > 0]), 1)

mds_df <- data.frame(
  Pop = rownames(mds$points),
  Dim1 = mds$points[,1],
  Dim2 = mds$points[,2]
)

# Make sure they are on the same order, expected
all(mds_df$Pop == regions$ind)

mds_df$Region <- regions$pop

pdf(file = paste0(plot_dir, "/Albanian_Population_MDS.pdf"), width = 6, height = 4)
mds_df %>% ggplot(aes(Dim1, Dim2, color = as.factor(Region))) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = pop_cols, name = "Region") + 
  labs(title = "Albanian Population MDS",
       x = paste0("MDS1 (", var_exp[1], "%)"),
       y = paste0("MDS2 (", var_exp[2], "%)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()


# Pairwise Population Differentiation
# Weir & Cockerham’s FST (θ)
# How genetically differentiated are each pair of populations, relative to the total genetic variance?
hf <- genind2hierfstat(gen)
fst <- pairwise.WCfst(hf)
fst

write.csv(fst, file = paste0(out_dir, "/fst_population_differences.csv"))

# Negatives are useless, so lets remove them
fst2 <- fst
fst2[fst2 < 0] <- 0
range(fst, na.rm=TRUE)
range(fst2, na.rm=TRUE)
fst2
write.csv(fst2, file = paste0(out_dir, "/fst_population_differences_no_negative.csv"))

# Let's run a bootstrap
library(parallel)
boot_pwfst_mclapply <- function(hf, nboot = 100, ncores = 4, seed = 1){
  set.seed(seed)
  L <- ncol(hf) - 1
  pops <- levels(as.factor(hf[,1]))
  
  boots <- mclapply(seq_len(nboot), function(i){
    cols <- c(1, 1 + sample.int(L, L, replace = TRUE))
    pairwise.WCfst(hf[, cols, drop = FALSE])
  }, mc.cores = ncores)
  
  arr <- array(NA_real_, dim = c(length(pops), length(pops), nboot),
               dimnames = list(pops, pops, NULL))
  for(b in seq_len(nboot)) arr[,,b] <- boots[[b]]
  arr
}

# This takes about 15 min, 4 cores
arr <- boot_pwfst_mclapply(hf, nboot = 500, ncores = 4, seed = 42)

# Let's Visualize this
pops <- dimnames(arr)[[1]]

# build long data frame of all off-diagonal pairs
long <- expand.grid(Pop1 = pops, Pop2 = pops, stringsAsFactors = FALSE) %>%
  filter(Pop1 < Pop2) %>% 
  rowwise() %>%
  mutate(Fst = list(as.numeric(arr[Pop1, Pop2, ]))) %>%
  unnest(Fst) %>%
  ungroup() %>%
  mutate(Pair = paste(Pop1, "vs", Pop2))

# summary table (median + CI) for annotation
sumtab <- long %>%
  group_by(Pair) %>%
  summarise(
    med = median(Fst, na.rm = TRUE),
    lo  = quantile(Fst, 0.025, na.rm = TRUE),
    hi  = quantile(Fst, 0.975, na.rm = TRUE),
    .groups = "drop")

pdf(file = paste0(plot_dir, "/FST_500_Boostrap_Regions.pdf"), width = 6, height = 5)
ggplot(long, aes(Fst)) +
  geom_density(color = "red") +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_segment(data = sumtab,
               aes(x = lo, xend = hi, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 1.2) +
  geom_point(data = sumtab,
             aes(x = med, y = 0),
             inherit.aes = FALSE, size = 2) +
  facet_wrap(~ Pair, scales = "free_y") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)) +
  labs(title = "Bootstrap FST with 95% CI (segment) and median (dot)",
       x = "Weir–Cockerham FST (bootstrap)", 
       y = "Density")
dev.off()

# Hardy-WeinbergEquilibrium per locus
# HWE per locus (Monte Carlo / permutation based)
# B = number of permutations (bigger = more stable p-values)
hwe_res <- pegas::hw.test(gen, B = 10000)

hwe_res

hwe_df <- as.data.frame(hwe_res)
hwe_df <- hwe_df %>% rownames_to_column("Locus") %>% 
  mutate(FDR = p.adjust(Pr.exact, method = "fdr")) %>% 
  arrange(Pr.exact)
write.csv(hwe_df, file = paste0(out_dir, "/HWE_per_locus_10k_permutations.csv"), row.names = FALSE)

pdf(file = paste0(plot_dir, "/Per_Locus_HWE_10k_Permutations.pdf"), width = 4.5, height = 4)
hwe_df %>% ggplot(aes(x = reorder(Locus, Pr.exact), y = -log10(Pr.exact))) +
  geom_col() +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 12)) +
  labs(title = "HWE test per locus",
       subtitle = "Vertical Line p = 0.05",
       x = NULL, 
       y = expression(-log[10](p)))
dev.off()

################################################################################
# Linkage Dusequilibrium within Markers
# library(graph4lg)
allele_to3 <- function(a) {
  a <- as.character(a)
  a <- gsub("_", ".", a)
  a <- gsub("\\s+", "", a)
  if (is.na(a) || a == "" || a == "NA") return(NA_character_)
  a_num <- gsub("\\.", "", a)
  if (!grepl("^[0-9]+$", a_num)) return(NA_character_)
  sprintf("%03d", as.integer(a_num))
}

geno_to_genepop_code <- function(g) {
  if (is.na(g) || g == "" || g == "NA") return("000000")
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return("000000")
  a1 <- allele_to3(parts[1])
  a2 <- allele_to3(parts[2])
  if (is.na(a1) || is.na(a2)) return("000000")
  paste0(a1, a2)
}

run_ld_region <- function(region,
                          geno1, ind, pop,
                          out_dir,
                          plot_dir,
                          dememorization = 5000,
                          batches = 100,
                          iterations = 10000,
                          plot_width = 4,
                          plot_height = 4) {
  
  # subset
  if (is.null(region) || toupper(region) %in% c("ALL", "FULL", "TOTAL")) {
    idx <- seq_along(pop)
    region_name <- "ALL"
  } else {
    region_name <- as.character(region)
    idx <- which(pop == region_name)
    if (length(idx) == 0)
      stop("No samples found for region = ", region_name)
  }
  
  geno1_sub <- geno1[idx, , drop = FALSE]
  ind_sub   <- ind[idx]
  
  # encode genepop
  gp_geno <- as.data.frame(
    lapply(geno1_sub, function(col) vapply(col, geno_to_genepop_code, "")),
    stringsAsFactors = FALSE
  )
  
  loci <- colnames(gp_geno)
  
  # write genepop
  gp_tmp  <- tempfile(fileext = ".txt")
  out_tmp <- tempfile(fileext = ".txt")
  
  con <- file(gp_tmp, open = "wt")
  writeLines(paste0("Genepop export - ", region_name), con)
  writeLines(paste(loci, collapse = ", "), con)
  writeLines("Pop", con)
  
  for (i in seq_len(nrow(gp_geno))) {
    id <- gsub("[^A-Za-z0-9_\\-]", "_", ind_sub[i])
    gline <- paste(gp_geno[i, ], collapse = " ")
    writeLines(paste0(id, " , ", gline), con)
  }
  close(con)
  
  # run LD
  genepop::test_LD(
    inputFile = gp_tmp,
    outputFile = out_tmp,
    dememorization = dememorization,
    batches = batches,
    iterations = iterations
  )
  
  # read + clean LD results
  ld_df <- suppressWarnings(straf::get_ld_gp(out_tmp)) %>%
    mutate(
      Locus_1 = trimws(as.character(Locus_1)),
      Locus_2 = trimws(as.character(Locus_2)),
      P_value = suppressWarnings(as.numeric(gsub("^<", "", as.character(P_value))))
    ) %>%
    filter(!is.na(Locus_1), !is.na(Locus_2), !is.na(P_value)) %>%
    filter(Locus_1 %in% loci, Locus_2 %in% loci)
  
  # matrix
  ld_mat <- matrix(NA_real_, length(loci), length(loci),
                   dimnames = list(loci, loci))
  diag(ld_mat) <- 1
  
  if (nrow(ld_df) > 0) {
    for (i in seq_len(nrow(ld_df))) {
      a <- ld_df$Locus_1[i]
      b <- ld_df$Locus_2[i]
      p <- ld_df$P_value[i]
      ld_mat[a, b] <- p
      ld_mat[b, a] <- p
    }
  }
  
  # save matrix
  matrix_csv <- file.path(out_dir,
                          paste0("LD_matrix_", region_name, "_", iterations, "iter.csv"))
  write.csv(ld_mat, matrix_csv, row.names = TRUE)
  
  ld_mat_keep <- ld_mat
  ld_mat_keep[upper.tri(ld_mat_keep, diag = TRUE)] <- NA
  
  # Plot
  ld_long_keep <- as.data.frame(ld_mat_keep) %>%
    tibble::rownames_to_column("Locus_1") %>%
    tidyr::pivot_longer(-Locus_1, names_to = "Locus_2", values_to = "P_value") %>%
    dplyr::filter(!is.na(P_value)) %>%
    dplyr::mutate(
      # Create factors with all loci in both axes
      # Y-axis reversed for diagonal
      Locus_1 = factor(Locus_1, levels = rev(loci)),
      # X-axis in original order
      Locus_2 = factor(Locus_2, levels = loci)
    )
  
  p_plot <- ggplot(ld_long_keep, aes(x = Locus_2, y = Locus_1, fill = -log10(P_value))) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                         na.value = "transparent",
                         expression(-log[10](p))) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.78, 0.78),
      legend.background = element_rect(
        fill = alpha("white", 0.75),
        color = "white",
        linewidth = 0.3)
    ) +
    labs(
      title = paste0("Linkage Disequilibrium: ", region_name), 
      fill = "p", 
      x = NULL, 
      y = NULL
    ) +
    coord_fixed() +
    # Ensure all loci appear on both axes
    scale_x_discrete(drop = FALSE, limits = loci) +
    scale_y_discrete(drop = FALSE, limits = rev(loci))
  
  plot_pdf <- file.path(plot_dir,
                        paste0("LD_heatmap_", region_name, "_", iterations, "iter.pdf"))
  
  ggsave(
    filename = plot_pdf,
    plot = p_plot,
    device = cairo_pdf,
    width = plot_width,
    height = plot_height
  )
  
  list(
    region = region_name,
    n = length(idx),
    ld_df = ld_df,
    ld_mat = ld_mat,
    matrix_csv = matrix_csv,
    plot_pdf = plot_pdf,
    plot = p_plot
  )
}

# Wrapper Function
run_ld_regions <- function(geno1, ind, pop,
                           regions = "AUTO",
                           include_all = TRUE,
                           out_dir,
                           plot_dir,
                           dememorization = 5000,
                           batches = 100,
                           iterations = 10000,
                           plot_width = 4,
                           plot_height = 4) {
  
  if (length(regions) == 1 && regions == "AUTO") {
    regions <- sort(unique(as.character(pop)))
  }
  
  res <- list()
  
  for (r in regions) {
    res[[r]] <- run_ld_region(
      region = r,
      geno1 = geno1, ind = ind, pop = pop,
      out_dir = out_dir,
      plot_dir = plot_dir,
      dememorization = dememorization,
      batches = batches,
      iterations = iterations,
      plot_width = plot_width,
      plot_height = plot_height
    )
  }
  
  if (isTRUE(include_all)) {
    res[["ALL"]] <- run_ld_region(
      region = "ALL",
      geno1 = geno1, ind = ind, pop = pop,
      out_dir = out_dir,
      plot_dir = plot_dir,
      dememorization = dememorization,
      batches = batches,
      iterations = iterations,
      plot_width = plot_width,
      plot_height = plot_height
    )
  }
  
  res
}

north <- run_ld_region(
  region = "NORTH",
  geno1 = geno1, ind = ind, pop = pop,
  out_dir = out_dir,
  plot_dir = plot_dir
)

south <- run_ld_region(
  region = "SOUTH",
  geno1 = geno1, ind = ind, pop = pop,
  out_dir = out_dir,
  plot_dir = plot_dir
)

center <- run_ld_region(
  region = "CENTER",
  geno1 = geno1, ind = ind, pop = pop,
  out_dir = out_dir,
  plot_dir = plot_dir
)
################################################################################
# Create Phylogenetic Tree to compare regions
pop(gen) <- factor(pop(gen))

# Collapse individuals into populations
gp <- genind2genpop(gen)

# Between-population genetic distance (Nei)
Dpop <- dist.genpop(gp, method = 1)

Dmat <- as.matrix(Dpop)
Dmat[upper.tri(Dmat)] <- NA

regions <- factor(unique(pop))

dfD <- as.data.frame(Dmat) %>%
  tibble::rownames_to_column("Pop1") %>%
  tidyr::pivot_longer(-Pop1, names_to = "Pop2", values_to = "NeiD") %>%
  dplyr::filter(!is.na(NeiD)) %>%
  dplyr::mutate(
    # Create factors with all loci in both axes
    # Y-axis reversed for diagonal
    Locus_1 = factor(Pop1, levels = rev(regions)),
    # X-axis in original order
    Locus_2 = factor(Pop2, levels = regions)
  )

pdf(file = paste0(plot_dir, "/NeI_Distance_Btw_Regions.pdf"), width = 3.5, height = 3.5)
ggplot(dfD, aes(x = Locus_2, y = Locus_1, fill = NeiD)) +
  geom_tile(color="white", linewidth=0.3) +
  # scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.82, 0.75),
    legend.background = element_rect(
      fill = alpha("white", 0.75),
      color = "white",
      linewidth = 0.3)) +
  labs(title="Nei distance between regions", x=NULL, y=NULL)
dev.off()


Dpop <- dist.genpop(gp, method = 1)

mds <- cmdscale(as.dist(Dpop), k = 2, eig = TRUE)

mds_df <- tibble(
  Pop = rownames(mds$points),
  Dim1 = mds$points[, 1],
  Dim2 = mds$points[, 2]
)

eig_pos <- mds$eig[mds$eig > 0]
var_exp <- round(100 * eig_pos / sum(eig_pos), 1)

pdf(file = paste0(plot_dir, "/Regions_Nei_Dist_MDS.pdf"), width = 4.5, height = 4.5)
ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Pop)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
  # scale_color_manual(values = c("black", "red"), guide = "none") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
  labs(
    title = "MDS of Regions",
    x = paste0("MDS1 (", var_exp[1], "%)"),
    y = paste0("MDS2 (", var_exp[2], "%)")
  )
dev.off()

# Run a bootsrap and visualize
boot_nei_loci <- function(gen, B = 500, seed = NULL, method = 1) {
  if (!is.null(seed)) set.seed(seed)
  
  pop(gen) <- droplevels(pop(gen))
  pops <- levels(pop(gen))
  comb <- t(combn(pops, 2))
  colnames(comb) <- c("A", "B")
  pair_names <- apply(comb, 1, paste, collapse = "_")
  
  # split into loci
  loc_list <- adegenet::seploc(gen)
  nL <- length(loc_list)
  
  # matrix: loci x pairs (per-locus distances)
  per_locus <- matrix(NA_real_, nrow = nL, ncol = nrow(comb),
                      dimnames = list(names(loc_list), pair_names))
  
  for (i in seq_along(loc_list)) {
    gp_i <- genind2genpop(loc_list[[i]])
    D_i  <- as.matrix(dist.genpop(gp_i, method = method))  # method=1 => Nei
    
    for (k in seq_len(nrow(comb))) {
      a <- comb[k, 1]; b <- comb[k, 2]
      per_locus[i, k] <- D_i[a, b]
    }
  }
  
  # bootstrap: resample loci indices with replacement, average distances
  out <- matrix(NA_real_, nrow = B, ncol = ncol(per_locus))
  colnames(out) <- colnames(per_locus)
  
  for (b in seq_len(B)) {
    idx <- sample.int(nL, size = nL, replace = TRUE)
    out[b, ] <- colMeans(per_locus[idx, , drop = FALSE], na.rm = TRUE)
  }
  
  list(
    boot = as.data.frame(out),
    per_locus = as.data.frame(per_locus),
    n_loci = nL,
    pairs = pair_names
  )
}

# run (Nei distance)
res <- boot_nei_loci(gen, B = 500, seed = 42, method = 1)
boot_df <- res$boot

summary(boot_df)

# Visualize
boot_long <- res$boot %>%
  pivot_longer(cols = everything(), names_to = "Pair", values_to = "NeiD")

pdf(file = paste0(plot_dir, "/Bootstrap_Nei_Dist_Violin.pdf"), width = 5, height = 5)
ggplot(boot_long, aes(Pair, NeiD, fill = Pair)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        legend.position = "none") +
  labs(title = "Bootstrap Nei distance (locus resampling)",
       x = NULL, y = "Nei distance")
dev.off()


sum_df <- boot_long %>%
  group_by(Pair) %>%
  summarise(
    median = median(NeiD, na.rm = TRUE),
    lo = quantile(NeiD, 0.025, na.rm = TRUE),
    hi = quantile(NeiD, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

pdf(file = paste0(plot_dir, "/Bootstrap_CI_Nei_Dist.pdf"), width = 5, height = 5)
ggplot(sum_df, aes(Pair, median)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
  labs(title = "Nei distance: median and 95% bootstrap CI",
       x = NULL, y = "Nei distance")
dev.off()

perloc_long <- res$per_locus %>%
  tibble::rownames_to_column("Locus") %>%
  pivot_longer(-Locus, names_to = "Pair", values_to = "NeiD")

pdf(file = paste0(plot_dir, "/Bootstrap_per_locus_nei_drivers_heatmap.pdf"), width = 5, height = 5)
ggplot(perloc_long, aes(Pair, Locus, fill = NeiD)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
  labs(title = "Per-locus Nei distance Drivers",
       x = NULL, y = NULL)
dev.off()

###############################################################################
# Population-Comparisons
# Download the (current) STRidER database
freqs <- read_STRidER_xml()

# See what populations exist
names(freqs)

# Grab the countries you care about
pops <- c("GERMANY", "FRANCE", "GREECE", "MONTENEGRO", "SLOVENIA", "HUNGARY", 
          "BOSNIA AND HERZEGOWINA")#, "POLAND")

# Some country labels may differ slightly; check names(freqs) and adjust
freqs_sub <- freqs[pops]

# list loci available for one population
freqs_sub[1]

# look at allele frequencies at TH01 for Germany
freqs_sub$GERMANY$TH01

##################################################################################
# Create Allele Frequencies for Albania
alb_path <- "albanian_data_by_regions.txt"

alb <- read.delim(alb_path, check.names = FALSE, stringsAsFactors = FALSE)

# keep only ALBANIA rows
north <- alb[alb$pop == "NORTH", ]
south <- alb[alb$pop == "SOUTH", ]
center <- alb[alb$pop == "CENTER", ]

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

north_freq_list <- setNames(lapply(loci_alb, function(L) locus_freq(north, L)), loci_alb)
south_freq_list <- setNames(lapply(loci_alb, function(L) locus_freq(south, L)), loci_alb)
center_freq_list <- setNames(lapply(loci_alb, function(L) locus_freq(center, L)), loci_alb)

freqs_sub$ALBANIA_NORTH <- north_freq_list
freqs_sub$ALBANIA_SOUTH <- south_freq_list
freqs_sub$ALBANIA_CENTER <- center_freq_list

# Shared loci across all populations (critical)
loci_per_pop <- lapply(freqs_sub, names)
common_loci <- Reduce(intersect, loci_per_pop)
stopifnot(length(common_loci) >= 2)

freqs_shared <- lapply(freqs_sub, function(pop) pop[common_loci])


# Let's see what we have
loci_per_pop <- lapply(freqs_sub, names)

# look at first few
loci_tbl <- tibble(
  Population = names(loci_per_pop),
  n_loci = sapply(loci_per_pop, length))

loci_tbl

all_loci <- Reduce(union, loci_per_pop)

missing_loci <- lapply(loci_per_pop, function(x) setdiff(all_loci, x))

missing_count <- sapply(missing_loci, length)
missing_count

length(common_loci)
common_loci

# Helpers: extract allele freqs + Nei distances + D matrix
# Convert any "locus frequency object" into (alleles, freqs)
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

nei_one_locus <- function(A, B){
  A <- af_extract(A); B <- af_extract(B)
  alleles <- union(A$alleles, B$alleles)
  p <- setNames(rep(0, length(alleles)), alleles)
  q <- p
  p[A$alleles] <- A$freqs
  q[B$alleles] <- B$freqs
  I <- sum(sqrt(p * q))
  I <- max(I, 1e-12)
  -log(I)
}

build_Dmat <- function(freqs_shared, loci_use){
  pops <- names(freqs_shared)
  n <- length(pops)
  D <- matrix(0, n, n, dimnames = list(pops, pops))
  
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      per_locus <- vapply(loci_use, function(loc){
        nei_one_locus(freqs_shared[[pops[i]]][[loc]], freqs_shared[[pops[j]]][[loc]])
      }, numeric(1))
      d <- mean(per_locus, na.rm = TRUE)
      D[i, j] <- d
      D[j, i] <- d
    }
  }
  D
}

Dmat <- build_Dmat(freqs_shared, common_loci)

write.csv(Dmat, file = paste0(out_dir, "/Between_populations_distance.csv"))

# Neighbor-Joining tree
tree <- nj(as.dist(Dmat))

pdf(file = paste0(plot_dir, "/Populations_NJ_Tree.pdf"), width = 6, height = 3.5)
plot(tree,
     main = "NJ tree (Nei distance)",
     use.edge.length = FALSE,
     cex = 0.7,
     font = 2,
     no.margin = TRUE)
dev.off()


tree_um <- chronos(tree)
pdf(file = paste0(plot_dir, "/Populations_NJ_Tree_Ultrametricized.pdf"), width = 6, height = 3.5)
plot(tree_um,
     main = "Ultrametricized NJ tree",
     use.edge.length = TRUE,
     cex = 0.7,
     font=2)
dev.off()

mds <- cmdscale(as.dist(Dmat), k = 2, eig = TRUE)
mds_df <- tibble(
  Pop = rownames(mds$points),
  Dim1 = mds$points[, 1],
  Dim2 = mds$points[, 2]
)

eig_pos <- mds$eig[mds$eig > 0]
var_exp <- round(100 * eig_pos / sum(eig_pos), 1)

pdf(file = paste0(plot_dir, "/Populations_Nei_Dist_MDS.pdf"), width = 5.5, height = 5)
ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Pop)) +
  geom_point(aes(color = Pop == "ALBANIA"), size = 3) +
  geom_text_repel(size = 4, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  theme_classic(base_size = 13) +
  labs(
    title = "MDS of STR populations",
    x = paste0("MDS1 (", var_exp[1], "%)"),
    y = paste0("MDS2 (", var_exp[2], "%)")
  )
dev.off()

# Bootstrap
bootstrap_nj_loci <- function(freqs_shared, B = 500, seed = 42, outgroup = NULL, return_boot_trees = FALSE){
  set.seed(seed)
  
  loci_per_pop <- lapply(freqs_shared, names)
  common_loci <- Reduce(intersect, loci_per_pop)
  if (length(common_loci) < 2) stop("Too few loci shared across populations.")
  
  D0 <- build_Dmat(freqs_shared, common_loci)
  ref <- nj(as.dist(D0))
  if (!is.null(outgroup) && outgroup %in% ref$tip.label) {
    ref <- root(ref, outgroup = outgroup, resolve.root = TRUE)
  }
  ref <- ladderize(ref)
  
  L <- length(common_loci)
  boot_trees <- vector("list", B)
  
  for (b in seq_len(B)) {
    loci_b <- sample(common_loci, size = L, replace = TRUE)
    Db <- build_Dmat(freqs_shared, loci_b)
    tb <- nj(as.dist(Db))
    if (!is.null(outgroup) && outgroup %in% tb$tip.label) {
      tb <- root(tb, outgroup = outgroup, resolve.root = TRUE)
    }
    boot_trees[[b]] <- ladderize(tb)
  }
  
  counts <- prop.clades(ref, boot_trees)
  support <- 100 * counts / B
  
  list(
    tree = ref,
    support = support,
    common_loci = common_loci,
    Dmat = D0,
    boot_trees = if (return_boot_trees) boot_trees else NULL
  )
}

res <- bootstrap_nj_loci(freqs_shared, B = 500, seed = 42, outgroup = NULL)

plot(res$tree, main = "NJ tree with locus-bootstrap support", cex = 0.9)
nodelabels(round(res$support, 1), frame = "none", cex = 0.8)

pdf(file = paste0(plot_dir, "/Populations_NJ_Tree_500_Bootsrap.pdf"), width = 6, height = 3.5)
plot(res$tree,
     main = "NJ tree (Nei distance)",
     use.edge.length = FALSE,
     cex = 0.7,
     font = 2,
     no.margin = TRUE)
dev.off()

tree_um <- chronos(res$tree)
pdf(file = paste0(plot_dir, "/Populations_NJ_Tree_Ultrametricized_500_Bootstrap.pdf"), width = 6, height = 3.5)
plot(tree_um,
     main = "Ultrametricized NJ tree",
     use.edge.length = TRUE,
     cex = 0.7,
     font=2)
dev.off()

# Distance heatmap (lower triangle) for the bootstrap reference Dmat
D_keep <- res$Dmat
D_keep[upper.tri(D_keep, diag = TRUE)] <- NA
labs_idx <- colnames(D_keep)

dm_long <- as.data.frame(D_keep) %>%
  rownames_to_column("Pop_1") %>%
  pivot_longer(-Pop_1, names_to = "Pop_2", values_to = "Dist") %>%
  filter(!is.na(Dist)) %>%
  mutate(
    Pop_1 = factor(Pop_1, levels = rev(labs_idx)),
    Pop_2 = factor(Pop_2, levels = labs_idx)
  )

pdf(file = paste0(plot_dir, "/Populations_NEI_Distance_Heatmap.pdf"), width = 5, height = 6)
ggplot(dm_long, aes(x = Pop_2, y = Pop_1, fill = Dist)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = 1, name = "Nei D") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.78, 0.78),
        legend.background = element_rect(
          fill = alpha("white", 0.75),
          color = "white",
          linewidth = 0.3)) +
  labs(title = "Nei distance", x = NULL, y = NULL) +
  coord_fixed() +
  scale_x_discrete(drop = FALSE, limits = labs_idx) +
  scale_y_discrete(drop = FALSE, limits = rev(labs_idx))
dev.off()

write.csv(res$Dmat, file = paste0(out_dir, "/Between_populations_distance_500_bootsrap.csv"))

cluster_data <- read.csv("ld_on_pca_unsupervised_clustering.csv")
head(cluster_data)

cluster_regions <- inner_join(regions, cluster_data, by = c("ind" = "Ind"))

table(cluster_regions$City, cluster_regions$Cluster)

pdf(file = paste0(plot_dir, "/Cluster_per_city_barplot.pdf"), width = 6, height = 5)
ggplot(cluster_regions, aes(x = City, fill = factor(Cluster))) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("steelblue", "firebrick", "darkorange")) +
  labs(
    x = "City",
    y = "Number of People",
    fill = "Cluster"
  ) +
  scale_y_continuous(trans = "log2")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
