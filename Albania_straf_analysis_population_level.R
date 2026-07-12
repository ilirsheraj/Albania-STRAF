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

# Working Here: "/home/ilir/EZaimi"
# Make two subdirectories
# ## One to save csv tables
# out_dir   <- file.path(getwd(), "Pop_Italy_outputs")
# ## Second, to save the figures
# plot_dir  <- file.path(getwd(), "Pop_Italy_plots")
## One to save csv tables
out_dir   <- file.path(getwd(), "Pop_Italy_Slovakia_outputs")
## Second, to save the figures
plot_dir  <- file.path(getwd(), "Pop_Italy_Slovakia_plots")

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Part 1: Data Preparation
dat <- read.table(
  "STRAF 2000 profile Full name Loci.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat[1:5, 1:6]
colnames(dat)

ind <- dat$ind
pop <- dat$pop

# Convert to Diploid Frequency Table
geno2 <- dat[, -(1:2)]

# One Column per Locus
loci <- colnames(geno2)[seq(1, ncol(geno2), by = 2)]
length(loci)

# Convert them to diploid frequencies table (Allele1/Allele2)
geno1 <- as.data.frame(
  sapply(seq_along(loci), function(i) {
    paste(geno2[[2*i - 1]], geno2[[2*i]], sep = "/")}),
  stringsAsFactors = FALSE
)

colnames(geno1) <- loci
rownames(geno1) <- ind


# Part 2: Building the genind object
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

# Names of Loci
locNames(gen)


# Part 3: EDA of geneInd
# Visualize allele frequencies for each locus (Histrograms)
af <- makefreq(gen, missing = "mean")

# Convert the matrix into a dataframe
af_df <- as.data.frame(af) %>%
  rownames_to_column("ID") %>% 
  pivot_longer(-ID, names_to = "AlleleID", values_to = "Freq") %>%
  filter(!is.na(Freq) & Freq > 0) %>%
  separate(AlleleID, into = c("Locus", "Allele"), sep = "\\.", extra = "merge") %>%
  mutate(Allele = gsub("_", ".", Allele))

pdf(file = paste0(plot_dir, "/Marker_Allele_Frequency.pdf"), width = 6, height = 6)
ggplot(af_df, aes(x = as.numeric(Allele), y = Freq)) +
  geom_col(width = 0.5, fill = "brown1") +
  facet_wrap(~ Locus, scales = "free_x") +
  theme_classic() +
  labs(title = "Allele frequencies per locus", 
       x = "Alleles", 
       y = "Frequency") +
  theme(axis.text.x = element_text(size = 7),
        strip.text = element_text(size = 9),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()

# Part 4: Summary stats (Ho, Hs, Fis, etc.)
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
perloc_df <- perloc_df[, !names(perloc_df) %in% c("Dstp", "Htp", "Fstp", "Dest")]
head(perloc_df)

write.csv(perloc_df, file = paste0(out_dir, "/Summary_statistics.csv"))

# Part 6: PCA & MDS
X <- scaleGen(gen, NA.method = "mean")
pca <- prcomp(X)

pve_percent <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

# Use ggplot to plot PCA
pdf(file = paste0(plot_dir, "/Albanian_Population_PCA.pdf"), width = 4, height = 4)
as.data.frame(pca$x[, 1:2]) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(color = "brown1") +
  labs(title = "Albanian Population PCA",
       x = paste0("PC1 (", pve_percent[1], "%)"),
       y = paste0("PC2 (", pve_percent[2], "%)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()

# # Classical Plot
# plot(pca$x[,1], pca$x[,2],
#      col = as.numeric(pop(gen)),
#      pch = 19,
#      xlab = "PC1", ylab = "PC2")
# legend("topright", legend = levels(pop(gen)),
#        col = 1:nlevels(pop(gen)), pch = 19)


# MDS: Make it Stable, otherwise it sucks
replen <- rep(4, nLoc(gen))
names(replen) <- locNames(gen)
if ("D22S1045" %in% names(replen)) replen["D22S1045"] <- 3

# Check 
length(replen)
nLoc(gen)
all(names(replen) == locNames(gen))

# At this point it is for the entire population, its kinda meaningless
d_ind <- bruvo.dist(gen, replen = replen)

mds <- cmdscale(d_ind, k = 2, eig = TRUE)

mds_df <- data.frame(
  Pop = rownames(mds$points),
  Dim1 = mds$points[,1],
  Dim2 = mds$points[,2]
)
var_exp <- round(100 * mds$eig / sum(mds$eig[mds$eig > 0]), 1)

# Same here: Use ggplot
pdf(file = paste0(plot_dir, "/Albanian_Population_MDS.pdf"), width = 4, height = 4)
mds_df %>% ggplot(aes(Dim1, Dim2)) +
  geom_point(color = "brown1") +
  labs(title = "Albanian Population MDS",
       x = paste0("MDS1 (", var_exp[1], "%)"),
       y = paste0("MDS2 (", var_exp[2], "%)")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))
dev.off()


# Pairwise Population Differentiation: Here this is uselesss
# hf <- genind2hierfstat(gen)
# pairwise.WCfst(hf)

# Part 7: Hardy-Weinberg Equilibrium per locus
# HWE per locus (Monte Carlo / permutation based)
# B = number of permutations (bigger = more stable p-values)
hwe_res <- pegas::hw.test(gen, B = 10000)

hwe_res

hwe_df <- as.data.frame(hwe_res) %>%
  rownames_to_column("Locus") %>%
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
# Linkage Disequilibrium (LD) within Markers
# library(graph4lg)

# A bunch of function to deal with problematic names: Pain in the ass
allele_to3 <- function(a) {
  # Convert allele like "9.3" -> "093" (3-digit), "32.2" -> "322"
  a <- as.character(a)
  a <- gsub("_", ".", a)
  a <- gsub("\\s+", "", a)
  
  if (is.na(a) || a == "" || a == "NA") return(NA_character_)
  
  # remove dot for microvariants: 9.3 -> 93, 32.2 -> 322
  a_num <- gsub("\\.", "", a)
  
  # pad to 3 digits if possible
  if (!grepl("^[0-9]+$", a_num)) return(NA_character_)
  sprintf("%03d", as.integer(a_num))
}

geno_to_genepop_code <- function(g) {
  # g like "9.3/9.3" or "28.3/32.2" etc
  if (is.na(g) || g == "" || g == "NA") return("000000")  # missing
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return("000000")
  
  a1 <- allele_to3(parts[1])
  a2 <- allele_to3(parts[2])
  if (is.na(a1) || is.na(a2)) return("000000")
  paste0(a1, a2)
}

# Apply to your geno1 data.frame (one column per locus, "a/b" strings)
gp_geno <- as.data.frame(lapply(geno1, function(col) vapply(col, geno_to_genepop_code, "")),
                         stringsAsFactors = FALSE)

head(gp_geno)


# write Genepop file with .txt extensions
gp_tmp <- tempfile(fileext = ".txt")
out_tmp <- tempfile(fileext = ".txt")

con <- file(gp_tmp, open = "wt")
writeLines("Genepop file exported from R", con)
writeLines(paste(colnames(gp_geno), collapse = ", "), con)
writeLines("Pop", con)

for (i in seq_len(nrow(gp_geno))) {
  id <- ind[i]
  id <- gsub("[^A-Za-z0-9_\\-]", "_", id)
  gline <- paste(gp_geno[i, ], collapse = " ")
  writeLines(paste0(id, " , ", gline), con)
}
close(con)

cat(readLines(gp_tmp, n = 6), sep = "\n")
file.exists(gp_tmp)

# Now calculate: 10000 permutations here
genepop::test_LD(
  inputFile = gp_tmp,
  outputFile = out_tmp,
  dememorization = 5000,
  batches = 100,
  iterations = 10000
)

ld_df <- straf::get_ld_gp(out_tmp)
head(ld_df)

# Convert the data frame into a matrix
loci <- colnames(gp_geno)

ld_mat <- matrix(NA_real_, length(loci), length(loci),
                 dimnames = list(loci, loci))
diag(ld_mat) <- 1

for (i in seq_len(nrow(ld_df))) {
  a <- as.character(ld_df$Locus_1[i])
  b <- as.character(ld_df$Locus_2[i])
  p <- ld_df$P_value[i]
  ld_mat[a, b] <- p
  ld_mat[b, a] <- p
}

ld_mat

write.csv(ld_mat, file = paste0(out_dir, "/Linkage_Disequilibrium_Matrix_10k_Permutations.csv"))

# Make it ready for visualization# Full Matrix
ld_long <- as.data.frame(ld_mat) %>%
  rownames_to_column("Locus_1") %>%
  pivot_longer(-Locus_1, names_to = "Locus_2", values_to = "P_value")

ld_long$Locus_1 <- factor(ld_long$Locus_1, levels = rev(colnames(ld_mat)))
ld_long$Locus_2 <- factor(ld_long$Locus_2, levels = colnames(ld_mat))

ggplot(ld_long, aes(Locus_1, Locus_2, fill = P_value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
  labs(
    title = "Linkage Disequilibrium Between Loci",
    fill = expression(pval),
    x = NULL, y = NULL
  )

# Triangle Matrix
# Keep lower triangle: Remove the diagonal
ld_mat_keep <- ld_mat
ld_mat_keep[upper.tri(ld_mat_keep, diag = TRUE)] <- NA

# Direct approach using your original lower triangular matrix
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

# Create the plot
pdf(file = paste0(plot_dir, "/Linkage_Disequilibrium_10k_Permutations.pdf"), width = 4, height = 4)
ggplot(ld_long_keep, aes(x = Locus_2, y = Locus_1, fill = -log10(P_value))) +
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
    title = "Linkage Disequilibrium Between Loci", 
    fill = "p", 
    x = NULL, 
    y = NULL
  ) +
  coord_fixed() +
  # Ensure all loci appear on both axes
  scale_x_discrete(drop = FALSE, limits = loci) +
  scale_y_discrete(drop = FALSE, limits = rev(loci))
dev.off()

##################################################################################
# Create Phylogenetic Tree: For the moment its bullshit because its on individual level
# tree_ind <- nj(as.dist(d_ind))
# 
# plot(tree_ind, cex = 0.7)
# tiplabels(pch = 19, col = as.numeric(pop(gen)), cex = 0.6)
# legend("topleft", legend = levels(pop(gen)),
#        col = 1:nlevels(pop(gen)), pch = 19, cex = 0.8)

###############################################################################
# Population-Comparisons
# Download the (current) STRidER database
freqs <- read_STRidER_xml()

# See what populations exist
names(freqs)

# Grab the countries you care about
pops <- c("GERMANY", "FRANCE", "GREECE", "MONTENEGRO", "SLOVENIA", "HUNGARY", 
          "BOSNIA AND HERZEGOWINA", "ITALY", "SLOVAKIA")#, "POLAND")

# Some country labels may differ slightly; check names(freqs) and adjust
freqs_sub <- freqs[pops]

# list loci available for one population
freqs_sub[1]

# look at allele frequencies at TH01 for Germany
freqs_sub$GERMANY$TH01

##################################################################################
# Create Allele Frequencies for Albania
alb_path <- "STRAF 2000 profile Full name Loci.txt"

alb <- read.delim(alb_path, check.names = FALSE, stringsAsFactors = FALSE)

# keep only ALBANIA rows (Can change later)
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

alb_freq_list <- setNames(lapply(loci_alb, function(L) locus_freq(alb, L)), loci_alb)
freqs_sub$ALBANIA <- alb_freq_list

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

write.csv(loci_tbl, file = paste0(out_dir, "/loci_per_country.csv"))

all_loci <- Reduce(union, loci_per_pop)

missing_loci <- lapply(loci_per_pop, function(x) setdiff(all_loci, x))

missing_count <- sapply(missing_loci, length)
missing_count

length(common_loci)
common_loci

################################################################################
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

# Neighbor-Joining tree
tree <- nj(as.dist(Dmat))

# plot(tree, main = "NJ tree (Nei distance)",
#      use.edge.length = FALSE)

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

################################################################################
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
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
  labs(
    title = "MDS of STR populations",
    x = paste0("MDS1 (", var_exp[1], "%)"),
    y = paste0("MDS2 (", var_exp[2], "%)")
  )
dev.off()

################################################################################
# Not Necessary
tree_nj <- nj(as.dist(Dmat)) %>% ladderize()

# Rooting is optional; only do it if your chosen outgroup is truly “outside”
outgroup <- "HUNGARY"
if (outgroup %in% tree_nj$tip.label) {
  tree_nj_root <- root(tree_nj, outgroup = outgroup, resolve.root = TRUE) %>% ladderize()
} else {
  tree_nj_root <- tree_nj
}

plot(tree_nj_root, main = "NJ tree (Nei distance)", cex = 0.9)
add.scale.bar()

hc <- hclust(as.dist(Dmat), method = "average")
tree_upgma <- as.phylo(hc) %>% ladderize()
if (outgroup %in% tree_upgma$tip.label) {
  tree_upgma <- root(tree_upgma, outgroup = outgroup, resolve.root = TRUE) %>% ladderize()
}
plot(tree_upgma, main = "UPGMA tree (Nei distance)", cex = 0.9)

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

# MDS
Dmat <- res$Dmat
mds <- cmdscale(as.dist(Dmat), k = 2, eig = TRUE)
mds_df <- tibble(
  Pop = rownames(mds$points),
  Dim1 = mds$points[, 1],
  Dim2 = mds$points[, 2]
)

eig_pos <- mds$eig[mds$eig > 0]
var_exp <- round(100 * eig_pos / sum(eig_pos), 1)

pdf(file = paste0(plot_dir, "/Populations_Nei_Dist_MDS_500_bootstrap.pdf"), width = 5.5, height = 5)
ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Pop)) +
  geom_point(aes(color = Pop == "ALBANIA"), size = 3) +
  geom_text_repel(size = 4, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
  labs(
    title = "MDS of STR populations",
    x = paste0("MDS1 (", var_exp[1], "%)"),
    y = paste0("MDS2 (", var_exp[2], "%)")
  )
dev.off()
#########################################
