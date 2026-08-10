# Prepared by Ilir Sheraj (Dr. Gjilpera)
# This script summarizes STRUCTURE outputs and generates images to interpret the data.
# FOr details of STRUCTURE, see the associated folder
library(pophelper)
library(dplyr)
library(tidyr)
library(ggplot2)

# Folder containing all 120 STRUCTURE output files
structure_dir <- "Albania_structure_full/primary_output"

out_dir   <- file.path(getwd(), "STRUCTURE_Tables")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)

structure_files <- list.files(
  structure_dir,
  pattern = "_f$",
  full.names = TRUE)

length(structure_files)

# Read STRUCTURE output files
qlist <- readQ(files = structure_files, filetype = "structure")

# Check that all files were imported
length(qlist)

# Extract all the summariezd values
run_table <- tabulateQ(qlist)
head(run_table)
# Make sure there are 20 per K
table(run_table$k)

# Summarize 
k_summary <- summariseQ(run_table)

k_summary
names(k_summary)

# Save it
write.csv(run_table, file = paste0(out_dir, "/STRUCTURE_run_table.csv"), 
          row.names = FALSE)
write.csv(k_summary, file = paste0(out_dir, "/STRUCTURE_K_summary.csv"), 
          row.names = FALSE)

# Single Plot, this is difficult to read
evanno_results <- evannoMethodStructure(
  data = k_summary,
  exportplot = TRUE,
  exportpath = out_dir,
  outputfilename = "STRUCTURE_Evanno")

evanno_results
names(evanno_results)

# Let's make it better
evanno_results <- evannoMethodStructure(
  data = k_summary,
  exportplot = TRUE,
  exportpath = out_dir,
  outputfilename = "STRUCTURE_Evanno",
  imgtype = "pdf",
  width = 18,
  height = 18,
  units = "cm",
  basesize = 8
)

# Extract all the values and create separate plots
evanno_table <- evannoMethodStructure(
  data = k_summary,
  returndata = TRUE,
  returnplot = FALSE,
  exportplot = FALSE)

str(evanno_table)
names(evanno_table)
evanno_table
write.csv(evanno_table, file = paste0(out_dir, "/evanno_table.csv"), 
          row.names = FALSE)

evanno_output <- evannoMethodStructure(
  data = k_summary,
  returndata = TRUE,
  returnplot = TRUE,
  exportplot = FALSE)

evanno_table <- evanno_output$data
evanno_plot  <- evanno_output$plot

names(evanno_output)

evanno_table <- evanno_output[vapply(evanno_output, is.data.frame, logical(1))][[1]]
names(evanno_table)

# Make a copy jus tin case you screw up
evanno_df <- evanno_table

#######################################################################
# First Set of Plots
# Mean estimated log probability
plot_A <- ggplot(evanno_df, aes(x = k, y = elpdmean)) +
  geom_errorbar(aes(ymin = elpdmin, ymax = elpdmax),
                width = 0.12, colour = "grey35") +
  geom_line(colour = "#377EB8", linewidth = 0.7) +
  geom_point(colour = "#377EB8", size = 2.5) +
  scale_x_continuous(breaks = evanno_df$k) +
  geom_errorbar(aes(
      ymin = elpdmean - elpdsd,
      ymax = elpdmean + elpdsd),
    width = 0.12,
    colour = "grey35") +
  labs(x = "Number of clusters (K)",
       y = expression("Mean ln P(D)")) +
  theme_classic(base_size = 11)

ggsave("Mean_Estimated_Log_Probability.pdf", 
       plot = plot_A, 
       width = 4, 
       height = 4, 
       units = "in",
       path = out_dir)

ggsave(
  "Mean_Estimated_Log_Probability.jpeg", 
  plot = plot_A, 
  width = 3, 
  height = 3, 
  units = "in",
  dpi = 600,
  path = out_dir)

# First-order rate of change
plot_B <- ggplot(subset(evanno_df, !is.na(lnk1)), aes(x = k, y = lnk1)) +
  geom_errorbar(aes(ymin = lnk1min, ymax = lnk1max),
    width = 0.12, colour = "grey35") +
  geom_line(colour = "#377EB8", linewidth = 0.7) +
  geom_point(colour = "#377EB8", size = 2.5) +
  scale_x_continuous(breaks = evanno_df$k) +
  labs(x = "Number of clusters (K)", y = expression(L^minute(K))) +
  theme_classic(base_size = 11)

ggsave("First_Order_Rate_of_Change.pdf", 
       plot = plot_B, 
       width = 4, 
       height = 4, 
       units = "in", 
       path = out_dir)

ggsave("First_Order_Rate_of_Change.jpeg",
       plot = plot_B, 
       width = 3, 
       height = 3, 
       units = "in",
       dpi = 600,
       path = out_dir)

# Absolute Second Order Rate of Change
plot_C <- ggplot(subset(evanno_df, !is.na(lnk2)), aes(x = k, y = lnk2)) +
  geom_errorbar(aes(ymin = lnk2min, ymax = lnk2max),
                width = 0.12, colour = "grey35") +
  geom_line(colour = "#377EB8", linewidth = 0.7) +
  geom_point(colour = "#377EB8", size = 2.5) +
  scale_x_continuous(breaks = evanno_df$k) +
  labs(x = "Number of clusters (K)",
       y = expression(abs(L^second(K)))) +
  theme_classic(base_size = 11)

ggsave("Absolute_Second_Order_Rate_of_Change.pdf", 
       plot = plot_C, 
       width = 4, 
       height = 4, 
       units = "in",
       path = out_dir)

ggsave("Absolute_Second_Order_Rate_of_Change.jpeg", 
       plot = plot_C, 
       width = 3, 
       height = 3, 
       units = "in",
       dpi = 600,
       path = out_dir)

# Evanno Delta_K
plot_D <- ggplot(subset(evanno_df, !is.na(deltaK)), aes(x = k, y = deltaK)) +
  geom_line(colour = "#377EB8", linewidth = 0.7) +
  geom_point(colour = "#377EB8", size = 2.5) +
  scale_x_continuous(breaks = evanno_df$k) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of clusters (K)",
       y = expression(Delta*K)) +
  theme_classic(base_size = 11)

ggsave("Evanno_Delta_K.pdf", 
       plot = plot_D, 
       width = 4, 
       height = 4, 
       units = "in",
       path = out_dir)

ggsave("Evanno_Delta_K.jpeg", 
       plot = plot_D, 
       width = 3, 
       height = 3, 
       units = "in",
       dpi = 600,
       path = out_dir)

evanno_table[, c("k", "elpdmean", "elpdsd", "lnk1", "lnk2", "deltaK")]

###################################################
# Best Supported Structure Model
best_model <- ggplot(run_table, aes(x = factor(k), y = elpd)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.65, size = 1.7) +
  labs(x = "Number of clusters (K)",
       y = "Estimated log Probability") +
  theme_classic()

ggsave("best_structure_model.pdf", 
       plot = best_model, 
       width = 5, 
       height = 4, 
       units = "in",
       path = out_dir)

ggsave("best_structure_model.jpeg", 
       plot = best_model, 
       width = 4, 
       height = 3, 
       units = "in",
       dpi = 600,
       path = out_dir)


# Create the heterogeneity plot
# Create a consensus Q matrix for each K
qlist_aligned <- alignK(qlist)

consensus_q <- lapply(2:6, function(selected_k) {
  
  runs_k <- qlist_aligned[
    vapply(qlist_aligned, ncol, integer(1)) == selected_k]
  
  if (length(runs_k) != 20) {
    warning("K = ", selected_k, " contains ", length(runs_k), " runs rather than 20.")
  }
  
  # Element-wise average across aligned replicates
  q_array <- simplify2array(lapply(runs_k, as.matrix))
  
  q_mean <- apply(q_array, c(1, 2), mean)
  
  # Correct tiny rounding deviations so rows sum exactly to 1
  q_mean <- q_mean / rowSums(q_mean)
  
  colnames(q_mean) <- paste0("Cluster_", seq_len(selected_k))
  
  q_mean
})

names(consensus_q) <- paste0("K", 2:6)

# Make a long dataframe for plotitng
q_long <- bind_rows(
  lapply(2:6, function(selected_k) {
    qmat <- consensus_q[[paste0("K", selected_k)]]
    
    as.data.frame(qmat) %>%
      mutate(
        Individual = seq_len(nrow(qmat)),
        K = paste0("K = ", selected_k)
      ) %>%
      pivot_longer(
        cols = starts_with("Cluster_"),
        names_to = "Cluster",
        values_to = "Q")
    })
  )

# Keep K panels in numerical order
q_long$K <- factor(q_long$K, levels = paste0("K = ", 2:6))

# Define 6 colors that are discernable from each other
cluster_colours <- c(
  "Cluster_1" = "#3B6FB6",
  "Cluster_2" = "#D95F59",
  "Cluster_3" = "#58A65C",
  "Cluster_4" = "#8C6BB1",
  "Cluster_5" = "#E6A13C",
  "Cluster_6" = "#56B4C2")

panel_a <- ggplot(q_long, aes(x = Individual, y = Q, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(K ~ ., switch = "y") +
  scale_fill_manual(
    values = cluster_colours,
    drop = FALSE) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(1, 500, 1000, 1500, 2000),
    labels = c("1", "0.5k", "1k", "1.5k", "2k")) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5)) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")) +
  labs(
    x = "Individuals",
    y = "Membership coefficient (Q)",
    fill = "Cluster") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    
    strip.text.y.left = element_text(
      size = 7.5,
      angle = 0,
      face = "bold",
      margin = margin(r = 4)),
    
    strip.background = element_blank(),
    # Increase separation between K panels
    panel.spacing.y = unit(0.22, "cm"),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 8)),
    
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    legend.position = "bottom")

ggsave("STRUCTURE_discrete_clusters.pdf", 
       plot = panel_a, 
       width = 6, 
       height = 5, 
       units = "in",
       path = out_dir)

# Modify to fit less space
panel_a <- ggplot(q_long, aes(x = Individual, y = Q, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(K ~ ., switch = "y") +
  scale_fill_manual(
    values = cluster_colours,
    drop = FALSE) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(1, 500, 1000, 1500, 2000),
    # use k for labels, shorter than 000
    labels = c("1", "0.5k", "1k", "1.5k", "2k")) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5)) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")) +
  labs(
    x = "Individuals",
    y = "Membership coefficient (Q)",
    fill = "Cluster"
  ) +
  theme_classic(base_size = 9) +
  scale_fill_manual(
    values = cluster_colours,
    breaks = paste0("Cluster_", 1:6),
    labels = 1:6,
    drop = FALSE) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE)) +
  
  theme(
    axis.text.y = element_text(size = 7),
    
    strip.text.y.left = element_text(
      size = 7.5,
      angle = 0,
      face = "bold",
      margin = margin(r = 4)),
    
    strip.background = element_blank(),
    panel.spacing.y = unit(0.22, "cm"),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 8)),
    
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.15, "cm"))

ggsave("STRUCTURE_discrete_clusters.jpeg", 
       plot = panel_a, 
       width = 4, 
       height = 3, 
       units = "in",
       dpi = 600,
       path = out_dir)

ggsave("STRUCTURE_discrete_clusters_new.pdf", 
       plot = panel_a, 
       width = 6, 
       height = 5, 
       units = "in",
       path = out_dir)

# Violin Plot for maximum membership
max_q_df <- bind_rows(
  lapply(2:6, function(selected_k) {
    
    qmat <- consensus_q[[paste0("K", selected_k)]]
    
    data.frame(
      K_number = selected_k,
      K = paste0("K = ", selected_k),
      Individual = seq_len(nrow(qmat)),
      Max_Q = apply(qmat, 1, max),
      Uniform_Q = 1 / selected_k)
  })
)

max_q_df$K <- factor(max_q_df$K, levels = paste0("K = ", 2:6))

uniform_lines <- max_q_df %>%
  distinct(K, Uniform_Q)

panel_b <- ggplot(max_q_df, aes(x = K, y = Max_Q)) +
  geom_violin(
    fill = "grey85",
    colour = "grey35",
    width = 0.8,
    scale = "width") +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white") +
  geom_segment(
    data = uniform_lines,
    aes(
      x = as.numeric(K) - 0.38,
      xend = as.numeric(K) + 0.38,
      y = Uniform_Q,
      yend = Uniform_Q),
    inherit.aes = FALSE,
    colour = "#C23B3B",
    linewidth = 0.8,
    linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(
    labels = c(
      "K = 2" = "2",
      "K = 3" = "3",
      "K = 4" = "4",
      "K = 5" = "5",
      "K = 6" = "6")) +
  
  labs(x = "Specified number of clusters (K)",
       y = "Max Membership Coeff per Individual (Qmax)") +
  theme_classic(base_size = 10)

ggsave("STRUCTURE_maximum_membership_coefficient_violin.pdf", 
       plot = panel_b, 
       width = 6, 
       height = 5, 
       units = "in",
       path = out_dir)

ggsave("STRUCTURE_maximum_membership_coefficient_violin.jpeg", 
       plot = panel_b, 
       width = 4, 
       height = 3, 
       units = "in",
       dpi = 600,
       path = out_dir)
