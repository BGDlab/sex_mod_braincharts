library(tidyverse)

clust_var <- "deriv_M_minus_F_z"

df <- full_join(LR.df_tot, sexdiffs.df.u) %>%
  filter(sig_fdr==TRUE & total==TRUE & weighted==FALSE & cv_sample=="B") %>%
  filter(!(pheno %in% vent_list)) %>%
  select(pheno, logAge_days, clust_var, pheno_cat)

# --- 1. Pivot to wide format: rows = pheno levels, cols = time points ---
ts_wide <- df %>%
  arrange(pheno, logAge_days) %>%
  select(pheno, logAge_days, clust_var) %>%
  pivot_wider(
    names_from  = logAge_days,
    values_from = clust_var
  )

# Store pheno labels, then convert to matrix
pheno_labels <- ts_wide$pheno
ts_matrix <- ts_wide %>%
  select(-pheno) %>%
  as.matrix()

rownames(ts_matrix) <- pheno_labels

# --- 2. Handle missing values (if time grids differ across pheno levels) ---
# Option A: remove time points with any NA
ts_matrix_complete <- ts_matrix[, colSums(is.na(ts_matrix)) == 0]

# Option B: impute NAs with row means (uncomment to use instead)
# ts_matrix_complete <- t(apply(ts_matrix, 1, function(x) {
#   x[is.na(x)] <- mean(x, na.rm = TRUE); x
# }))

# --- 3. Compute Euclidean distance between time series ---
dist_matrix <- dist(ts_matrix_complete, method = "euclidean")

# --- 4. Hierarchical clustering ---
hc <- hclust(dist_matrix, method = "ward.D2")  # try also "complete", "average"

# --- 5. Plot dendrogram ---
plot(hc,
     main  = "Hierarchical Clustering of Pheno Levels",
     xlab  = "Pheno",
     ylab  = "Euclidean Distance",
     hang  = -1,
     cex   = 0.8)

# Elbow / scree plot
wss <- sapply(1:10, function(k) {
  ct <- cutree(hc, k)
  sum(sapply(unique(ct), function(cl) {
    rows <- ts_matrix_complete[ct == cl, , drop = FALSE]
    sum(apply(rows, 2, var, na.rm = TRUE)) * (sum(ct == cl) - 1)
  }))
})
plot(1:10, wss, type = "b", xlab = "Number of clusters", ylab = "WSS")

# --- 6. Cut tree into k clusters (choose k) ---
k <- 5  # adjust as needed
clusters <- cutree(hc, k = k)

# Add cluster assignments back to original data
cluster_df <- tibble(
  pheno   = names(clusters),
  cluster = clusters
)

data_with_clusters <- df %>%
  left_join(cluster_df, by = "pheno")

# --- 7. Visualise time series coloured by cluster ---
data_with_clusters %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = logAge_days, y = !!sym(clust_var),
             group = pheno, colour = cluster)) +
  geom_line(alpha = 0.7) +
  facet_wrap(pheno_cat~cluster) +
  labs(
    title  = "Pheno Time Series by Cluster",
    x      = "log(Age in days)",
    y      = "Centile M minus F z-score",
    colour = "Cluster"
  ) +
  theme_minimal()

# --- 8. Plot average trajectory per cluster ---
data_with_clusters %>%
  mutate(cluster = factor(cluster)) %>%
  group_by(cluster, logAge_days) %>%
  summarise(
    mean_trajectory = mean(!!sym(clust_var), na.rm = TRUE),
    se_centile   = sd(!!sym(clust_var), na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = logAge_days, y = mean_trajectory, colour = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = mean_trajectory - se_centile,
                  ymax = mean_trajectory + se_centile),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~cluster, scales='free_y') +
  labs(
    title    = "Average Trajectory per Cluster",
    subtitle = "Ribbon = ±1 SE across pheno levels",
    x        = "log(Age in days)",
    colour   = "Cluster",
    fill     = "Cluster"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
