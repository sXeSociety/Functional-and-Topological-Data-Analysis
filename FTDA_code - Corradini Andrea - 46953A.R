# Load the libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(fda)
library(zoo)
library(lubridate)
library(cluster)
library(factoextra)
library(patchwork)
library(stringr)
library(depthTools)
library(tibble)
library(fda.usc)
library(TDA)
library(rgl)

# Define base data directory (relative to project root)
data_dir <- file.path("Replace-BG_Dataset", "Data Tables")
# Load datasets
cgm <- read_delim(file.path(data_dir, "HDeviceCGM.txt"), delim = "|", show_col_types = FALSE)
insulin <- read_delim(file.path(data_dir, "HDeviceBolus.txt"), delim = "|", show_col_types = FALSE)
HbA1c <- read_delim(file.path(data_dir, "HLocalHbA1c.txt"), delim = "|", show_col_types = FALSE)
Quests2 <- read_delim(file.path(data_dir, "HQuestHypoFear.txt"), delim = "|", show_col_types = FALSE)

# Convert device timestamps to full POSIX datetime
cgm <- cgm %>%
  mutate(
    seconds_from_enroll = DeviceDtTmDaysFromEnroll * 86400,
    seconds_from_midnight = as.numeric(DeviceTm),
    datetime = as.POSIXct(seconds_from_enroll + seconds_from_midnight, origin = "2017-01-01", tz = "UTC")
  ) %>%
  filter(!is.na(GlucoseValue))  # Remove missing values

# Keep only days with at least 260 CGM entries
cgm <- cgm %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(PtID, date) %>%
  filter(n() >= 260) %>%
  ungroup()

# Identify patients with at least 7 valid days
valid_days <- cgm %>%
  group_by(PtID, date) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  filter(n_obs >= 260) %>%
  group_by(PtID) %>%
  summarise(n_days = n(), .groups = "drop") %>%
  filter(n_days >= 7)

# Filter for selected patients
selected_ids <- valid_days$PtID
cgm_sub <- cgm %>% filter(PtID %in% selected_ids)

# Convert time to minutes since midnight
cgm_sub <- cgm_sub %>%
  mutate(time_min = as.numeric(DeviceTm) / 60)

# Define 288-point time grid (5-minute resolution over 24h)
grid_minutes <- seq(0, 1440, length.out = 288)

# Interpolate glucose values on the fixed time grid
#curves_df <- cgm_sub %>%
  #group_by(PtID, date) %>%
  #arrange(time_min, .by_group = TRUE) %>%
  #distinct(time_min, .keep_all = TRUE) %>%
  #summarise(glucose_interp = list(approx(time_min, GlucoseValue, xout = grid_minutes)$y), .groups = "drop") %>%
  #mutate(curve_id = paste(PtID, date, sep = "_"))
#saveRDS(curves_df, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/curves_df.rds")
curves_df <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/curves_df.rds")

# Create raw glucose matrix
glucose_matrix <- do.call(cbind, curves_df$glucose_interp)

# Keep only curves with ≤ 10% missing values
na_counts <- apply(glucose_matrix, 2, function(col) sum(is.na(col)))
keep_columns <- which(na_counts <= 28)
glucose_matrix_partial <- glucose_matrix[, keep_columns]
curves_df_partial <- curves_df[keep_columns, ]

# --- FUNCTIONS ---

# Function to impute missing values using linear interpolation
impute_missing_values <- function(mat, grid) {
  apply(mat, 2, function(col) na.approx(col, x = grid, rule = 2))
}

# Function to smooth glucose curves
smooth_curves <- function(y_matrix, grid, nbasis = 25, rangeval = range(grid)) {
  basis <- create.bspline.basis(rangeval = rangeval, nbasis = nbasis)
  smooth.basis(argvals = grid, y = y_matrix, fdParobj = basis)$fd
}

# Function to plot 9 smoothed curves with time of day labels
plot_sample_curves <- function(fd_obj, metadata, time_pos, time_labs) {
  par(mfrow = c(3, 3))
  for (i in 1:9) {
    plot(fd_obj[i],
         main = paste0("Curve ", i, "\nPtID: ", metadata$PtID[i], "\nDate: ", metadata$date[i]),
         ylab = "Glucose (mg/dL)", xlab = "Time of day",
         col = "steelblue", lwd = 2, xaxt = "n")
    axis(1, at = time_pos, labels = time_labs)
  }
}

# --- PROCESSING ---

# Normalize time from [0, 1440] to [0, 1]
normalized_time <- grid_minutes / 1440
# Impute missing values
#glucose_matrix_imputed <- impute_missing_values(glucose_matrix_partial, normalized_time)
#saveRDS(glucose_matrix_imputed, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/glucose_matrix_imputed.rds")
glucose_matrix_imputed <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/glucose_matrix_imputed.rds")
glucose_matrix_imputed <- matrix(glucose_matrix_imputed, nrow = 288)

# Smooth curves on normalized time
#smoothed_fd <- smooth_curves(glucose_matrix_imputed, grid = normalized_time)
#saveRDS(smoothed_fd, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/smoothed_fd.rds")
smoothed_fd <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/smoothed_fd.rds")

# Smooth again on real-time scale [0, 86400]
real_time <- seq(0, 86400, length.out = 288)
smoothed_fd_real <- smooth_curves(glucose_matrix_imputed, grid = real_time, nbasis = 21, rangeval = c(0, 86400))

# Select 3 patients × 3 days for visualization
unique_ids <- unique(curves_df$PtID)
selected_ids_plot <- unique_ids[1:3]
selected_indexes <- unlist(lapply(selected_ids_plot, function(id) which(curves_df_partial$PtID == id)[1:3]))

# Time labels for the X-axis
time_labels <- c("00:00", "06:00", "12:00", "18:00", "24:00")
time_positions <- c(0, 21600, 43200, 64800, 86400)

# Save plot
png(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/1. Curves.png",
    width = 10, height = 8, units = "in", res = 300)
# Plot the selected 9 curves
plot_sample_curves(smoothed_fd_real[selected_indexes], curves_df_partial[selected_indexes, ], time_positions, time_labels)
dev.off()

# --- FUNCTIONAL PCA ---

# Perform Functional Principal Component Analysis on the smoothed functional data
# We retain the first 10 harmonics (adjustable if needed)
fpca_res <- pca.fd(fdobj = smoothed_fd, nharm = 10, centerfns = TRUE)

# Build dataframe with percentage of variance explained
fpca_df <- data.frame(
  Component = factor(1:length(fpca_res$varprop)),
  Variance = fpca_res$varprop * 100  # Convert to %
)

# Plot the scree plot with bars and line
ggplot(fpca_df, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.6) +
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  geom_point(color = "black", size = 2.5) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)),
            vjust = -0.5, hjust = -0.3, size = 3) +
  labs(title = "Scree Plot – Functional PCA",
       x = "Functional Principal Components",
       y = "Variance Explained (%)") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  ylim(0, max(fpca_df$Variance) + 12)
# Save plot as PNG
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/2. Scree Plot.png", 
       width = 10, height = 6, dpi = 300)

# Extract FPCA scores and add metadata
scores_df <- as.data.frame(fpca_res$scores)
colnames(scores_df) <- paste0("Comp.", seq_len(ncol(scores_df)))
scores_df <- scores_df %>%
  mutate(curve_id = curves_df_partial$curve_id,
         PtID = curves_df_partial$PtID,
         date = as.Date(curves_df_partial$date))

# Compute enrollment date per patient (earliest CGM reading)
enroll_dates <- cgm %>%
  group_by(PtID) %>%
  summarise(enroll_date = min(as.Date(datetime)), .groups = "drop")

# Merge and compute day difference between curve and HbA1c test
scores_df <- scores_df %>%
  left_join(enroll_dates, by = "PtID") %>%
  mutate(days_after_enroll = as.integer(date - enroll_date))

# Clean HbA1c test results (exclude "not done")
HbA1c_clean <- HbA1c %>%
  filter(is.na(HbA1cNotDone) | HbA1cNotDone != 1) %>%
  dplyr::select(PtID, HbA1cTestDtDaysAfterEnroll, HbA1cTestRes)

# Join FPCA scores with HbA1c (many-to-many match by PtID)
scores_merged <- scores_df %>%
  left_join(HbA1c_clean, by = "PtID", relationship = "many-to-many") %>%
  mutate(day_gap = days_after_enroll - HbA1cTestDtDaysAfterEnroll) %>%
  filter(day_gap >= 0)

# Keep the HbA1c test closest in time (but before the curve)
scores_final <- scores_merged %>%
  group_by(curve_id) %>%
  slice_min(day_gap, with_ties = FALSE) %>%
  ungroup()


# --- FPCA SCATTERPLOT ---

# Sample 500 curves for plotting
set.seed(123)
scores_sample <- scores_final %>% sample_n(500)

# FPCA biplot: PC1 vs PC2 colored by HbA1c
ggplot(scores_sample, aes(x = Comp.1, y = Comp.2, color = HbA1cTestRes)) +
  geom_point(size = 2.5, alpha = 1) +
  scale_color_gradientn(
    colours = c("#08306B", "#2171B5", "#6BAED6", "#FD8D3C", "#BD0026"),
    name = "HbA1c"
  ) +
  labs(
    title = "FPCA Scores – PC1 vs PC2 colored by HbA1c",
    x = "Functional Principal Component 1",
    y = "Functional Principal Component 2"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right"
  )

# Save plot as PNG
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/3. Scatterplot.png", 
       width = 10, height = 6, dpi = 300)

# Print selected columns for visual check
scores_sample %>%
  dplyr::select(PtID, date, Comp.1, Comp.2, HbA1cTestRes) %>%
  print(n = 500)

# --- HARMONICS PLOTTING ---

# Evaluate mean function and harmonics on normalized grid
mean_eval <- eval.fd(normalized_time, fpca_res$meanfd)
harmonics_eval <- eval.fd(normalized_time, fpca_res$harmonics)

# Reuse already defined real-time grid and labels
# (time_real and time_labels_custom already defined above)

# Save harmonics plot
png(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/4. Harmonics.png",
    width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(3, 2))  # Layout: 3 rows × 2 columns

# Plot mean glucose curve
plot(real_time, mean_eval, type = "l", col = "black", lwd = 2,
     xaxt = "n", main = "Mean glucose curve",
     xlab = "Time of day", ylab = "Glucose (mg/dL)")
axis(1, at = time_positions, labels = time_labels)

# Plot first 5 harmonics
for (k in 1:5) {
  plot(real_time, harmonics_eval[, k], type = "l", col = "darkred", lwd = 2,
       xaxt = "n", main = paste("Harmonic", k),
       xlab = "Time of day", ylab = expression(phi(t)))
  abline(h = 0, lty = 2)
  axis(1, at = time_positions, labels = time_labels)
}

dev.off()

# --- CLUSTERING ON FPCA SCORES ---

# Use the first 3 FPCA components and scale them for clustering
fpca_scores <- scores_final %>%
  dplyr::select(Comp.1:Comp.3) %>%
  as.matrix() %>%
  scale()

# Sample 500 observations for elbow and silhouette analysis
set.seed(123)
fpca_sample <- fpca_scores[sample(1:nrow(fpca_scores), 500), ]

# Elbow method to choose optimal number of clusters
fviz_nbclust(fpca_sample, kmeans, method = "wss") +
  labs(title = "WSS Plot (Elbow Method)",
       x = "Number of Clusters",
       y = "Within-Cluster Sum of Squares") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
# Saving the plot as an image file
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/5. WSS.png", 
       plot = last_plot(), width = 10, height = 6, dpi = 300)

# Silhouette method to assess optimal k
fviz_nbclust(fpca_sample, kmeans, method = "silhouette") +
  labs(title = "Silhouette Plot",
       x = "Number of Clusters",
       y = "Average Silhouette Width") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
# Saving the plot as an image file
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/6. Silhouette.png", 
       plot = last_plot(), width = 10, height = 6, dpi = 300)

# Apply k-means clustering with optimal k
k <- 4
set.seed(123)
kmeans_result <- kmeans(fpca_scores, centers = k, nstart = 50)

# Add cluster assignments to FPCA scores
scores_clustered <- scores_final %>%
  mutate(cluster = factor(kmeans_result$cluster))
# Silhouette analysis: use same sample
sample_idx <- sample(1:nrow(fpca_scores), size = 500)  
# Ensure it's a matrix with correct dimensions
scores_sample <- fpca_scores[sample_idx, , drop = FALSE]
cluster_sample <- kmeans_result$cluster[sample_idx]
# Compute distance
diss_sample <- dist(scores_sample)
# Compute silhouette object
silhouette_obj <- cluster::silhouette(cluster_sample, diss_sample)
# Define a color palette
custom_colors <- c("#08306B", "#6BAED6", "#FD8D3C", "#BD0026")

# Silhouette plot
fviz_silhouette(silhouette_obj) +  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Silhouette Plot – K-Means (Sample of 500)",
       y = "Silhouette Width") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
# Saving the plot as an image file
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/7. Silhouette2.png", 
       plot = last_plot(), width = 10, height = 6, dpi = 300)

# FPCA biplot with clusters
ggplot(scores_clustered, aes(x = Comp.1, y = Comp.2, color = cluster)) +
  geom_point(alpha = 0.55, size = 2) +
  scale_color_manual(values = custom_colors) +
  labs(title = "K-Means Clustering on FPCA Scores",
       x = "FPCA Component 1",
       y = "FPCA Component 2",
       color = "Cluster") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom")
# Saving the plot as an image file
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/8. Kmeans.png", 
       plot = last_plot(), width = 10, height = 6, dpi = 300)

# Cluster-wise mean HbA1c
scores_clustered %>%
  group_by(cluster) %>%
  summarise(avg_HbA1c = mean(HbA1cTestRes, na.rm = TRUE),
            count = n())

# --- CLUSTER-WISE GLUCOSE PROFILES ---

# Evaluate smoothed curves (from FPCA) on normalized time grid
glucose_matrix_clustered <- eval.fd(normalized_time, smoothed_fd)
colnames(glucose_matrix_clustered) <- curves_df_partial$curve_id

# Keep only curves assigned to clusters
valid_ids <- intersect(colnames(glucose_matrix_clustered), scores_clustered$curve_id)
glucose_matrix_clustered <- glucose_matrix_clustered[, valid_ids]

# Build long dataframe for plotting
glucose_df <- as.data.frame(glucose_matrix_clustered)
glucose_df$time <- normalized_time

glucose_long <- pivot_longer(glucose_df, cols = -time,
                             names_to = "curve_id", values_to = "glucose") %>%
  left_join(scores_clustered %>% dplyr::select(curve_id, cluster), by = "curve_id")

# Compute mean and 95% CI for each cluster
glucose_summary <- glucose_long %>%
  group_by(cluster, time) %>%
  summarise(
    mean = mean(glucose, na.rm = TRUE),
    sd = sd(glucose, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_low = mean - qt(0.975, df = n - 1) * se,
    ci_high = mean + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

# Plot cluster-wise average glucose curves
ggplot(glucose_summary, aes(x = time * 24, y = mean, color = cluster, fill = cluster)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(0, 24, 6),
                     labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Cluster-wise Mean Glucose Curves",
       x = "Time of Day", y = "Glucose (mg/dL)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
# Saving the plot as an image file
ggsave(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/9. Curves-Cluster.png", 
       plot = last_plot(), width = 10, height = 6, dpi = 300)

# --- FUNCTIONAL BOXPLOT WITH MBD AND BD2 ---

# Save side-by-side functional boxplots as high-resolution PNG
png(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/10. FunctionalBoxplots.png",
    width = 1200, height = 600, res = 150)

# Set 1 row, 2 columns layout 
par(mfrow = c(1, 2),
    mar = c(5, 5, 4, 2),
    oma = c(0, 0, 0, 0),
    mgp = c(1, 0.2, 0),
    cex.main = 1,
    cex.lab = 0.7,
    cex.axis = 0.55,
    las = 0.8)

# Functional Boxplot (MBD)
boxplot(smoothed_fd, method = "MBD",
        main = "",
        xlab = "Time of Day", ylab = "Glucose (mg/dL)",
        col = "lightsteelblue1",
        frame.plot = TRUE,
        xaxt = "n")
axis(1, at = seq(0, 1, length.out = 5),
     labels = c("00:00", "06:00", "12:00", "18:00", "24:00"))
box(lwd = 0.5)
# Add title closer to plot
title(main = "Functional Boxplot (MBD)", line = 1)

# Functional Boxplot (BD2)
boxplot(smoothed_fd, method = "BD2",
        main = "",
        xlab = "Time of Day", ylab = "Glucose (mg/dL)",
        col = "lightsteelblue1",
        frame.plot = TRUE,
        xaxt = "n")
axis(1, at = seq(0, 1, length.out = 5),
     labels = c("00:00", "06:00", "12:00", "18:00", "24:00"))
box(lwd = 0.5)
# Add title closer to plot
title(main = "Functional Boxplot (BD2)", line = 1)

dev.off()

# --- DEPTH MEASURES (MBD) ---

# Evaluate smoothed functional data on [0,1] grid 
fd_matrix <- t(eval.fd(normalized_time, smoothed_fd))  # rows = curves
n_curves <- nrow(fd_matrix)

# Compute MBD on all curves
#mbd_all <- MBD(fd_matrix)$MBD
#saveRDS(mbd_all, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/mbd_all.rds")
mbd_all <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/mbd_all.rds")
mbd_vec <- (mbd_all - min(mbd_all)) / (max(mbd_all) - min(mbd_all))  # normalize

# Save histogram as PNG
png(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/11. Histogram_mbd_all.png",
    width = 1000, height = 600, res = 150)

# Set consistent graphical parameters
par(mar = c(5, 5, 4, 2),            # margins: bottom, left, top, right
    cex.main = 1,                   # title size
    cex.lab = 0.7,                  # axis labels
    cex.axis = 0.55,                # axis ticks
    mgp = c(1, 0.2, 0),             # label positioning
    las = 0.8)                      # horizontal y-axis

# Draw histogram
hist(mbd_vec,
     breaks = 30,
     col = "steelblue",
     border = "black",             
     main = "",
     xlab = "MBD Value",
     ylab = "Frequency",
     xaxt = "n")                   
axis(1, lwd = 0.5, cex.axis = 0.55)
box(lwd = 0.5) 
title(main = "Histogram – Modified Band Depth (MBD)", line = 1)

# Close device
dev.off()

# --- WILCOXON / KRUSKAL-WALLIS TEST: MBD VS MULTI HbA1c GROUPS ---

# Use curve IDs corresponding to partial matrix
curve_ids_all <- curves_df_partial$curve_id
if (length(curve_ids_all) != n_curves) {
  stop("Number of curve IDs does not match number of curves in fd_matrix")
}
ptid_all <- sapply(strsplit(curve_ids_all, "_"), `[`, 1)

# Map HbA1c values to PtIDs
hb_values <- HbA1c %>%
  group_by(PtID) %>%
  summarise(HbA1c = mean(HbA1cTestRes, na.rm = TRUE)) %>%
  deframe()

hb_all <- hb_values[ptid_all]

# Define HbA1c group categories
group_labels_all <- cut(
  hb_all,
  breaks = c(-Inf, 6.5, 7.5, 8.5, Inf),
  labels = c("≤6.5", "(6.5,7.5]", "(7.5,8.5]", ">8.5"),
  right = TRUE
)

# Filter valid data
valid_idx <- which(!is.na(group_labels_all) & !is.na(mbd_vec))
group_valid <- group_labels_all[valid_idx]
depth_valid <- mbd_vec[valid_idx]

# Kruskal-Wallis test
kruskal_result <- kruskal.test(depth_valid ~ group_valid)
cat("Kruskal-Wallis Test (MBD vs HbA1c groups):\n")
print(kruskal_result)

# Pairwise Wilcoxon tests
pairwise_result <- pairwise.wilcox.test(depth_valid, group_valid, p.adjust.method = "BH")
cat("\nPairwise Wilcoxon Tests (BH corrected):\n")
print(pairwise_result)

# Create the data frame for plotting
df_plot <- data.frame(
  MBD = depth_valid,
  HbA1cGroup = group_valid
)

ggplot(df_plot, aes(x = HbA1cGroup, y = MBD, fill = HbA1cGroup)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  labs(title = "MBD by HbA1c Group",
       x = "HbA1c Group",
       y = "Normalized MBD") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 9),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("C:/Users/andre/OneDrive/Desktop/FTDA/Plots/12. MBD_violin_boxplot.png",
       width = 10, height = 6, dpi = 300)

# ----- REGRESSION -----

# Keep only complete cases for Comp.1 and HbA1c
lm_data <- scores_clustered %>%
  filter(!is.na(Comp.1), !is.na(HbA1cTestRes))

# Fit linear model: PC1 as response, HbA1c as predictor
lm_fpca <- lm(Comp.1 ~ HbA1cTestRes, data = lm_data)

# Summary of the linear regression
summary(lm_fpca)

# Compute number of boluses per day 
boli_per_day <- insulin %>%
  mutate(date = as.Date(DeviceDtTmDaysFromEnroll, origin = "2017-01-01")) %>%
  group_by(PtID, date) %>%
  summarise(n_boli = sum(!is.na(Normal) & Normal > 0), .groups = "drop")

# Compute average fear score per PtID
fear_avg <- Quests2 %>%
  dplyr::select(PtID, starts_with("Worry"), starts_with("Avoid"), starts_with("Keep"), 
                starts_with("Test"), starts_with("Eat"), starts_with("Red"), 
                starts_with("Carry"), starts_with("Ck")) %>%
  rowwise() %>%
  mutate(mean_fear = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(PtID) %>%
  summarise(mean_fear = mean(mean_fear, na.rm = TRUE), .groups = "drop")

# Merge everything into a complete dataset
lm_data <- scores_clustered %>%
  left_join(boli_per_day, by = c("PtID", "date")) %>%
  left_join(fear_avg, by = "PtID") %>%
  filter(!is.na(HbA1cTestRes), !is.na(n_boli), !is.na(mean_fear), 
         !is.na(Comp.1), !is.na(Comp.2), !is.na(Comp.3))

# Standardize covariates
lm_data <- lm_data %>%
  mutate(
    HbA1c_z = scale(HbA1cTestRes)[,1],
    n_boli_z = scale(n_boli)[,1],
    mean_fear_z = scale(mean_fear)[,1]
  )

# --- Fit regressions for PC1, PC2, PC3 ---
lm_PC1 <- lm(Comp.1 ~ HbA1c_z + n_boli_z + mean_fear_z, data = lm_data)
lm_PC2 <- lm(Comp.2 ~ HbA1c_z + n_boli_z + mean_fear_z, data = lm_data)
lm_PC3 <- lm(Comp.3 ~ HbA1c_z + n_boli_z + mean_fear_z, data = lm_data)

# Show results
summary(lm_PC1)
summary(lm_PC2)
summary(lm_PC3)

# --- FUNCTIONAL REGRESSION ---

# Prepare metadata for regression
# Merge HbA1c, number of boluses, and fear scores with curve metadata
metadata_fd <- curves_df_partial %>%
  dplyr::select(curve_id, PtID, date) %>%
  left_join(HbA1c_clean, by = "PtID", relationship = "many-to-many") %>%
  mutate(
    enroll_date = enroll_dates$enroll_date[match(PtID, enroll_dates$PtID)],
    days_after_enroll = as.integer(date - enroll_date),
    day_gap = days_after_enroll - HbA1cTestDtDaysAfterEnroll
  ) %>%
  filter(day_gap >= 0) %>%
  group_by(curve_id) %>%
  slice_min(day_gap, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(boli_per_day, by = c("PtID", "date")) %>%
  left_join(fear_avg, by = "PtID") %>%
  filter(!is.na(HbA1cTestRes), !is.na(n_boli), !is.na(mean_fear))

# Standardized variables (HbA1c, n_boli, fear)
metadata_fd <- metadata_fd %>%
  mutate(
    HbA1c_z = scale(HbA1cTestRes),
    n_boli_z = scale(n_boli),
    fear_z = scale(mean_fear)
  )

# Match curves used in FPCA and regression
idx_keep <- which(curves_df_partial$curve_id %in% metadata_fd$curve_id)
fd_regression <- smoothed_fd[idx_keep]
metadata_fd <- metadata_fd[match(curves_df_partial$curve_id[idx_keep], metadata_fd$curve_id), ]
stopifnot(all.equal(curves_df_partial$curve_id[idx_keep], metadata_fd$curve_id))

# Define normalized time grid for evaluation
normalized_time <- seq(0, 1, length.out = 288)

# --- MODEL 1: Scalar-on-function REGRESSION – Only glucose curve ---

# Response: standardized HbA1c
HbA1c_z <- as.vector(metadata_fd$HbA1c_z)
# Predictor: functional glucose curve only
xfdlist1 <- list(fd_regression)
betalist1 <- list(fdPar(fd_regression$basis, lambda = 1e-2))
#freg1 <- fRegress(y = HbA1c_z, xfdlist = xfdlist1, betalist = betalist1)
#saveRDS(freg1, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/freg1.rds")
freg1 <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/freg1.rds")

# Evaluate beta₁(t)
beta1_fd <- freg1$betaestlist[[1]]$fd
beta1_eval <- eval.fd(normalized_time, beta1_fd)

# Plot β₁(t)
beta1_df <- data.frame(
  time = normalized_time * 24,
  beta = beta1_eval)

# Plot the functional coefficient with ggplot
ggplot(beta1_df, aes(x = time, y = beta)) +
  geom_line(color = "steelblue", linewidth = 1.2) +  # Smooth line
  geom_hline(yintercept = 0, linetype = "dashed") +  # Horizontal reference
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 24),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00")
  ) +
  labs(
    title = expression("Functional Coefficient " * beta[1](t) * " (HbA1c ~ curve only)"),
    x = "Time of Day",
    y = expression(beta(t))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("C:/Users/andre/OneDrive/Desktop/FTDA/Plots/13. Beta1_only_curve.png", 
       width = 10, height = 6, dpi = 300)

# --- MODEL 2: Scalar-on-function REGRESSION – Curve + n_boli + fear ---

# Create constant basis for scalar predictors
const_basis <- create.constant.basis(c(0, 1))

# Convert scalar predictors to functional data objects
n_boli_fd <- fd(matrix(metadata_fd$n_boli_z, nrow = 1), const_basis)
fear_fd   <- fd(matrix(metadata_fd$fear_z, nrow = 1), const_basis)

# Define functional predictor list: [glucose_curve_fd, n_boli_fd, fear_fd]
xfdlist2 <- list(fd_regression, n_boli_fd, fear_fd)

# Define beta list as before (1 smooth, 2 constant)
betalist2 <- list(
  fdPar(fd_regression$basis, lambda = 1e-2),  # for glucose curve
  fdPar(const_basis),                         # for n_boli
  fdPar(const_basis)                          # for fear
)

# Now fit the model
#freg2 <- fRegress(y = HbA1c_z, xfdlist = xfdlist2, betalist = betalist2)
#saveRDS(freg2, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/freg2.rds")
freg2 <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/freg2.rds")

# Evaluate functional coefficient β₁(t)
beta2_fd <- freg2$betaestlist[[1]]$fd
beta2_eval <- eval.fd(normalized_time, beta2_fd)

# Plot β₁(t) from the model with covariates
beta2_df <- data.frame(
  time = normalized_time * 24,
  beta = beta2_eval
)

# Plot the functional coefficient with ggplot
ggplot(beta2_df, aes(x = time, y = beta)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +  # Smooth line
  geom_hline(yintercept = 0, linetype = "dashed") +  # Horizontal reference
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 24),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00")
  ) +
  labs(
    title = expression("Functional Coefficient " * beta[1](t) * " (HbA1c ~ curve + covariates)"),
    x = "Time of Day",
    y = expression(beta(t))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("C:/Users/andre/OneDrive/Desktop/FTDA/Plots/14. Beta1_with_covariates.png", 
       width = 10, height = 6, dpi = 300)

# --- MODEL 3: Function-on-scalar REGRESSION – Curve ~ n_boli + fear ---

# Define design matrix for scalar predictors
X_mat <- model.matrix(~ n_boli_z + fear_z, data = metadata_fd)
# Define the functional response
Y_fd <- fd_regression  # glucose curves as response

# Set up the basis and penalization
beta_basis_list <- list(
  fdPar(fd_regression$basis, lambda = 1e-2),  # intercept
  fdPar(fd_regression$basis, lambda = 1e-2),  # n_boli_z
  fdPar(fd_regression$basis, lambda = 1e-2)   # fear_z
  )

# Convert design matrix (excluding intercept column) to a list of numeric vectors
xfdlist_fos <- list(
  as.numeric(X_mat[, 1]),  # Intercept
  as.numeric(X_mat[, 2]),  # n_boli_z
  as.numeric(X_mat[, 3])   # fear_z
  )

# Fit function-on-scalar model
#fos_reg <- fRegress(y = Y_fd, xfdlist = xfdlist_fos, betalist = beta_basis_list)
#saveRDS(fos_reg, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/fos_reg.rds")
fos_reg <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/fos_reg.rds")

# Evaluate β₂(t) and β₃(t)
beta_nboli <- eval.fd(normalized_time, fos_reg$betaestlist[[2]]$fd)
beta_fear <- eval.fd(normalized_time, fos_reg$betaestlist[[3]]$fd)

# Plot both β₂(t) and β₃(t)
beta_fos_df <- data.frame(
  time = normalized_time * 24,
  n_boli = beta_nboli,
  fear = beta_fear)

# Plot the two functional coefficients with ggplot
ggplot(beta_fos_df, aes(x = time)) +
  geom_line(aes(y = n_boli, color = "n_boli_z"), linewidth = 1.2) +  # Coef for n_boli
  geom_line(aes(y = fear, color = "fear_z"), linewidth = 1.2) +      # Coef for fear
  geom_hline(yintercept = 0, linetype = "dashed") +                  # Horizontal reference
  scale_color_manual(values = c("n_boli_z" = "#2171B5", "fear_z" = "#D7301F")) +
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 24),
    labels = c("00:00", "06:00", "12:00", "18:00", "24:00")
  ) +
  labs(
    title = "Function-on-scalar regression: effects over time",
    x = "Time of Day",
    y = expression(beta(t)),
    color = "Covariate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("C:/Users/andre/OneDrive/Desktop/FTDA/Plots/15. Function_on_scalar.png", 
       width = 10, height = 6, dpi = 300)

# FANOVA

# Ensure the clusters are available for each curve
# fd_regression: functional data object for selected curves
# metadata_fd$cluster: factor containing the cluster label for each curve

# Merge cluster info into metadata if not already done
metadata_fd <- metadata_fd %>%
  left_join(scores_clustered %>% dplyr::select(curve_id, cluster), by = "curve_id") %>%
  filter(!is.na(cluster))

# Evaluate smoothed curves on normalized grid 
fd_matrix <- eval.fd(normalized_time, fd_regression)

# Create fdata object: rows = curves, columns = time points
fdata_obj <- fdata(t(fd_matrix), argvals = normalized_time)

# --- FANOVA by Cluster ---

# Grouping factor by cluster
cluster_factor <- as.factor(metadata_fd$cluster)
# Run FANOVA
#fanova_cluster <- fanova.onefactor(fdata_obj, group = cluster_factor)
#saveRDS(fanova_cluster, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/fanova_cluster.rds")
fanova_cluster <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/fanova_cluster.rds")

# Extract F-statistics and compute threshold
f_stat_cluster <- fanova_cluster$wm
time_hours <- seq(0, 24, length.out = length(f_stat_cluster))
threshold_cluster <- quantile(f_stat_cluster, 0.95)

# Prepare data for plotting
fanova_df_cluster <- data.frame(
  time = time_hours,
  F_stat = f_stat_cluster)

# Plot F-statistic curve (clusters)
ggplot(fanova_df_cluster, aes(x = time, y = F_stat)) +
  geom_line(color = "steelblue", linewidth = 1.1) +
  geom_hline(yintercept = threshold_cluster, linetype = "dashed", color = "darkred") +
  scale_x_continuous(breaks = seq(0, 24, by = 6),
                     labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  labs(title = "FANOVA – F-statistic by Cluster",
       x = "Time of Day", y = "F-statistic") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

ggsave("C:/Users/andre/OneDrive/Desktop/FTDA/Plots/16. Fanova.png", 
       width = 10, height = 6, dpi = 300)

# --- FANOVA by HbA1c multi-category ---

# Create multi-category factor for HbA1c
hba1c_factor_multi2 <- cut(
  metadata_fd$HbA1cTestRes,
  breaks = c(-Inf, 6.5, 7.5, 8.5, Inf),
  labels = c("≤6.5", "(6.5,7.5]", "(7.5,8.5]", ">8.5")
)

# Run FANOVA
#fanova_hba1c_multi2 <- fanova.onefactor(fdata_obj, group = hba1c_factor_multi2)
#saveRDS(fanova_hba1c_multi2, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/fanova_hba1c_multi2.rds")
fanova_hba1c_multi2 <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/fanova_hba1c_multi2.rds")

# Extract F-statistics and threshold
f_stat_hba1c_multi2 <- fanova_hba1c_multi2$wm
threshold_hba1c_multi2 <- quantile(f_stat_hba1c_multi2, 0.95)

# Prepare data for plotting
fanova_df_hba1c_multi2 <- data.frame(
  time = time_hours,
  F_stat = f_stat_hba1c_multi2
)

# Plot F-statistic curve (HbA1c multi-group)
ggplot(fanova_df_hba1c_multi2, aes(x = time, y = F_stat)) +
  geom_line(color = "steelblue", linewidth = 1.1) +
  geom_hline(yintercept = threshold_hba1c_multi2, linetype = "dashed", color = "darkred") +
  scale_x_continuous(breaks = seq(0, 24, by = 6),
                     labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  labs(title = "FANOVA – Glucose Curves by HbA1c Multi-Group",
       x = "Time of Day", y = "F-statistic") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

ggsave("C:/Users/andre/OneDrive/Desktop/FTDA/Plots/17. Fanova_Multi.png", 
       width = 10, height = 6, dpi = 300)

# --- CROSS-VALIDATION (10-fold) ON A SUBSAMPLE ---

set.seed(42)
# Select 1000 random curves
sample_idx <- sample(seq_along(HbA1c_z), size = 1000)
# Subsample response and functional predictor
HbA1c_z_sub <- HbA1c_z[sample_idx]
fd_regression_sub <- fd_regression[sample_idx]

# --- CV for Model 1: curve only ---

xfdlist1_sub <- list(fd_regression_sub)
betalist1_sub <- list(fdPar(fd_regression$basis, lambda = 1e-2))
#cv1_sub <- fRegress.CV(HbA1c_z_sub, xfdlist1_sub, betalist1_sub, nfold = 10)
#saveRDS(cv1_sub, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/cv1_sub.rds")
cv1_sub <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/cv1_sub.rds")

# Compute R² from cross-validated SSE
cv1_r2 <- 1 - cv1_sub$SSE.CV / sum((HbA1c_z_sub - mean(HbA1c_z_sub))^2)
cat("Model 1 – CV R² (subset):", round(cv1_r2, 4), "\n")

# --- CV for Model 2: curve + n_boli + fear ---

# Subsample scalar covariates
n_boli_z_sub <- metadata_fd$n_boli_z[sample_idx]
fear_z_sub   <- metadata_fd$fear_z[sample_idx]

# Recreate functional data objects for scalar predictors
n_boli_fd_sub <- fd(matrix(n_boli_z_sub, nrow = 1), const_basis)
fear_fd_sub   <- fd(matrix(fear_z_sub, nrow = 1), const_basis)

xfdlist2_sub <- list(fd_regression_sub, n_boli_fd_sub, fear_fd_sub)
betalist2_sub <- list(
  fdPar(fd_regression$basis, lambda = 1e-2),
  fdPar(const_basis),
  fdPar(const_basis))

#cv2_sub <- fRegress.CV(HbA1c_z_sub, xfdlist2_sub, betalist2_sub, nfold = 10)
#saveRDS(cv2_sub, file = "C:/Users/andre/OneDrive/Desktop/FTDA/Objects/cv2_sub.rds")
cv2_sub <- readRDS("C:/Users/andre/OneDrive/Desktop/FTDA/Objects/cv2_sub.rds")

cv2_r2 <- 1 - cv2_sub$SSE.CV / sum((HbA1c_z_sub - mean(HbA1c_z_sub))^2)
cat("Model 2 – CV R² (subset):", round(cv2_r2, 4), "\n")

# --- TOPOLOGICAL DIFFERENCES – HbA1c 4 CATEGORIES ---

# Identify indices for each group
idx_1 <- which(metadata_fd$HbA1cTestRes <= 6.5)
idx_2 <- which(metadata_fd$HbA1cTestRes > 6.5 & metadata_fd$HbA1cTestRes <= 7.5)
idx_3 <- which(metadata_fd$HbA1cTestRes > 7.5 & metadata_fd$HbA1cTestRes <= 8.5)
idx_4 <- which(metadata_fd$HbA1cTestRes > 8.5)

# Select up to 100 curves per group
curve_1 <- t(glucose_matrix_imputed[, idx_1[1:min(100, length(idx_1))]])
curve_2 <- t(glucose_matrix_imputed[, idx_2[1:min(100, length(idx_2))]])
curve_3 <- t(glucose_matrix_imputed[, idx_3[1:min(100, length(idx_3))]])
curve_4 <- t(glucose_matrix_imputed[, idx_4[1:min(100, length(idx_4))]])

# Compute persistence diagrams
Diag_1 <- alphaShapeDiag(curve_1, library = "GUDHI")
Diag_2 <- alphaShapeDiag(curve_2, library = "GUDHI")
Diag_3 <- alphaShapeDiag(curve_3, library = "GUDHI")
Diag_4 <- alphaShapeDiag(curve_4, library = "GUDHI")

png(filename = "C:/Users/andre/OneDrive/Desktop/FTDA/Plots/18. Persistence_Diagrams_Multi.png",
    width = 1800, height = 600, res = 150)

# Set 4 panels
par(mfrow = c(1, 4),
    mar = c(5, 5, 4, 2),
    cex.main = 1,
    cex.lab = 0.7,
    cex.axis = 0.55,
    mgp = c(0.5, 0.5, 0),
    las = 0.8)

# Function to plot a diagram if present
plot_if_not_empty <- function(Diag, title_text) {
  if (!is.null(Diag$diagram) && nrow(Diag$diagram) > 0) {
    plot(Diag$diagram,
         main = "",
         diagLim = c(0, max(Diag$diagram[, 2:3], na.rm = TRUE)),
         barcode = FALSE)
    axis(1, lwd = 0.5, cex.axis = 0.55)
    box(lwd = 0.5)
    title(main = title_text, line = 1)
  } else {
    plot.new()
    title(main = paste(title_text, "\nNo topological features"), line = -1)
  }
}

# Plot for each group
plot_if_not_empty(Diag_1, "Persistence Diagram – HbA1c ≤ 6.5")
plot_if_not_empty(Diag_2, "Persistence Diagram – 6.5 < HbA1c ≤ 7.5")
plot_if_not_empty(Diag_3, "Persistence Diagram – 7.5 < HbA1c ≤ 8.5")
plot_if_not_empty(Diag_4, "Persistence Diagram – HbA1c > 8.5")

dev.off()

# Support function to calculate distance only if both diagrams have features
compute_bn_dist <- function(DiagA, DiagB, nameA, nameB) {
  if (nrow(DiagA$diagram) > 0 && nrow(DiagB$diagram) > 0) {
    dist <- bottleneck(DiagA$diagram, DiagB$diagram, dimension = 0)
    cat("Bottleneck distance (H₀) between", nameA, "and", nameB, ":", dist, "\n")
    return(dist)
  } else {
    cat("Bottleneck distance (H₀) between", nameA, "and", nameB, ": not computed (one or both empty)\n")
    return(NA)
  }
}

# Calculate the distances between all pairs
bn_12 <- compute_bn_dist(Diag_1, Diag_2, "HbA1c ≤ 6.5", "6.5 < HbA1c ≤ 7.5")
bn_13 <- compute_bn_dist(Diag_1, Diag_3, "HbA1c ≤ 6.5", "7.5 < HbA1c ≤ 8.5")
bn_14 <- compute_bn_dist(Diag_1, Diag_4, "HbA1c ≤ 6.5", "HbA1c > 8.5")
bn_23 <- compute_bn_dist(Diag_2, Diag_3, "6.5 < HbA1c ≤ 7.5", "7.5 < HbA1c ≤ 8.5")
bn_24 <- compute_bn_dist(Diag_2, Diag_4, "6.5 < HbA1c ≤ 7.5", "HbA1c > 8.5")
bn_34 <- compute_bn_dist(Diag_3, Diag_4, "7.5 < HbA1c ≤ 8.5", "HbA1c > 8.5")