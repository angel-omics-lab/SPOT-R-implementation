library("readxl")
library("stringr")
library("LaplacesDemon")
library("tidyr")
library("plyr")
library("Matrix")
library("ggplot2")
library("RColorBrewer")
library("gridExtra")
library("gridGraphics")
library("igraph")
library("ggdendro")

#--spatial packages----
library("sf")
library("spatstat")
library("spdep")
library("spatialEco")
library("spatialreg")

#### some functions #####
range_std <- function(x, lowest, highest){(x - lowest)/(highest - lowest)}
bw_calculator <- function(x){
  0.9 * min(sd(x, na.rm = T), IQR(x, na.rm = T)/1.34) / length(x)^(1/5)
}


#### sample A5 #####
ROI_names = str_remove(excel_sheets(path = "/Users/sealso/Research/Collaboration/Peggi/Data/DCIS065_A5_HE_total.xlsx"), "ROI")

#### other samples ####
#ROI_names2 = str_remove(excel_sheets(path = "/Users/sealso/Research/Collaboration/Peggi/Data/20230816_DCIS091_A7_Export.xlsx"), "ROI")
#ROI_names3 = str_remove(excel_sheets(path = "/Users/sealso/Research/Collaboration/Peggi/Data/DCIS_003_B5_C3_total.xlsx"), "ROI")
#ROI_names4 = str_remove(excel_sheets(path = "/Users/sealso/Research/Collaboration/Peggi/Data/DCIS_024_A29_C3_total.xlsx"), "ROI")

protein_groups<-readxl::read_excel("/Users/sealso/Research/Collaboration/Peggi/Data/DCIS C3 Proteomic  IDs_pruned.xlsx")
protein_groups_two<-as.data.frame(protein_groups)[,-2]


full_dat <- check_zeros <- NULL
dat_frames <- list()
ROI_labels <- data.frame(ROI = ROI_names, 
                lesion = c("DCIS", "DCIS", "DCIS" , "IBC", "DCIS", "IBC", 
                "DCIS", rep("IBC", 13), rep("DCIS", 4), "IBC", "DCIS",
                "DCIS", "Normal"))

for(ROI in 1:length(ROI_names)){
  
  # loading the data for one ROI
  my_data <- readxl::read_excel("/Users/sealso/Research/Collaboration/Peggi/Data/DCIS065_A5_HE_total.xlsx",
                                sheet = ROI, skip = 8)
  # removing two unneeded columns
  df <- as.data.frame(my_data[, -c(1, 4)]) 
  
  # extracting m/z values as peptide names
  peptide_names <- as.character(round(as.numeric(names(df[, -c(1, 2)])), digits = 9))
  colnames(df) <- c("x", "y", peptide_names)

  # checking the proportion of NA values for each peptide
  NA_prop <- apply(df[, -c(1, 2)], 2, function(x) mean(is.na(x))) 
  
  # removing an ROI if it has more than 25% peptides with more than 25% NA values (slightly strong threshold for safety)
  if(mean(NA_prop > 0.25) > 0.25){
    dat.frame <- NULL  
  }else{
    df[is.na(df)] <- 0      # forcing NA's to be 0     
    dat_frames[[ROI]] <- df # raw ROI-level data
    dat.frame <- data.frame(ROI = paste0(ROI_labels[ROI, 2], 
                 "_A5_", ROI_names[ROI]), df[, 1:2], log1p(df[, -c(1:2)])) # ROI-level data after log-transformation
    zerof <- data.frame(ROI = paste0(ROI_labels[ROI, 2], "_A5_", ROI_names[ROI]), 
              t(apply(df[, -c(1, 2)], 2, function(x) mean(x==0)))) # storing proportion of zero values per peptide
    
  }
  full_dat <- plyr::rbind.fill(full_dat, dat.frame) # appending ROI-level transformed data
  check_zeros <- plyr::rbind.fill(check_zeros, zerof) # appending proportion of zero values per peptide per ROI
}

colnames(full_dat) <- c("ROI", "x", "y", peptide_names)
ComplexHeatmap::pheatmap(cor(full_dat[, -c(1:3)]), cluster_rows = F, cluster_cols = F)


# checking sparsity per peptide in the full data
sparsity <- array(NA, dim = c(length(peptide_names), 1))
for(g in 1:length(peptide_names)){
  sparsity[g] <- mean(check_zeros[, (g + 1)] > 0.2) # proportion of ROIs with % 0 of peptide g > 20%
}
good_peptides <- which(sparsity < 0.1) # selecting peptides with less than 10% ROIs having > 20% 0 values

# median value per ROI per peptide
ROI_all <- unique(full_dat$ROI)
median_vals <- array(NA, dim = c(length(peptide_names), length(ROI_all)))
for(ROI in 1:length(ROI_all)){
  for(g in 1:length(peptide_names)){
    pep_df <-  full_dat[full_dat$ROI == ROI_all[ROI], c("x", "y", peptide_names[g])]
    median_vals[g, ROI] <-  median(pep_df[, - c(1:2)])
  }
}
rownames(median_vals) <-  peptide_names
colnames(median_vals) <-  ROI_all

# differential intensity test using Kruskal-Wallis test (ignoring the "normal" ROI)
KT <- array(NA, dim = c(length(peptide_names), 2))
med_diff <- array(NA, dim = c(length(peptide_names), 1))
ROI_types <- sub("_.*", "", ROI_all)
which(ROI_types == "Normal")

for(g in  1:length(peptide_names)){
  ANOVA_data <- data.frame(peptide = median_vals[g, ], group = lesion_types)
  KTest <- kruskal.test(peptide ~ group, data = ANOVA_data[-c(28), ])
  KT[g, ] <- c(KTest$statistic, KTest$p.value)
  med_diff[g] <- median(median_vals[g, ][which(lesion_types == "DCIS")]) - 
    median(median_vals[g, ][which(lesion_types == "IBC")])
}
colnames(KT) <- c("test.statistic", "p.value")

median_KT_res <- data.frame(mean_diff = med_diff[good_peptides], KT[good_peptides, ], 
                            m.z = peptide_names[good_peptides])
median_KT_res <- merge(median_KT_res, protein_groups, by.x = "m.z", by.y = "m/z")
median_KT_res <- median_KT_res[!duplicated(median_KT_res$mean_diff), ] 
median_KT_res$adj_pvalue <- p.adjust(median_KT_res$p.value, "BH")


# filter significant peptides
sig_peptides <- median_KT_res$m.z[median_KT_res$adj_pvalue < 0.05]

# create long-format data for plotting
plot_data <- data.frame()

for (pep in sig_peptides) {
  temp <- data.frame(
    peptide = pep,
    value = median_vals[pep, ][-28], 
    group = lesion_types[ - 28]
  )
  plot_data <- rbind(plot_data, temp)
}

# convert factors
plot_data$group <- factor(plot_data$group, levels = c("DCIS", "IBC"))

col_pal <- scales::hue_pal()(3)[c(3, 1, 2)]



png(paste0("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Diff_exp_ROI_level",
           ".png"), height = 1200, width = 1600, res = 120)
# boxplots side by side for each peptide
ggplot(plot_data, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) + scale_fill_manual(values = col_pal[1:2]) +
  facet_wrap(~ peptide, scales = "free_y", nrow = 6, ncol = 6) +
  theme_bw(base_size = 12) + #ylim(c(4, 10)) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 10), plot.title = element_text(hjust = 0.5),
  ) +
  labs(
    x = "Group",
    y = "Mean peptide intensity",
    title = "Significant peptides based on an ROI-level test (adj. p-value < 0.05)"
  )

dev.off()



write.csv(median_KT_res, row.names = F,
          "/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Differential_expression.csv")



######################################
###### Spatial median plots ##########

g <- 1  
chk <- med_coords <- spatial_median <- NULL
spatial_median <- array(NA, dim = c(length(ROI_all), (length(peptide_names) + 2)))

for(ROI in 1:length(ROI_all)){
  pep_df <-  as.data.frame(full_dat[full_dat$ROI == ROI_all[ROI], -1])
  spatial_median[ROI, ] <-  colMeans(pep_df) # centroid location of each ROI
}
colnames(spatial_median) <- c("x", "y",  paste0("peptide: ", peptide_names))

df <- as.data.frame(spatial_median[, c("x", "y", paste0("peptide: ", sig_peptides))])
df[, -c(1:2)] <- apply(df[, -c(1:2)], 2, scale)

library(ggplot2)
library(tidyverse)

df_long <- df %>% pivot_longer(cols = starts_with("peptide: "),
                               names_to = "peptide", values_to = "intensity")

expression_min <- -1.96
expression_max <- 1.96
myPalette <- viridis::turbo(10)
sc_LH <- scale_colour_gradientn(colours = myPalette, 
                                limits = c(expression_min, expression_max),
                                breaks = seq(expression_min, expression_max, length.out = 9))
df_long$intensity <- pmax(pmin(df_long$intensity, expression_max), expression_min)
df_long$peptide <- gsub("^peptide:", "", df_long$peptide)
df_long$peptide <- factor(df_long$peptide, levels = unique(df_long$peptide))

png(paste0("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Mean_peptide_",
           ".png"), height = 1400, width = 1600, res = 140)
ggplot(df_long, aes(x = x, y = y, color = intensity)) +
  geom_point(size = 2) + 
  sc_LH + theme_minimal() + facet_wrap(~ peptide, ncol = 6) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean peptide intensities in each ROI", color = "Scaled \nIntensity")
dev.off()


############################################################
# Hiearchical clustering based on only significant peptides
############################################################

mat <- t(median_vals[good_peptides[which(median_KT_res$p.value < 0.05)], ])
rownames(mat) <-   gsub("A5_", "", gsub("^(DCIS|IBC|Normal)_", "", ROI_all))
dist_median <- dist(mat, method = "euclidean")
hc_median <- hclust(dist_median, "ward.D2")


library(dendextend)

# Convert to dendrogram object
dend <- as.dendrogram(hc_median)
# Map disease groups to colors
col_pal <- scales::hue_pal()(3)[c(3, 1, 2)]
label_colors <- ifelse(lesion_types == "IBC", col_pal[2], 
                       ifelse(lesion_types == "DCIS", col_pal[1], col_pal[3]))

# Apply colors to labels
labels_colors(dend) <- label_colors[order.dendrogram(dend)]

# Cut into 4 clusters and assign colors
dend <- color_branches(dend, k = 3, col = c("red4",  "red2", "green4"))


png(paste0("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/ROI_clusters_onlyA5",
           ".png"), height = 900, width = 1000, res = 140)
plot(dend, main = "Clustering ROIs based on median peptide expression")
legend("topright", legend = c("DCIS", "IBC", "Normal"), title = "Leaves", 
       fill = col_pal, border = NA, inset = c(0.00, 0.00))
legend("topright",
       legend = paste("Cluster", 1:3), title = "Branches",
       fill = c("red4",  "red2", "green4"),
       border = NA, inset = c(0.00, 0.22))
dev.off()



###### Random forest ########
library(randomForest)
library(datasets)
library(caret)
mat <- t(median_vals[good_peptides[which(median_KT_res$p.value < 0.05)], ])
rownames(mat) <-  ROI_all

rf_data <- as.data.frame(mat)
rf_data$labels <- lesion_types
rf_data <- rf_data[rf_data$labels != "Normal", ]
rf_data$labels <- factor(rf_data$labels)

colnames(rf_data)[1:(ncol(rf_data) - 1)] <- paste0("var", 1:(ncol(rf_data) - 1))

rf <- randomForest(labels~., data = rf_data, proximity= TRUE, ntree = 1000) 
vimp <- varImp(rf)
row.names(vimp) <-  colnames(mat)

library(tibble)
# Convert to data frame and strip "peptide:" from row names
vimp_df <- as.data.frame(vimp) %>% rownames_to_column(var = "peptide") 

# Plot
png(paste0("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Peptide_importance_RF_",
           ".png"), height = 800, width = 900, res = 140)

ggplot(vimp_df, aes(x = reorder(peptide, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +   theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Peptide (m/z)", y = "Variable Importance", 
       title = paste0("Peptide importance in an RF model with OOB error rate: ",
                      round(100*(sum(rf$confusion) - sum(diag(rf$confusion)))/sum(rf$confusion), 2), " %")
  )

dev.off()


############ Cardinal package ##############
############################################

library(Cardinal)
combined_intensities <- apply(full_dat[, -c(1:3)], 2, scale)
combined_coords <- apply(full_dat[, c(2, 3)], 2, as.numeric)
colnames(combined_coords) <- c("x", "y")


combined_intensities <- combined_intensities[, median_KT_res$m.z[median_KT_res$p.value< 0.05]]
sample_labels <-  gsub("^(DCIS|IBC|Normal)_", "", full_dat[, 1])
g_labels <-  sub("_.*", "", full_dat[, 1])


# Step 2: Create an MSImagingExperiment object
msi_exp <- MSImagingExperiment(spectraData = t(combined_intensities), 
                               pixelData = PositionDataFrame(coord = combined_coords, 
                                                             run = factor(sample_labels), type = factor(g_labels)))

png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Sanity_check2.png", 
    height = 2000, width = 2000, res = 100)  # Double width for side-by-side
image(msi_exp, i = 1 ,free = "xy")
dev.off()

#  the total ion current
msi_exp <- summarizePixels(msi_exp, stat=c(TIC="sum"))
image(msi_exp, "TIC", free = "xy")

# PCA
msi_pca <- PCA(msi_exp, ncomp = 2)
image(msi_pca, type="x", free = "xy", scale = TRUE)
plot(msi_pca, type="x", groups = msi_exp$type, shape=20)

# NMF
msi_nmf <- NMF(msi_exp, ncomp = 2, niter = 5)
image(msi_nmf, type="x", free = "xy", scale = TRUE)

png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/PCA.png", 
    height = 600, width = 600, res = 100)  # Double width for side-by-side
plot(msi_pca, type = "x", groups = msi_exp$type, shape = 20)
title("PCA")
dev.off()

png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/NMF.png", 
    height = 600, width = 600, res = 100)  # Double width for side-by-side
plot(msi_nmf, type = "x", groups = msi_exp$type, shape = 20)
title("NMF")
dev.off()


png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/NMF_im.png", 
    height = 3200, width = 3200, res = 200)  # Double width for side-by-side
image(msi_nmf, type="x", free = "xy", scale = TRUE, 
      col = hcl.colors(50, "Inferno", rev = TRUE)[c(50, 1)])
dev.off()

# differential intensity test using cardinal (based on mean)
rcc_mtest <- meansTest(msi_exp, ~type)
rcc_mtest
rcc_mtest_top <- topFeatures(rcc_mtest)
plot(rcc_mtest)
png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Mean_diff.png", 
    height = 600, width = 900, res = 200)  # Double width for side-by-side
plot(rcc_mtest, i=c("m/z = 1098.5902" = 1, "m/z = 1588.7813" = 2,  
                    "m/z = 1458.7006" = 3, "m/z = 1286.66" = 4, 
                    "m/z = 1226.6124" = 5, "m/z = 1767.9236" = 6))
dev.off()


########## SingleCellExperiment package ############
###################################################

library(SingleCellExperiment)
combined_intensities <- full_dat[, -c(1:3)]
combined_intensities <- apply(combined_intensities, 2, scale)
combined_coords <- full_dat[, c(2, 3)]
colnames(combined_coords) <- c("x", "y")
combined_intensities <- combined_intensities[, median_KT_res$m.z[median_KT_res$p.value < 0.05]]

# Step 1: Transpose intensities: features (peptides) x pixels
expr_mat <- t(as.matrix(combined_intensities))  # peptides x pixels
sample_labels <-  gsub("^(DCIS|IBC|Normal)_", "", full_dat[, 1])
g_labels <-  sub("_.*", "", full_dat[, 1])

# Step 2: Create pixel-level metadata
pixel_metadata <- DataFrame(
  sample = sample_labels,
  group = g_labels,           
  x = combined_coords[, "x"],         # x-coordinates
  y = combined_coords[, "y"]          # y-coordinates
)

# Step 3: Create peptide metadata (optional)
peptide_metadata <- DataFrame(
  peptide = rownames(expr_mat)        # peptide names as rownames
)

# Step 4: Construct SingleCellExperiment
sce <- SingleCellExperiment(
  assays = list(counts = expr_mat),
  colData = pixel_metadata,
  rowData = peptide_metadata
)

# Optional: Store spatial coordinates as a reduced dimension
reducedDims(sce) <- SimpleList(spatial = as.matrix(combined_coords))
colLabels(sce) <- sce$sample  # or sce$group, depending on your goal

library(scater)
# PCA at the pixel level
sce <- scater::runPCA(sce, exprs_values = "counts",  ncomponents = 5)
reducedDimNames(sce)
plotPCA(sce, colour_by = "group")  # or "sample", etc.

library(dplyr)
plot_sce_pca <- function(sce, colour_by = "group", pcs = c(1,2), legend_title = "Lesion type") {
  # Extract PCA coordinates
  pca_coords <- reducedDim(sce, "PCA")
  
  # Make main data.frame
  df <- data.frame(
    PC1 = pca_coords[, pcs[1]],
    PC2 = pca_coords[, pcs[2]],
    ColourVar = factor(colData(sce)[[colour_by]])
  )
  
  custom_colors <- c("#377EB8", "#FF8C00", "#4DAF4A")  # blue, orange, green

  # PCA variance for axes
  pca_var <- attr(reducedDim(sce, "PCA"), "percentVar")
  pc1_var <- round(pca_var[pcs[1]], 1)
  pc2_var <- round(pca_var[pcs[2]], 1)
  
  ggplot(df, aes(x = PC1, y = PC2, color = ColourVar)) +
    geom_point(size = 1) +
    labs(
      x = paste0("PC", pcs[1], " (", pc1_var, "% variance)"),
      y = paste0("PC", pcs[2], " (", pc2_var, "% variance)"),
      colour = legend_title
    ) +
    theme_minimal() +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18)
    ) +
    guides(color = guide_legend(override.aes = list(size = 6)))
}

# UMAP at the pixel level
sce <- scater::runUMAP(sce, exprs_values = "counts", dimred = "PCA", n_dimred  = 2)
colLabels(sce) <- sce$sample  # or sce$group, depending on your goal
by.cluster <- aggregateAcrossCells(sce, ids=colLabels(sce))
centroids <- reducedDim(by.cluster, "PCA")

# MST based pseudotime by TSCAN package
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst
col_pal <- scales::hue_pal()(3)[c(3, 1, 2)]
names(col_pal) <- c("DCIS", "IBC", "Normal")  # <-- Make sure these match levels(sce$group)

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")

library("ggrepel")
png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/PCA_MST_UMAP.png", 
    height = 1000, width = 1200, res = 100)  # Double width for side-by-side

umap_coords <- reducedDim(sce, "UMAP")
umap_df <- data.frame(
  sample = sce$sample, group = sce$group, UMAP1 = umap_coords[, 1], 
  UMAP2 = umap_coords[, 2]
)

# Compute centroid per sample
centroids <- aggregate(cbind(UMAP1, UMAP2) ~ sample + group, data = umap_df, FUN = mean)
line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP")

# Plot
plotUMAP(sce, colour_by = "group", point_size = 0.2) + 
  geom_point(data = centroids, aes(x = UMAP1, y = UMAP2, color = group),
             size = 6) + 
  geom_line(data = line.data, aes(x = UMAP1, y = UMAP2, group = edge),
            color = "black", linewidth = 1) +
  geom_text_repel(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, label = sample),
    size = 6, fontface = "bold", color = "black"
  ) +   scale_x_reverse() +   # <-- mirror the x-axis
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) + guides(color = guide_legend(title = "Lesion type", override.aes = list(size = 6)))
ggtitle("UMAP with MST connecting the ROI centroids") 
dev.off()

map.tscan <- mapCellsToEdges(sce, mst=mst, use.dimred="PCA")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)
common.pseudo <- averagePseudotime(tscan.pseudo)
png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/PCA_MST_UMAP_pseudotime.png", 
    height = 1400, width = 2000, res = 100)  # Double width for side-by-side


plotUMAP(sce, colour_by = I(common.pseudo), point_size = 0.2) +
  geom_point(data = centroids, aes(x = UMAP1, y = UMAP2), size = 6) +
  geom_line(data = line.data, aes(x = UMAP1, y = UMAP2, group = edge),
            color = "black", linewidth = 1)  +
  geom_text_repel(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, label = sample),
    size = 6, fontface = "bold", color = "black"
  ) + 
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()


# pseudotime estimation using slingshot package (could be slow)
library(slingshot)
sce.sling <- slingshot(sce, reducedDim='PCA')
head(sce.sling$slingPseudotime_1)
embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])

png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/Slingshot_UMAP_pseudotime.png", 
    height = 1000, width = 1200, res = 100)  # Double width for side-by-side

custom_colors <- c("#377EB8", "#FF8C00", "#4DAF4A")  # blue, orange, green

plotUMAP(sce.sling, colour_by = "slingPseudotime_1", point_size = 0.2) +
  geom_path(data = embedded, aes(x=UMAP1, y=UMAP2), size = 1) + 
  geom_point(data = centroids, aes(x = UMAP1, y = UMAP2), size = 3) +
  geom_text_repel(
    data = centroids,
    aes(x = UMAP1, y = UMAP2, label = sample),  # map group to color
    size = 6, color = "black",
    fontface = "bold") +  scale_x_reverse() + 
  scale_color_viridis_c(name = "Pseudotime") +  # <--- set legend title here
  ggtitle("UMAP with Slingshot Pseudotime and ROI centroids") +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))


dev.off()



### SingleCellExperiment object using summary data (ROI-level)
##############################################################

require(SingleCellExperiment)
median_coords <- array(NA, dim = c(length(ROI_all), 2))

for(ROI in 1:length(ROI_all)){
  pep_df <-  full_dat[full_dat$ROI == ROI_all[ROI], c("x", "y")]
  median_coords[ROI, ] <-  apply(pep_df[,  c(1:2)], 2, median)
}

# Hierachical clustering based on only significant peptides
# Step 1: Transpose intensities: features (peptides) x pixels
expr_mat <- t(mat)  # peptides x pixels

# Step 2: Create pixel-level metadata
pixel_metadata <- DataFrame(
  sample = factor(ROI_labels[, 1]),     # sample ID for each pixel
  group = factor(ROI_labels[, 2]),           # group label (e.g., SSL, TVA/VA)
  x =  median_coords[, 1],         # x-coordinates
  y =  median_coords[, 2]          # y-coordinates
)

# Step 3: Create peptide metadata (optional)
peptide_metadata <- DataFrame(
  peptide = rownames(expr_mat)        # peptide names as rownames
)

# Step 4: Construct SingleCellExperiment
sce2 <- SingleCellExperiment(
  assays = list(counts = expr_mat),
  colData = pixel_metadata,
  rowData = peptide_metadata
)

# Optional: Store spatial coordinates as a reduced dimension
# reducedDims(sce) <- SimpleList(spatial = as.matrix(combined_coords))

colLabels(sce2) <- sce2$sample  # or sce$group, depending on your goal

library(scater)
library(ggplot2)
library(ggrepel)  # for better text labeling

sce2 <- logNormCounts(sce2)
sce2 <- scater::runPCA(sce2, exprs_values = "counts",  ncomponents = 5)
reducedDimNames(sce2)

png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/PCA_on_mean_rotated.png", 
    height = 800, width = 1000, res = 150)  # Double width for side-by-side

plot_pca <- function(sce, colour_by = "group", pcs = c(1,2), legend_title = "Lesion type") {
  # Extract PCA coordinates
  pca_coords <- reducedDim(sce, "PCA")
  
  # Make main data.frame
  df <- data.frame(
    PC1 = pca_coords[, pcs[1]],
    PC2 = pca_coords[, pcs[2]],
    ColourVar = factor(colData(sce)[[colour_by]]),
    sample = colData(sce)[["sample"]]
  )
  
  custom_colors <- c("#377EB8", "#FF8C00", "#4DAF4A")  # blue, orange, green
  
  # PCA variance for axes
  pca_var <- attr(reducedDim(sce, "PCA"), "percentVar")
  pc1_var <- round(pca_var[pcs[1]], 1)
  pc2_var <- round(pca_var[pcs[2]], 1)
  
  ggplot(df, aes(x = -PC2, y = -PC1, color = ColourVar)) +
    geom_point(size = 3) + geom_text_repel(aes(label = sample), size = 6, max.overlaps = 20) +
    labs(
      x = paste0("PC", pcs[2], " (", pc2_var, "% variance)"),
      y = paste0("PC", pcs[1], " (", pc1_var, "% variance)"),
      colour = legend_title
    ) + ggtitle("PCA on the mean peptide expression") +
    theme_minimal() +
    scale_color_manual(values = custom_colors) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12), 
          axis.title = element_text(size = 15)) + 
    guides(color = guide_legend(override.aes = list(size = 4)))
}

plot_pca(sce2)
dev.off()

sce2 <- scater::runUMAP(sce2, exprs_values = "counts", dimred = "PCA", n_dimred  = 5)

png("/Users/sealso/Research/Collaboration/Peggi/Results/Plots/DCIS065_A5/Aug7/UMAP_on_mean_based_PCs.png", 
    height = 800, width = 1000, res = 150)  # Double width for side-by-side

plot_umap <- function(sce, colour_by = "group", legend_title = "Lesion type") {
  # Extract UMAP coordinates
  umap_coords <- reducedDim(sce, "UMAP")
  
  # Make main data.frame
  df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    ColourVar = factor(colData(sce)[[colour_by]]),
    sample = colData(sce)[["sample"]]
  )
  
  custom_colors <- c("#377EB8", "#FF8C00", "#4DAF4A")  # blue, orange, green
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = ColourVar)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = sample), size = 6, max.overlaps = 20) +
    labs(
      x = "UMAP1",
      y = "UMAP2",
      colour = legend_title
    ) + 
    ggtitle("UMAP on the mean peptide expression") +
    theme_minimal() +
    scale_color_manual(values = custom_colors) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12), 
      axis.title = element_text(size = 15)
    ) + 
    guides(color = guide_legend(override.aes = list(size = 4)))
}

plot_umap(sce2)
dev.off()

