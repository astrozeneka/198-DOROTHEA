#BiocManager::install("dorothea")
#BiocManager::install("viper")
library(dorothea)
library(dplyr)
library(Seurat)

larc_seurat <- readRDS("../193-SEURAT-HPC/larc_datasets/larc_merged_6k.rds")

# ── Filter to annotated cells & tag source sample ─────────────────────────────
annotation_files <- list.files(
  "../191-SPLIT-COSMIX-SAMPLES/splitted/roi_annotated",
  pattern = "\\.csv$", full.names = TRUE
)

cell_sample_list <- lapply(annotation_files, function(f) {
  sname <- basename(f)
  sname <- sub("^LARC_[AB]_", "", sname)
  sname <- sub("_annotated\\.csv$", "", sname)
  df    <- read.csv(f, row.names = 1)
  data.frame(
    cell      = rownames(df)[df$selection_mask == 1],
    sample_id = sname,
    stringsAsFactors = FALSE
  )
})

cell_sample_df <- do.call(rbind, cell_sample_list)
larc_seurat    <- subset(larc_seurat, cells = cell_sample_df$cell)

# Add sample_id to each cell's metadata
sample_vec  <- setNames(cell_sample_df$sample_id, cell_sample_df$cell)
larc_seurat <- AddMetaData(larc_seurat, metadata = sample_vec[Cells(larc_seurat)], col.name = "sample_id")

# Step 0: Regulon Coverage Check
# Load regulons - confidence tiers A and B recommended
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))  # optionally add "C"

# Genes present in your CosMX panel
cosmx_genes <- rownames(larc_seurat)

# Check per-TF coverage
regulon_coverage <- regulons %>%
  group_by(tf) %>%
  summarise(
    n_targets_total    = n(),
    n_targets_in_panel = sum(target %in% cosmx_genes),
    pct_covered        = round(100 * n_targets_in_panel / n_targets_total, 1)
  ) %>%
  arrange(desc(n_targets_in_panel))

# Recommended filter: keep TFs with >=5 targets covered AND >=20% coverage
reliable_tfs <- regulon_coverage %>%
  dplyr::filter(n_targets_in_panel >= 5, pct_covered >= 20)

cat("Total TFs in A+B regulons:", n_distinct(regulons$tf), "\n")
cat("TFs passing coverage filter:", nrow(reliable_tfs), "\n")

# Filter regulons to only reliable TFs
regulons_filtered <- regulons %>%
  dplyr::filter(tf %in% reliable_tfs$tf,
                target %in% cosmx_genes)

# Barplot of target coverage with TF names on x-axis
png("coverage_histogram_annotated_samples.png", width = max(3600, nrow(regulon_coverage) * 55), height = 2700, res = 330)
par(mar = c(8, 4, 4, 2))
barplot(regulon_coverage$n_targets_in_panel,
        names.arg = regulon_coverage$tf,
        main      = "TF Target Coverage in CosMX Panel",
        xlab      = "",
        ylab      = "Number of Targets Covered",
        las       = 2,
        cex.names = 0.6)
dev.off()

# ── Step 1: Run VIPER to Compute TF Activity ──────────────────────────────────
library(viper)

# Extract normalized expression matrix (genes x cells)
expr_matrix <- GetAssayData(larc_seurat, assay = "RNA", layer = "data")

# Convert dorothea regulons to VIPER-compatible format
regulons_viper <- df2regulon(regulons_filtered)

# Run VIPER — output: matrix of TFs x cells (NES scores)
# >>> IMPOSSIBLE TO RUN DUE TO THE LARGE AMOUNT OF DATA >>>
tf_activity <- viper(
  eset       = as.matrix(expr_matrix),
  regulon    = regulons_viper,
  minsize    = 5,      # matches coverage filter above
  method     = "none", # use normalized expression as-is
  pleiotropy = TRUE,
  verbose    = TRUE
)
saveRDS(tf_activity, "rdata/tf_activity.rds")

# ── Reload from here to skip re-running VIPER ─────────────────────────────────
# tf_activity <- readRDS("rdata/tf_activity.rds")

# ── Step 2: Scale TF Activity Matrix ──────────────────────────────────────────
# Z-score each TF across cells (equivalent to ScaleData), used for differential testing
tf_activity_scaled <- t(scale(t(tf_activity)))

# ── Step 3: Differential TF Activity — Responders vs. Non-Responders ──────────
# Map sample_id → response via clinical data (same source as 401_volcano_with_annotated.R)
clinical     <- read.csv(
  "../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv",
  stringsAsFactors = FALSE
)
response_map <- setNames(clinical$pCR, clinical$SampleId)

cell_meta   <- larc_seurat@meta.data
cell_groups <- response_map[cell_meta$sample_id]
names(cell_groups) <- rownames(cell_meta)

unmatched <- names(cell_groups)[is.na(cell_groups)]
if (length(unmatched) > 0)
  message("Cells with no clinical match (excluded): ", length(unmatched))

responder_cells     <- names(cell_groups)[cell_groups == "Responder"]
non_responder_cells <- names(cell_groups)[cell_groups == "Non-responder"]

mat_resp    <- tf_activity_scaled[, responder_cells,     drop = FALSE]
mat_nonresp <- tf_activity_scaled[, non_responder_cells, drop = FALSE]

# Per-TF Wilcoxon test; avg_log2FC here is the difference in mean NES (appropriate for NES scores)
tf_diff <- do.call(rbind, lapply(rownames(tf_activity_scaled), function(tf) {
  x  <- mat_resp[tf, ]
  y  <- mat_nonresp[tf, ]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.frame(
    tf         = tf,
    avg_log2FC = mean(x) - mean(y),  # mean NES difference (responder - non-responder)
    p_val      = wt$p.value,
    stringsAsFactors = FALSE
  )
}))

tf_diff$p_val_adj <- p.adjust(tf_diff$p_val, method = "BH")
rownames(tf_diff) <- tf_diff$tf

tf_diff_sig <- tf_diff %>%
  dplyr::filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

print(tf_diff_sig)

# ── Step 4: Volcano Plot of Differential TF Activity ──────────────────────────
library(ggplot2)
library(ggrepel)

tf_diff$neg_log10_p <- -log10(tf_diff$p_val_adj + 1e-300)
tf_diff$significant  <- tf_diff$p_val_adj < 0.05 & abs(tf_diff$avg_log2FC) > 0.5
tf_diff$label        <- ifelse(tf_diff$significant, tf_diff$tf, NA)

ggplot(tf_diff, aes(x = avg_log2FC, y = neg_log10_p)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("grey70", "firebrick")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "grey40") +
  labs(
    title    = "Differential TF Activity: Responders vs. Non-Responders",
    subtitle = "DoRothEA + VIPER (A+B regulons)",
    x        = "Mean NES difference (responder - non-responder)",
    y        = "-log10(adjusted p-value)"
  ) +
  theme_classic() +
  theme(legend.position = "none")