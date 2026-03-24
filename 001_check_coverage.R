#BiocManager::install("dorothea")
library(dorothea)
library(dplyr)

larc_seurat <- readRDS("../193-SEURAT-HPC/larc_datasets/larc_merged_6k.rds")

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
png("coverage_histogram.png", width = max(3600, nrow(regulon_coverage) * 55), height = 2700, res = 330)
par(mar = c(8, 4, 4, 2))
barplot(regulon_coverage$n_targets_in_panel,
        names.arg = regulon_coverage$tf,
        main      = "TF Target Coverage in CosMX Panel",
        xlab      = "",
        ylab      = "Number of Targets Covered",
        las       = 2,
        cex.names = 0.6)
dev.off()