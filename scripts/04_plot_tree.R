#!/usr/bin/env Rscript
# ============================================================================
# 04_plot_tree.R — Publication-quality circular phylogenetic tree with
#                  concentric heatmap rings showing protein identity
#
# Reads:  data/taxonomy.csv, data/metadata.csv
# Writes: output/IS621_phylo_tree.pdf
#
# Dependencies: ggtree, ape, ggplot2, ggnewscale, treeio, dplyr, tidyr
# ============================================================================

suppressPackageStartupMessages({
  library(ggtree)
  library(ape)
  library(ggplot2)
  library(ggnewscale)
  library(treeio)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------
# 1. Read data
# ------------------------------------------------------------------
taxonomy <- read.csv("data/taxonomy.csv", stringsAsFactors = FALSE)
metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE,
                      na.strings = c("NA", ""))

cat("Loaded", nrow(taxonomy), "taxa and", nrow(metadata), "metadata rows\n")

# ------------------------------------------------------------------
# 2. Build tree from taxonomy hierarchy
# ------------------------------------------------------------------
# Ensure all taxonomy columns are factors (required by as.phylo.formula)
tax_for_tree <- taxonomy %>%
  mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species),
                as.factor))

tree <- as.phylo.formula(
  ~Domain/Phylum/Class/Order/Family/Genus/Species,
  data = tax_for_tree
)

cat("Tree tips:", length(tree$tip.label), "\n")

# ------------------------------------------------------------------
# 3. Set up phylum coloring via groupOTU
# ------------------------------------------------------------------
# Map tip labels to phyla
tip_phylum <- setNames(metadata$Phylum, metadata$Tip_ID)

# Build phylum group list (only for tips present in tree)
tips_in_tree <- tree$tip.label
phyla_present <- unique(tip_phylum[tips_in_tree])
phyla_present <- phyla_present[!is.na(phyla_present)]

phylum_groups <- lapply(phyla_present, function(p) {
  names(tip_phylum[tip_phylum == p & names(tip_phylum) %in% tips_in_tree])
})
names(phylum_groups) <- phyla_present

tree_grouped <- groupOTU(tree, phylum_groups, group_name = "Phylum")

# ------------------------------------------------------------------
# 4. Curated phylum color palette
# ------------------------------------------------------------------
phylum_colors <- c(
  "Pseudomonadota"          = "#1f77b4",  # blue
  "Bacillota"               = "#ff7f0e",  # orange
  "Actinomycetota"          = "#2ca02c",  # green
  "Bacteroidota"            = "#d62728",  # red
  "Cyanobacteriota"         = "#17becf",  # teal
  "Spirochaetota"           = "#9467bd",  # purple
  "Thermodesulfobacteriota" = "#8c564b",  # brown
  "Deinococcota"            = "#e377c2",  # pink
  "Chloroflexota"           = "#7f7f7f",  # gray
  "Planctomycetota"         = "#bcbd22",  # olive
  "Verrucomicrobiota"       = "#aec7e8",  # light blue
  "Fusobacteriota"          = "#ffbb78",  # light orange
  "Chlamydiota"             = "#98df8a",  # light green
  "Euryarchaeota"           = "#ff9896",  # light red
  "Thermoproteota"          = "#c5b0d5",  # light purple
  "Asgardarchaeota"         = "#c49c94"   # light brown
)
# Add a fallback for any phyla not in the palette
all_phyla <- levels(attr(tree_grouped, "group"))
missing_phyla <- setdiff(all_phyla, c(names(phylum_colors), "0"))
if (length(missing_phyla) > 0) {
  extra_colors <- rainbow(length(missing_phyla))
  names(extra_colors) <- missing_phyla
  phylum_colors <- c(phylum_colors, extra_colors)
}
# The "0" group (unassigned) gets dark grey
phylum_colors <- c("0" = "#444444", phylum_colors)

# ------------------------------------------------------------------
# 5. Create circular cladogram
# ------------------------------------------------------------------
p <- ggtree(tree_grouped,
            layout = "circular",
            branch.length = "none",
            aes(color = Phylum),
            size = 0.4) +
  scale_color_manual(
    values = phylum_colors,
    name = "Phylum",
    breaks = names(phylum_colors)[names(phylum_colors) != "0"],
    guide = guide_legend(ncol = 1, override.aes = list(size = 3))
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# ------------------------------------------------------------------
# 6. Prepare heatmap data
# ------------------------------------------------------------------
heatmap_df <- metadata[, c("Tip_ID", "EF_Tu_Identity", "SecY_Identity",
                           "MuA_Identity", "IS621_Identity")]
rownames(heatmap_df) <- heatmap_df$Tip_ID
heatmap_df$Tip_ID <- NULL

# Rename columns for display
colnames(heatmap_df) <- c("EF-Tu", "SecY", "MuA", "IS621")

# Ensure only tree tips are in the heatmap
heatmap_df <- heatmap_df[rownames(heatmap_df) %in% tips_in_tree, , drop = FALSE]

# ------------------------------------------------------------------
# 7. Add concentric heatmap rings
# ------------------------------------------------------------------
p2 <- gheatmap(
  p,
  heatmap_df,
  offset = 1,
  width = 0.5,
  colnames_angle = 95,
  colnames_offset_y = 0.4,
  colnames_position = "top",
  color = NA,     # no cell borders
  hjust = 0,
  font.size = 3
) +
  scale_fill_gradient(
    low = "white",
    high = "#08306B",
    na.value = "grey80",
    name = "% AA Identity\nto E. coli",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 6,
      title.position = "top"
    )
  )

# ------------------------------------------------------------------
# 8. Highlight E. coli
# ------------------------------------------------------------------
ecoli_tip <- metadata$Tip_ID[metadata$Ecoli_Reference == "TRUE"]
if (length(ecoli_tip) > 0) {
  ecoli_label <- ecoli_tip[1]
  p2 <- p2 +
    geom_tippoint(
      aes(subset = (label == ecoli_label)),
      color = "red",
      size = 3,
      shape = 18  # diamond
    ) +
    geom_tiplab2(
      aes(subset = (label == ecoli_label), label = "E. coli"),
      color = "red",
      size = 3,
      fontface = "bold.italic",
      offset = 8,
      hjust = 0
    )
}

# ------------------------------------------------------------------
# 9. Final styling and export
# ------------------------------------------------------------------
p_final <- p2 +
  ggtitle("Distribution of IS621 and Conserved Proteins\nAcross the Prokaryotic Tree of Life") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.box = "vertical",
    legend.spacing.y = unit(0.3, "cm")
  )

# Save
dir.create("output", showWarnings = FALSE)
ggsave(
  "output/IS621_phylo_tree.pdf",
  plot = p_final,
  width = 14,
  height = 14,
  device = cairo_pdf
)

cat("Saved: output/IS621_phylo_tree.pdf\n")

# Also save PNG for quick preview
ggsave(
  "output/IS621_phylo_tree.png",
  plot = p_final,
  width = 14,
  height = 14,
  dpi = 300
)

cat("Saved: output/IS621_phylo_tree.png\n")
