#!/usr/bin/env Rscript
# =============================================================================
# microbiome_plots.R
# Comprehensive microbiome visualization for amptk amplicon data (16S / ITS)
#
# Usage:
#   Rscript microbiome_plots.R \
#     --input  BEEHONEY_16S_ASV \
#     --metadata metadata/16S_metadata.tsv \
#     --outdir  results/16S_ASV_plots
#
# Outputs:
#   *_sample_depth.pdf/png          -- per-sample read counts
#   *_rarefaction.pdf/png           -- rarefaction curves
#   *_alpha_diversity.pdf/png       -- Shannon/Simpson/Chao1/Observed by host
#   *_alpha_by_material.pdf/png     -- Shannon by honey vs pollen
#   *_barplot_Phylum/Class/Order/Family/Genus.pdf/png
#   *_microshade_phylum_genus.pdf/png  -- microshade sub-grouped bars
#   *_beta_bray_pcoa.pdf/png        -- Bray-Curtis PCoA
#   *_beta_jaccard_pcoa.pdf/png     -- Jaccard PCoA
#   *_beta_bray_nmds.pdf/png        -- NMDS (Bray-Curtis)
#   *_beta_wunifrac_pcoa.pdf/png    -- Weighted UniFrac (if tree present)
#   *_beta_unifrac_pcoa.pdf/png     -- Unweighted UniFrac (if tree present)
#   *_beta_combined.pdf             -- combined panel
#   *_heatmap_genus.pdf/png         -- top-genus abundance heatmap
#   *_honey_vs_pollen_MA.pdf/png    -- MA plot honey vs pollen
#   *_phyloseq.rds                  -- full phyloseq object
#   *_phyloseq_rarefied.rds         -- rarefied phyloseq object
# =============================================================================

# ---- package loading --------------------------------------------------------
suppressPackageStartupMessages({
  required_pkgs <- c(
    "optparse", "phyloseq", "vegan", "ggplot2",
    "dplyr", "tidyr", "tibble", "patchwork",
    "RColorBrewer", "scales", "ggrepel"
  )
  missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))
})

# microshade: optional, from GitHub
HAS_MICROSHADE <- requireNamespace("microshade", quietly = TRUE)
if (!HAS_MICROSHADE) {
  message("microshade not found – attempting install from GitHub (KasperSkytte/microshade)...")
  tryCatch({
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", quiet = TRUE)
    remotes::install_github("KasperSkytte/microshade", quiet = TRUE)
    library(microshade)
    HAS_MICROSHADE <- TRUE
  }, error = function(e) {
    message("microshade unavailable; will use manual shaded palette instead.")
  })
} else {
  library(microshade)
}

# ---- argument parsing -------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"),      type = "character", default = NULL,
              help = "Amptk output folder (e.g. BEEHONEY_16S_ASV)", metavar = "DIR"),
  make_option(c("-m", "--metadata"),   type = "character", default = NULL,
              help = "Metadata TSV (sample-id, Country, host, material, ...)", metavar = "FILE"),
  make_option(c("-o", "--outdir"),     type = "character", default = "results/plots",
              help = "Output directory [default: results/plots]", metavar = "DIR"),
  make_option(c("-t", "--type"),       type = "character", default = NULL,
              help = "Marker: 16S or ITS (auto-detected from folder name)", metavar = "TYPE"),
  make_option(c("--table_type"),       type = "character", default = NULL,
              help = "Table type: ASV or cluster (auto-detected if omitted)", metavar = "TYPE"),
  make_option(c("--min_reads"),        type = "integer",   default = 1000,
              help = "Min reads per sample to retain [default: 1000]", metavar = "INT"),
  make_option(c("--top_taxa"),         type = "integer",   default = 15,
              help = "Top taxa shown in bar/heatmap plots [default: 15]", metavar = "INT"),
  make_option(c("--rarefy"),           type = "integer",   default = NULL,
              help = "Rarefaction depth (default: auto = min passing sample)", metavar = "INT"),
  make_option(c("--width"),            type = "numeric",   default = 12,
              help = "Plot width in inches [default: 12]", metavar = "NUM"),
  make_option(c("--height"),           type = "numeric",   default = 8,
              help = "Plot height in inches [default: 8]", metavar = "NUM")
)

opt_parser <- OptionParser(option_list = option_list,
                            description = "Microbiome plots from amptk amplicon results")
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$metadata)) {
  print_help(opt_parser)
  stop("--input and --metadata are required.", call. = FALSE)
}

# auto-detect marker type
if (is.null(opt$type)) {
  opt$type <- dplyr::case_when(
    grepl("16S", opt$input, ignore.case = TRUE) ~ "16S",
    grepl("ITS", opt$input, ignore.case = TRUE) ~ "ITS",
    TRUE ~ "amplicon"
  )
  message("Auto-detected marker type: ", opt$type)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- rank definitions -------------------------------------------------------
TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
RANK_PREFIXES <- c(k = "Kingdom", p = "Phylum", c = "Class",
                   o = "Order",   f = "Family", g = "Genus", s = "Species")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Parse a vector of amptk taxonomy strings into a data.frame
#' Format: "SS|score|accession;k:X,p:X,c:X,o:X,f:X,g:X"
parse_amptk_taxonomy <- function(tax_strings) {
  parse_one <- function(s) {
    result <- setNames(rep(NA_character_, length(TAX_RANKS)), TAX_RANKS)
    s <- trimws(as.character(s))
    # Split on first ; to get the taxonomy portion
    parts <- strsplit(s, ";", fixed = TRUE)[[1]]
    if (length(parts) < 2) return(result)
    tax_part <- paste(parts[-1], collapse = ";")  # rejoin in case of nested ;
    if (grepl("No hit|no hit", tax_part)) return(result)
    fields <- strsplit(trimws(tax_part), ",")[[1]]
    for (f in fields) {
      f <- trimws(f)
      if (grepl("^[kpcofgs]:", f)) {
        prefix <- substr(f, 1, 1)
        val    <- sub("^[kpcofgs]:", "", f)
        val    <- gsub('"', '', val)   # some amptk names include quotes
        val    <- trimws(val)
        if (prefix %in% names(RANK_PREFIXES) && nchar(val) > 0) {
          result[RANK_PREFIXES[prefix]] <- val
        }
      }
    }
    result
  }
  tax_mat <- do.call(rbind, lapply(tax_strings, parse_one))
  rownames(tax_mat) <- names(tax_strings)
  as.data.frame(tax_mat, stringsAsFactors = FALSE)
}

#' Load amptk OTU/ASV taxonomy table; returns list(otu_mat, tax_df)
load_amptk_table <- function(filepath, marker_type = "16S") {
  message("Reading OTU table: ", filepath)
  raw <- read.table(filepath, header = TRUE, sep = "\t",
                    comment.char = "", check.names = FALSE, quote = "")
  # first col = OTU ID, last col = Taxonomy, middle = sample counts
  otu_ids   <- as.character(raw[[1]])
  tax_col   <- as.character(raw[[ncol(raw)]])
  count_mat <- raw[, 2:(ncol(raw) - 1), drop = FALSE]

  rownames(count_mat) <- otu_ids
  names(tax_col)      <- otu_ids

  count_mat <- as.matrix(count_mat)
  storage.mode(count_mat) <- "integer"

  tax_df <- parse_amptk_taxonomy(tax_col)
  list(otu_mat = count_mat, tax_df = tax_df)
}

#' Load metadata, skipping the optional QIIME2 type-row
load_metadata <- function(filepath) {
  raw <- read.table(filepath, header = TRUE, sep = "\t",
                    comment.char = "", check.names = FALSE, quote = "")
  # QIIME2 puts a row of "categorical" / "numeric" as row 1 after header
  if (nrow(raw) > 0 &&
      any(vapply(raw[1, ], function(x) x %in% c("categorical","numeric","boolean"), logical(1)))) {
    raw <- raw[-1, , drop = FALSE]
  }
  id_col <- grep("^sample.?id$", colnames(raw), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(id_col)) rownames(raw) <- raw[[id_col]]
  # coerce obvious factor columns
  for (col in c("Country", "host", "material")) {
    if (col %in% colnames(raw)) raw[[col]] <- as.character(raw[[col]])
  }
  raw
}

#' Save both PDF and PNG versions of a ggplot
save_plot <- function(p, base_path, w = opt$width, h = opt$height, dpi = 150) {
  ggsave(paste0(base_path, ".pdf"), p, width = w, height = h)
  ggsave(paste0(base_path, ".png"), p, width = w, height = h, dpi = dpi)
  message("  Saved: ", basename(base_path))
}

# =============================================================================
# LOCATE DATA FILES
# =============================================================================
input_dir <- normalizePath(opt$input, mustWork = TRUE)
prefix    <- basename(input_dir)

# Determine table suffix (ASVs vs cluster)
if (!is.null(opt$table_type)) {
  tbl_suffix <- if (tolower(opt$table_type) == "asv") "ASVs" else "cluster"
} else {
  # auto-detect: prefer ASVs file if present
  if (file.exists(file.path(input_dir,
                             paste0(prefix, ".ASVs.otu_table.taxonomy.txt")))) {
    tbl_suffix <- "ASVs"
  } else if (file.exists(file.path(input_dir,
                                    paste0(prefix, ".cluster.otu_table.taxonomy.txt")))) {
    tbl_suffix <- "cluster"
  } else {
    tbl_suffix <- "ASVs"  # fallback; will error with a clear message
  }
  message("Auto-detected table suffix: ", tbl_suffix)
}

otu_file <- file.path(input_dir, paste0(prefix, ".", tbl_suffix, ".otu_table.taxonomy.txt"))
if (!file.exists(otu_file)) {
  candidates <- list.files(input_dir, pattern = "otu_table\\.taxonomy\\.txt$", full.names = TRUE)
  if (length(candidates) == 0)
    stop("Cannot find OTU/ASV table in: ", input_dir)
  otu_file <- candidates[1]
  message("Falling back to: ", otu_file)
}

tree_file <- file.path(input_dir, paste0(prefix, ".", tbl_suffix, ".tree.phy"))
if (!file.exists(tree_file)) tree_file <- NULL

# =============================================================================
# LOAD DATA & BUILD PHYLOSEQ
# =============================================================================
amptk_data <- load_amptk_table(otu_file, opt$type)
meta        <- load_metadata(opt$metadata)

# Align samples between OTU table and metadata
common_samples <- intersect(colnames(amptk_data$otu_mat), rownames(meta))
if (length(common_samples) == 0)
  stop("No samples overlap between OTU table and metadata. Check sample IDs.")

message(sprintf("Samples — OTU table: %d  |  metadata: %d  |  common: %d",
                ncol(amptk_data$otu_mat), nrow(meta), length(common_samples)))

otu_mat  <- amptk_data$otu_mat[, common_samples, drop = FALSE]
meta_sub <- meta[common_samples, , drop = FALSE]
tax_df   <- amptk_data$tax_df

# Remove low-depth samples
samp_sums <- colSums(otu_mat)
keep_samp <- names(samp_sums[samp_sums >= opt$min_reads])
dropped   <- setdiff(common_samples, keep_samp)
if (length(dropped) > 0)
  message(sprintf("Dropping %d samples with < %d reads: %s",
                  length(dropped), opt$min_reads, paste(dropped, collapse = ", ")))
otu_mat  <- otu_mat[, keep_samp, drop = FALSE]
meta_sub <- meta_sub[keep_samp, , drop = FALSE]

# Remove empty OTUs after sample filtering
otu_mat <- otu_mat[rowSums(otu_mat) > 0, , drop = FALSE]
tax_df  <- tax_df[rownames(otu_mat), , drop = FALSE]

# Replace NA strings in metadata
for (col in c("Country", "host", "material")) {
  if (col %in% colnames(meta_sub)) {
    meta_sub[[col]][is.na(meta_sub[[col]]) | meta_sub[[col]] == "NA"] <- "Unknown"
  }
}

if (!"description" %in% colnames(meta_sub)) {
  meta_sub$description <- rownames(meta_sub)
}
meta_sub$description <- trimws(as.character(meta_sub$description))
missing_description <- is.na(meta_sub$description) |
  meta_sub$description == "" |
  meta_sub$description == "NA"
meta_sub$description[missing_description] <- rownames(meta_sub)[missing_description]

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_df)),
  sample_data(meta_sub)
)

# Load tree if present
if (!is.null(tree_file) && file.exists(tree_file)) {
  tryCatch({
    tree <- read_tree(tree_file)
    ps   <- merge_phyloseq(ps, tree)
    message("Phylogenetic tree loaded.")
  }, error = function(e) message("Warning: tree load failed – ", e$message))
}

message(sprintf("Phyloseq: %d taxa x %d samples", ntaxa(ps), nsamples(ps)))

build_sample_label_map <- function(ps_obj) {
  sdat <- as.data.frame(sample_data(ps_obj))
  sample_ids <- rownames(sdat)

  desc <- if ("description" %in% colnames(sdat)) {
    trimws(as.character(sdat$description))
  } else {
    sample_ids
  }
  missing_desc <- is.na(desc) | desc == "" | desc == "NA"
  desc[missing_desc] <- sample_ids[missing_desc]

  species <- if ("host" %in% colnames(sdat)) {
    trimws(as.character(sdat$host))
  } else {
    rep(NA_character_, length(sample_ids))
  }

  has_species <- !is.na(species) & species != "" & species != "NA" & species != "Unknown"
  sample_labels <- ifelse(has_species, paste(species, desc, sep = " | "), desc)
  stats::setNames(sample_labels, sample_ids)
}

apply_sample_labels <- function(plot_obj, ps_obj) {
  plot_obj + scale_x_discrete(labels = build_sample_label_map(ps_obj))
}

# =============================================================================
# COLOUR / SHAPE PALETTES
# =============================================================================
hosts     <- sort(unique(na.omit(sample_data(ps)$host)))
countries <- sort(unique(na.omit(sample_data(ps)$Country)))
materials <- sort(unique(na.omit(sample_data(ps)$material)))

n_hosts     <- length(hosts)
n_countries <- length(countries)

# Host colour palette
if (n_hosts <= 8) {
  host_pal <- setNames(brewer.pal(max(3, n_hosts), "Set2")[seq_len(n_hosts)], hosts)
} else if (n_hosts <= 12) {
  host_pal <- setNames(brewer.pal(12, "Paired")[seq_len(n_hosts)], hosts)
} else {
  host_pal <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(n_hosts), hosts)
}

# Country shapes: use shapes 21-25 (support independent fill)
base_shapes  <- c(21L, 22L, 23L, 24L, 25L)
extra_shapes <- c(3L, 4L, 7L, 8L, 9L, 10L)
country_shapes <- setNames(
  c(base_shapes, extra_shapes)[seq_len(n_countries)],
  countries
)

# Material fill label (for legend)
material_labels <- c(honey = "filled", pollen = "open", unknown = "half", Unknown = "half")

# =============================================================================
# RAREFACTION / NORMALISATION
# =============================================================================
rarefy_depth <- if (is.null(opt$rarefy)) {
  d <- min(sample_sums(ps))
  message("Auto rarefaction depth: ", format(d, big.mark = ","))
  d
} else {
  opt$rarefy
}

set.seed(42)
ps_rare <- suppressMessages(
  rarefy_even_depth(ps, sample.size = rarefy_depth, rngseed = 42,
                    replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
)
message(sprintf("After rarefaction: %d taxa", ntaxa(ps_rare)))

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# =============================================================================
# SECTION 1: SAMPLE READ DEPTH
# =============================================================================
message("\n--- Sample read depths ---")

depth_df <- as.data.frame(sample_data(ps))
depth_df$Reads  <- sample_sums(ps)
depth_df$Sample <- rownames(depth_df)
depth_df <- depth_df[order(depth_df$material, depth_df$host), ]
depth_df$Sample <- factor(depth_df$Sample, levels = depth_df$Sample)

p_depth <- ggplot(depth_df, aes(x = Sample, y = Reads, fill = host)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = rarefy_depth, linetype = "dashed", colour = "red", linewidth = 0.8) +
  scale_fill_manual(values = host_pal, name = "Host") +
  scale_y_continuous(labels = comma, expand = c(0, 0)) +
  facet_grid(~material, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        strip.background = element_rect(fill = "grey88")) +
  labs(title    = paste(opt$type, "Sample Read Depths"),
       subtitle = paste0("Red dashed line = rarefaction depth (",
                         format(rarefy_depth, big.mark = ","), ")"),
       x = NULL, y = "Read Count")

p_depth <- apply_sample_labels(p_depth, ps)

save_plot(p_depth,
          file.path(opt$outdir, paste0(opt$type, "_sample_depth")),
          h = 6)

# =============================================================================
# SECTION 2: RAREFACTION CURVES
# =============================================================================
message("\n--- Rarefaction curves ---")

rare_tidy <- tryCatch({
  rarecurve(as.matrix(t(otu_table(ps))),
            step = max(200, round(rarefy_depth / 50)),
            sample = rarefy_depth, tidy = TRUE)
}, error = function(e) { message("Rarefaction error: ", e$message); NULL })

if (!is.null(rare_tidy)) {
  sdat <- as.data.frame(sample_data(ps))
  rare_tidy$host     <- sdat[rare_tidy$Site, "host"]
  rare_tidy$material <- sdat[rare_tidy$Site, "material"]
  rare_tidy$Country  <- sdat[rare_tidy$Site, "Country"]

  p_rare <- ggplot(rare_tidy,
                   aes(x = Sample, y = Species, group = Site,
                       colour = host, linetype = material)) +
    geom_line(alpha = 0.65, linewidth = 0.6) +
    geom_vline(xintercept = rarefy_depth, linetype = "dashed",
               colour = "red", alpha = 0.7) +
    scale_colour_manual(values = host_pal, name = "Host") +
    scale_linetype_manual(
      values = c(honey = "solid", pollen = "dashed",
                 unknown = "dotted", Unknown = "dotted"),
      name = "Material") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right") +
    labs(title = paste(opt$type, "Rarefaction Curves"),
         x = "Reads Sampled", y = "Taxa Observed")

  save_plot(p_rare,
            file.path(opt$outdir, paste0(opt$type, "_rarefaction")),
            h = opt$height - 1)
}

# =============================================================================
# SECTION 3: ALPHA DIVERSITY
# =============================================================================
message("\n--- Alpha diversity ---")

alpha_div <- estimate_richness(ps_rare,
                               measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_div$Sample <- rownames(alpha_div)

sdat_rare <- as.data.frame(sample_data(ps_rare))
alpha_div  <- merge(alpha_div, sdat_rare,
                    by.x = "Sample", by.y = "row.names", all.x = TRUE)

# Order hosts by median Shannon
host_order <- alpha_div |>
  dplyr::group_by(host) |>
  dplyr::summarise(med = median(Shannon, na.rm = TRUE)) |>
  dplyr::arrange(dplyr::desc(med)) |>
  dplyr::pull(host)

alpha_long <- tidyr::pivot_longer(alpha_div,
                                   cols = c("Observed", "Chao1", "Shannon", "Simpson"),
                                   names_to = "metric", values_to = "value")
alpha_long$host <- factor(alpha_long$host, levels = host_order)

p_alpha <- ggplot(alpha_long,
                  aes(x = host, y = value, fill = host)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(aes(shape = material, fill = host),
              width = 0.2, size = 2.5, colour = "grey25", alpha = 0.85) +
  scale_fill_manual(values = host_pal, name = "Host") +
  scale_shape_manual(
    values = c(honey = 21, pollen = 24, unknown = 23, Unknown = 23),
    name = "Material") +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        strip.background = element_rect(fill = "grey88")) +
  labs(title = paste(opt$type, "Alpha Diversity"), x = NULL, y = "Value")

save_plot(p_alpha,
          file.path(opt$outdir, paste0(opt$type, "_alpha_diversity")))

# Shannon by material
p_alpha_mat <- ggplot(alpha_div,
                      aes(x = material, y = Shannon, fill = Country)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(colour = Country), width = 0.2, size = 2.2, alpha = 0.8) +
  scale_fill_brewer(palette = "Set1")   +
  scale_colour_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(title = paste(opt$type, "Shannon Diversity by Material"),
       x = "Material", y = "Shannon Index")

save_plot(p_alpha_mat,
          file.path(opt$outdir, paste0(opt$type, "_alpha_by_material")),
          w = 8, h = 6)

# =============================================================================
# SECTION 4: TAXONOMY BARPLOTS (standard + microshade)
# =============================================================================
message("\n--- Taxonomy barplots ---")

#' Generate a taxon label used in barplots at a given rank
label_taxa <- function(tax_table_df, rank) {
  lbl <- tax_table_df[[rank]]
  lbl[is.na(lbl) | lbl == "NA"] <- paste0("Unclassified")
  lbl
}

#' Standard stacked barplot at one taxonomic rank
plot_taxa_bar <- function(ps_obj, rank, top_n = 15, marker_type = "16S") {
  # Agglomerate
  ps_g <- tryCatch(tax_glom(ps_obj, taxrank = rank, NArm = FALSE),
                   error = function(e) ps_obj)

  otu_sums <- taxa_sums(ps_g)
  top_names <- names(sort(otu_sums, decreasing = TRUE))[seq_len(min(top_n, length(otu_sums)))]

  melt_df <- psmelt(ps_g)

  # Labels
  melt_df$TaxLabel <- as.character(melt_df[[rank]])
  melt_df$TaxLabel[is.na(melt_df$TaxLabel) | melt_df$TaxLabel == "NA"] <- "Unclassified"
  melt_df$TaxLabel[!melt_df$OTU %in% top_names] <- "Other"

  # Sample order: material > host > Country
  samp_ord <- melt_df |>
    dplyr::select(Sample, material, host, Country) |>
    dplyr::distinct() |>
    dplyr::arrange(material, host, Country) |>
    dplyr::pull(Sample)
  melt_df$Sample <- factor(melt_df$Sample, levels = samp_ord)

  agg_df <- melt_df |>
    dplyr::group_by(Sample, TaxLabel, material, host, Country) |>
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

  special <- c("Other", "Unclassified")
  tax_levels <- c(setdiff(unique(agg_df$TaxLabel), special), special)
  agg_df$TaxLabel <- factor(agg_df$TaxLabel, levels = rev(tax_levels))

  n_taxa <- length(tax_levels) - length(intersect(tax_levels, special))
  tax_pal <- c(
    if (n_taxa > 0) {
      if (n_taxa <= 12) brewer.pal(max(3, n_taxa), "Paired")[seq_len(n_taxa)]
      else colorRampPalette(brewer.pal(12, "Paired"))(n_taxa)
    },
    "grey65", "grey88"
  )
  names(tax_pal) <- c(
    setdiff(unique(as.character(agg_df$TaxLabel)), special),
    special
  )

  p <- ggplot(agg_df, aes(x = Sample, y = Abundance, fill = TaxLabel)) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = tax_pal, name = rank) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    facet_grid(~material, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      legend.text      = element_text(size = 9),
      legend.key.size  = unit(0.4, "cm"),
      strip.background = element_rect(fill = "grey88"),
      panel.spacing    = unit(0.15, "lines")
    ) +
    labs(title = paste(marker_type, rank, "Composition"),
         x = NULL, y = "Relative Abundance")

  apply_sample_labels(p, ps_obj)
}

for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  tryCatch({
    p_bar <- plot_taxa_bar(ps_rel, rank, top_n = opt$top_taxa, marker_type = opt$type)
    save_plot(p_bar,
              file.path(opt$outdir, paste0(opt$type, "_barplot_", rank)),
              w = opt$width + 2, h = opt$height)
  }, error = function(e) message("  Barplot ", rank, " skipped: ", e$message))
}

# ---- microshade barplot -----------------------------------------------------
# microshade creates shaded sub-colours within each top-level group (Phylum)
# so related taxa at Genus level share a hue family

if (HAS_MICROSHADE) {
  message("\n--- microshade barplots ---")
  tryCatch({
    # microshade_prep returns a list with long-format data and colour palette
    ms_prep <- microshade_prep(ps_rel,
                               top_group = "Phylum",
                               sub_group = "Genus",
                               top_n     = 5)

    # Sample order
    sdat_ord <- as.data.frame(sample_data(ps_rel))
    samp_ord  <- rownames(sdat_ord)[order(sdat_ord$material, sdat_ord$host)]

    p_ms <- microshade_plot(ms_prep,
                            sample_order = samp_ord,
                            legend_text_size = 7) +
      facet_grid(~ sdat_ord[as.character(.) , "material"],
                 scales = "free_x", space = "free_x") +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        strip.background = element_rect(fill = "grey88")
      ) +
      labs(title = paste(opt$type, "Phylum > Genus Composition (microshade)"))

    p_ms <- apply_sample_labels(p_ms, ps_rel)

    save_plot(p_ms,
              file.path(opt$outdir, paste0(opt$type, "_microshade_phylum_genus")),
              w = opt$width + 2, h = opt$height)
  }, error = function(e) message("  microshade plot failed: ", e$message))
} else {
  # Fallback: manual shaded colour scheme using HSV within each phylum group
  message("\n--- Shaded barplot (microshade fallback) ---")
  tryCatch({
    ps_g <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)

    # Top 5 phyla by abundance
    ps_p <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
    top5_phyla <- names(sort(taxa_sums(ps_p), decreasing = TRUE))[1:5]
    top5_names <- as.character(tax_table(ps_p)[top5_phyla, "Phylum"])

    melt_g <- psmelt(ps_g)
    melt_g$Phylum[is.na(melt_g$Phylum) | melt_g$Phylum == "NA"] <- "Unclassified"
    melt_g$Genus [is.na(melt_g$Genus)  | melt_g$Genus  == "NA"] <- "Unclassified_Genus"

    # Keep only genera in top phyla; collapse rest to "Other"
    melt_g$GLabel <- ifelse(melt_g$Phylum %in% top5_names,
                            paste(melt_g$Phylum, melt_g$Genus, sep = " | "),
                            "Other")

    # Build shaded palette: each phylum gets a base hue; genera get lighter shades
    phylum_hues <- c(0.12, 0.55, 0.35, 0.72, 0.04)
    names(phylum_hues) <- top5_names

    genus_pal <- c()
    for (ph in top5_names) {
      genera_in <- unique(melt_g$Genus[melt_g$Phylum == ph & melt_g$GLabel != "Other"])
      n_g <- length(genera_in)
      if (n_g == 0) next
      shades <- hsv(phylum_hues[ph],
                    seq(1.0, 0.25, length.out = max(n_g, 1)),
                    seq(0.5, 0.95, length.out = max(n_g, 1)))
      lbl <- paste(ph, genera_in, sep = " | ")
      genus_pal <- c(genus_pal, setNames(shades, lbl))
    }
    genus_pal["Other"] <- "grey80"

    # Sample order
    sdat_ms <- as.data.frame(sample_data(ps_rel))
    samp_ord <- rownames(sdat_ms)[order(sdat_ms$material, sdat_ms$host)]
    melt_g$Sample <- factor(melt_g$Sample, levels = samp_ord)

    agg_ms <- melt_g |>
      dplyr::group_by(Sample, GLabel) |>
      dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
    agg_ms <- merge(agg_ms,
                    melt_g[, c("Sample", "material"), drop = FALSE] |> dplyr::distinct(),
                    by = "Sample")
    agg_ms$Sample <- factor(agg_ms$Sample, levels = samp_ord)

    glabels_ord <- c(names(genus_pal)[names(genus_pal) != "Other"], "Other")
    agg_ms$GLabel <- factor(agg_ms$GLabel, levels = rev(glabels_ord))

    p_ms_fallback <- ggplot(agg_ms, aes(x = Sample, y = Abundance, fill = GLabel)) +
      geom_bar(stat = "identity", width = 0.85) +
      scale_fill_manual(values = genus_pal, name = "Phylum | Genus") +
      scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
      facet_grid(~material, scales = "free_x", space = "free_x") +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.4, "cm"),
        strip.background = element_rect(fill = "grey88")
      ) +
      labs(title = paste(opt$type, "Phylum > Genus (shaded palette)"),
           x = NULL, y = "Relative Abundance")

    p_ms_fallback <- apply_sample_labels(p_ms_fallback, ps_rel)

    save_plot(p_ms_fallback,
              file.path(opt$outdir, paste0(opt$type, "_shaded_phylum_genus")),
              w = opt$width + 2, h = opt$height)
  }, error = function(e) message("  Shaded barplot failed: ", e$message))
}

# =============================================================================
# SECTION 5: BETA DIVERSITY
# =============================================================================
message("\n--- Beta diversity ---")

#' Build PCoA or NMDS beta diversity plot
#' Aesthetics: colour = host, shape = country, fill = host (honey) / white (pollen)
plot_beta_div <- function(ps_obj, dist_method = "bray", ord_type = "PCoA",
                           marker_type = "16S") {
  # Distance matrix
  dist_mat <- tryCatch(
    phyloseq::distance(ps_obj, method = dist_method),
    error = function(e) stop("Distance failed (", dist_method, "): ", e$message)
  )

  # Ordination
  ord <- ordinate(ps_obj, method = ord_type, distance = dist_mat)

  # Extract coordinates + axis labels
  if (ord_type == "PCoA") {
    coords <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
    colnames(coords) <- c("Axis1", "Axis2")
    pct <- round(ord$values$Relative_eig[1:2] * 100, 1)
    ax_lab <- c(sprintf("PCoA1 (%.1f%%)", pct[1]), sprintf("PCoA2 (%.1f%%)", pct[2]))
  } else {  # NMDS
    sc <- tryCatch(scores(ord, display = "sites"), error = function(e) ord$points)
    coords  <- as.data.frame(sc[, 1:2, drop = FALSE])
    colnames(coords) <- c("Axis1", "Axis2")
    stress  <- tryCatch(round(ord$stress, 4), error = function(e) NA)
    ax_lab  <- c(sprintf("NMDS1 (stress=%.4f)", stress), "NMDS2")
  }

  # Merge with metadata
  sdat <- as.data.frame(sample_data(ps_obj))
  pdat <- cbind(coords, sdat[rownames(coords), , drop = FALSE])

  # PERMANOVA  (suppress error if model is rank-deficient)
  perm <- tryCatch({
    adonis2(dist_mat ~ host + Country + material,
            data = sdat[labels(dist_mat), , drop = FALSE],
            permutations = 999, by = "margin")
  }, error = function(e) NULL)

  subtitle <- if (!is.null(perm)) {
    r2 <- round(perm[c("host","Country","material"), "R2"], 3)
    paste0("PERMANOVA (margin): host R²=", r2["host"],
           ", country R²=", r2["Country"],
           ", material R²=", r2["material"])
  } else ""

  # Key aesthetic:
  #  • shape  = Country (21-25 supports independent fill)
  #  • colour = host (point border colour)
  #  • fill   = host for honey; "white" for pollen/unknown
  pdat$pt_fill <- ifelse(pdat$material == "honey",
                         host_pal[pdat$host], "white")

  ggplot(pdat, aes(x = Axis1, y = Axis2)) +
    # confidence ellipses per host
    stat_ellipse(aes(colour = host, group = host),
                 type = "norm", linetype = 2, level = 0.75, alpha = 0.45) +
    # all points with shape = country, stroke = host colour
    geom_point(aes(shape = Country, colour = host, fill = pt_fill),
               size = 3.5, stroke = 1.3, alpha = 0.88) +
    scale_shape_manual(values = country_shapes, name = "Country") +
    scale_colour_manual(values = host_pal, name = "Host") +
    scale_fill_identity() +   # fill is pre-computed
    # manual legend for material (open/filled)
    annotate("text",
             x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = "Filled circle = honey   |   Open/white = pollen",
             size = 2.8, colour = "grey45") +
    theme_bw(base_size = 12) +
    theme(legend.position  = "right",
          legend.key.size  = unit(0.5, "cm"),
          legend.text      = element_text(size = 9)) +
    labs(title    = paste(marker_type, ord_type, toupper(dist_method)),
         subtitle = subtitle,
         x = ax_lab[1], y = ax_lab[2])
}

# Define which methods to run (UniFrac added only if tree present)
beta_cfg <- list(
  list(method = "bray",    type = "PCoA", tag = "bray_pcoa"),
  list(method = "jaccard", type = "PCoA", tag = "jaccard_pcoa"),
  list(method = "bray",    type = "NMDS", tag = "bray_nmds")
)
if (!is.null(phy_tree(ps_rare, errorIfNULL = FALSE))) {
  beta_cfg <- c(beta_cfg, list(
    list(method = "wunifrac", type = "PCoA", tag = "wunifrac_pcoa"),
    list(method = "unifrac",  type = "PCoA", tag = "unifrac_pcoa")
  ))
}

beta_plots <- list()
for (cfg in beta_cfg) {
  message("  ", cfg$tag)
  tryCatch({
    p <- plot_beta_div(ps_rare, cfg$method, cfg$type, opt$type)
    beta_plots[[cfg$tag]] <- p
    save_plot(p,
              file.path(opt$outdir, paste0(opt$type, "_beta_", cfg$tag)),
              h = opt$height - 1)
  }, error = function(e) message("  SKIPPED (", cfg$tag, "): ", e$message))
}

# Combined 2-panel (Bray PCoA + NMDS)
if (length(beta_plots) >= 2) {
  tryCatch({
    p_comb <- wrap_plots(beta_plots[1:min(2L, length(beta_plots))], ncol = 2) +
      plot_layout(guides = "collect") &
      theme(legend.position = "right")
    ggsave(file.path(opt$outdir, paste0(opt$type, "_beta_combined.pdf")),
           p_comb, width = opt$width * 1.8, height = opt$height)
    message("  Saved: ", opt$type, "_beta_combined.pdf")
  }, error = function(e) message("  Combined panel failed: ", e$message))
}

# =============================================================================
# SECTION 6: TOP-GENUS HEATMAP
# =============================================================================
message("\n--- Genus heatmap ---")

tryCatch({
  n_heat <- min(opt$top_taxa * 2, 40)
  ps_gen <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
  top_g  <- names(sort(taxa_sums(ps_gen), decreasing = TRUE))[seq_len(min(n_heat, ntaxa(ps_gen)))]
  ps_top <- prune_taxa(top_g, ps_gen)

  heat_df <- as.data.frame(t(otu_table(ps_top)))
  ttbl    <- as.data.frame(tax_table(ps_top))
  colnames(heat_df) <- ifelse(is.na(ttbl$Genus) | ttbl$Genus == "NA",
                              paste0("Unclassified_", seq_len(ncol(heat_df))),
                              ttbl$Genus)

  sdat_h <- as.data.frame(sample_data(ps_top))
  heat_df$material <- sdat_h[rownames(heat_df), "material"]
  heat_df$host     <- sdat_h[rownames(heat_df), "host"]

  samp_ord_h <- rownames(sdat_h)[order(sdat_h$material, sdat_h$host)]

  heat_long <- heat_df |>
    tibble::rownames_to_column("Sample") |>
    tidyr::pivot_longer(cols = -c(Sample, material, host),
                        names_to = "Genus", values_to = "RelAbund")
  heat_long$Sample <- factor(heat_long$Sample, levels = samp_ord_h)

  p_heat <- ggplot(heat_long, aes(x = Sample, y = Genus, fill = RelAbund)) +
    geom_tile(colour = NA) +
    scale_fill_gradientn(
      colours = c("white", "#deebf7", "#9ecae1", "#3182bd", "#08306b"),
      labels  = percent_format(accuracy = 0.1),
      name    = "Relative\nAbundance",
      na.value = "white"
    ) +
    facet_grid(~material, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 10) +
    theme(axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y      = element_text(size = 8),
          strip.background = element_rect(fill = "grey88"),
          panel.spacing    = unit(0.1, "lines")) +
    labs(title = paste(opt$type, "Top", n_heat, "Genus Abundance Heatmap"),
         x = NULL, y = NULL)

  p_heat <- apply_sample_labels(p_heat, ps_top)

  save_plot(p_heat,
            file.path(opt$outdir, paste0(opt$type, "_heatmap_genus")),
            w = opt$width + 2, h = max(8, n_heat * 0.3 + 2))
}, error = function(e) message("  Heatmap failed: ", e$message))

# =============================================================================
# SECTION 7: HONEY vs POLLEN MA-PLOT
# =============================================================================
message("\n--- Honey vs Pollen comparison ---")

tryCatch({
  ps_gen2 <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  keep_m  <- sample_data(ps_gen2)$material %in% c("honey", "pollen")
  ps_mat  <- prune_samples(keep_m, ps_gen2)
  ps_mr   <- transform_sample_counts(ps_mat, function(x) x / sum(x))

  ma_df <- psmelt(ps_mr) |>
    dplyr::group_by(Genus, material) |>
    dplyr::summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = material, values_from = mean_abund,
                       values_fill = 0) |>
    dplyr::filter(!is.na(Genus), Genus != "NA") |>
    dplyr::mutate(
      log2FC    = log2((honey + 1e-6) / (pollen + 1e-6)),
      avg_abund = (honey + pollen) / 2,
      enriched  = dplyr::case_when(
        log2FC >  1 ~ "Honey",
        log2FC < -1 ~ "Pollen",
        TRUE        ~ "Neither"
      )
    ) |>
    dplyr::arrange(dplyr::desc(abs(log2FC)))

  top_ma <- head(ma_df, 20)

  p_ma <- ggplot(ma_df, aes(x = log2(avg_abund + 1e-6), y = log2FC)) +
    geom_point(data = subset(ma_df, enriched == "Neither"),
               colour = "grey70", size = 1.5, alpha = 0.5) +
    geom_point(data = subset(ma_df, enriched != "Neither"),
               aes(colour = enriched), size = 2.5, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = top_ma |> dplyr::filter(enriched != "Neither"),
      aes(label = Genus, colour = enriched),
      size = 3, max.overlaps = 20, show.legend = FALSE
    ) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, linetype = "solid",  colour = "black", linewidth = 0.4) +
    scale_colour_manual(values = c(Honey = "#d73027", Pollen = "#4575b4", Neither = "grey70"),
                        name = "Enriched in") +
    theme_bw(base_size = 12) +
    labs(title    = paste(opt$type, "Genus: Honey vs Pollen (MA-plot)"),
         subtitle = "Dashed lines: |log2FC| = 1",
         x = "log2(Mean Relative Abundance + 1e-6)",
         y = "log2(Honey / Pollen)")

  save_plot(p_ma,
            file.path(opt$outdir, paste0(opt$type, "_honey_vs_pollen_MA")),
            w = 10, h = 8)
}, error = function(e) message("  MA-plot failed (ggrepel available?): ", e$message))

# =============================================================================
# SAVE PHYLOSEQ OBJECTS
# =============================================================================
message("\n--- Saving phyloseq objects ---")
saveRDS(ps,      file.path(opt$outdir, paste0(opt$type, "_phyloseq.rds")))
saveRDS(ps_rare, file.path(opt$outdir, paste0(opt$type, "_phyloseq_rarefied.rds")))

message("\n========================================")
message("DONE. Output written to: ", opt$outdir)
message("========================================\n")
