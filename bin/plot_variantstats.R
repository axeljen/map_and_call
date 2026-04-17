#!/usr/bin/env Rscript

# Variant Statistics Report Generator
# This script creates a comprehensive variant statistics PDF report from combine_stats TSV files
# Usage: Rscript plot_variantstats.R <input_directory>
#   input_directory: Directory containing *_ab_dp.tsv, *_qual_fmiss_maf_dp.tsv, 
#                    *_record_counts.tsv, and *_sample_sumstats.tsv files (outputs from combine_stats.py)

suppressPackageStartupMessages({
  library(tidyverse)
  library(knitr)
  library(rmarkdown)
  library(scales)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_variantstats.R <input_directory>\n  Reads combine_stats.py output files and generates PDF report")
}

input_dir <- args[1]
# Convert to absolute path
input_dir <- normalizePath(input_dir, mustWork = TRUE)
output_format <- "html"  # Always output as HTML for now

if (!dir.exists(input_dir)) {
  stop(paste("Input directory does not exist:", input_dir))
}

# Find input files from combine_stats.py output
ab_file <- list.files(input_dir, pattern = "_ab\\.tsv$", full.names = TRUE)
dp_file <- list.files(input_dir, pattern = "_dp\\.tsv$", full.names = TRUE)  # Use same file for DP if available
qual_file <- list.files(input_dir, pattern = "_qual_fmiss_maf_dp\\.tsv$", full.names = TRUE)
counts_file <- list.files(input_dir, pattern = "_record_counts\\.tsv$", full.names = TRUE)
sample_file <- list.files(input_dir, pattern = "_sample_stats\\.tsv$", full.names = TRUE)

# Extract dataset name from directory or file names
dir_basename <- basename(normalizePath(input_dir))
if (dir_basename == "." || dir_basename == "") {
  dataset_name <- "Variant Statistics"
  file_prefix <- "variant_stats"
} else {
  dataset_name <- gsub("_", " ", dir_basename)
  dataset_name <- tools::toTitleCase(dataset_name)
  file_prefix <- dir_basename
}

# Create RMarkdown content
rmd_content <- sprintf('---
title: "%s Report"
date: "`r format(Sys.time(), \'%%B %%d, %%Y\')`"
output:
  %s:
    toc: true
    toc_depth: 3
    number_sections: true
    %s
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  fig.width = 10,
  fig.height = 6,
  dpi = 300
)

library(tidyverse)
library(scales)
library(knitr)

theme_set(theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  ))
```

', 
dataset_name,
ifelse(output_format == "html", "html_document", "pdf_document"),
ifelse(output_format == "html", "theme: cosmo\n    code_folding: hide", ""))

# Add Summary section
if (length(sample_file) > 0) {
  rmd_content <- paste0(rmd_content, sprintf('
# Summary

```{r summary-stats}
sample_data <- read_tsv("%s", show_col_types = FALSE)
n_samples <- nrow(sample_data)
# Get total sites from record_counts file
sample_file_path <- "%s"
record_counts_files <- list.files(dirname(sample_file_path), pattern = "_record_counts\\\\.tsv$", full.names = TRUE)
if (length(record_counts_files) > 0) {
  record_counts <- read_tsv(record_counts_files[1], show_col_types = FALSE)
  total_sites <- sum(record_counts$count)
} else {
  # Fallback: use max num_records if record_counts not available
  total_sites <- max(sample_data$num_records)
}
cat(sprintf("**Number of samples:** %%s\\n\\n", comma(n_samples)))
cat(sprintf("**Sites in VCF:** %%s\\n\\n", comma(total_sites)))
```

---

', sample_file[1], sample_file[1]))
} else {
  rmd_content <- paste0(rmd_content, '
# Summary

No sample statistics available.

---

')
}

# Add Record Counts section
if (length(counts_file) > 0) {
  rmd_content <- paste0(rmd_content, sprintf('
# Variant Type Distribution

```{r load-counts}
counts <- read_tsv("%s", show_col_types = FALSE)
```

## Genotypes Summary

```{r counts-plot, fig.height=5}
counts %%>%%
  filter(count > 0) %%>%%
  mutate(
    record_type = factor(record_type, levels = record_type[order(count)]),
    label_text = paste0(comma(count), "\\n(", percent(count/sum(count), accuracy = 0.1), ")")
  ) %%>%%
  ggplot(aes(x = record_type, y = count, fill = record_type)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = label_text), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribution of Variant Record Types",
    x = NULL,
    y = "Number of Records"
  )
```

', counts_file[1]))
}

# Add Sample Statistics section
if (length(sample_file) > 0) {
  rmd_content <- paste0(rmd_content, sprintf('
# Per-Sample Genotype Statistics

```r
```{r load-samples}
sample_stats <- read_tsv("%s", show_col_types = FALSE) %%>%%
  mutate(across(everything(), ~replace_na(., 0)))
```

## Genotype Distribution by Sample

```{r genotype-dist, fig.height=max(4, nrow(sample_stats) * 0.4)}
sample_stats %%>%%
  pivot_longer(
    cols = all_of(c("num_hom_ref", "num_het", "num_hom_alt", "num_missing")),
    names_to = "genotype",
    values_to = "count"
  ) %%>%%
  mutate(
    genotype = case_when(
      genotype == "num_hom_ref" ~ "Homozygous Reference",
      genotype == "num_het" ~ "Heterozygous",
      genotype == "num_hom_alt" ~ "Homozygous Alternate",
      genotype == "num_missing" ~ "Missing",
      genotype == "num_records" ~ "Total Genotypes",
      TRUE ~ NA_character_
    )
  ) %%>%%
  filter(genotype != "Total Genotypes") %%>%%
  mutate(
    genotype = fct_reorder(genotype, rev(count), .desc = TRUE),
    sample = fct_rev(factor(sample))
  ) %%>%%
  group_by(sample) %%>%%
  mutate(percentage = count / sum(count) * 100) %%>%%
  ungroup() %%>%%
  ggplot(aes(x = sample, y = percentage, fill = genotype)) +
  geom_col(position = "stack") +
  geom_text(
    aes(label = ifelse(percentage > 5, paste0(round(percentage, 1), "%%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Homozygous Reference" = "#2166ac",
      "Heterozygous" = "#4daf4a",
      "Homozygous Alternate" = "#d62728",
      "Missing" = "#999999"
    )
  ) +
  labs(
    title = "Genotype Distribution Across Samples",
    x = NULL,
    y = "Percentage (%%)",
    fill = "Genotype"
  )
```

## Missing Genotypes

```{r missing-plot, fig.height=max(4, nrow(sample_stats) * 0.3)}
sample_stats %%>%%
  mutate(
    missing_rate = num_missing / num_records * 100,
    sample = fct_reorder(factor(sample), missing_rate)
  ) %%>%%
  ggplot(aes(x = sample, y = missing_rate)) +
  geom_col(fill = "#999999") +
  geom_text(aes(label = sprintf("%%.1f%%%%", missing_rate)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) paste0(x, "%%%%"),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Missing Genotype Rate by Sample",
    x = NULL,
    y = "Missing Rate (%%%%)"
  )
```

', sample_file[1]))
}

# Add Allele Balance section
if (length(ab_file) > 0) {
  rmd_content <- paste0(rmd_content, sprintf('
# Allele Balance

```{r load-ab}
ab_data <- read_tsv("%s", show_col_types = FALSE)
```

## Allele Balance by Sample

```{r ab-by-sample, fig.height=max(6, length(unique(ab_data$sample)) * 0.5)}
ab_data %%>%%
  ggplot(aes(x = minor_allele_support, y = sample)) +
  geom_violin(fill = "#984ea3", alpha = 0.6) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(
    title = "Allele Balance Distribution by Sample",
    x = "Minor Allele Support",
    y = NULL
  )
```

', ab_file[1]))
}

# Add Quality/Missingness/MAF section
if (length(qual_file) > 0) {
  rmd_content <- paste0(rmd_content, sprintf('
# Variant Quality Metrics

```{r load-qual}
qual_data <- read_tsv("%s", show_col_types = FALSE)
```

## Quality Distribution

```{r qual-hist}
ggplot(qual_data, aes(x = qual)) +
  geom_histogram(bins = 50, fill = "#377eb8", color = "white") +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Distribution of Variant Quality Scores",
    x = "Quality (QUAL, log10 scale)",
    y = "Number of Variants"
  )
```

## Minor Allele Frequency Distribution

```{r maf-hist}
qual_data %%>%%
  filter(maf > 0) %%>%%
  ggplot(aes(x = maf)) +
  geom_histogram(bins = 50, fill = "#4daf4a", color = "white") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Minor Allele Frequency (MAF) Distribution",
    subtitle = "Sites with MAF > 0",
    x = "Minor Allele Frequency",
    y = "Number of Variants"
  )
```

', qual_file[1]))
}

# Add Depth section (if file is not too large)
if (length(dp_file) > 0) {
  file_size <- file.info(dp_file[1])$size / (1024^2)  # Size in MB
  
  if (file_size < 100) {  # Only process if less than 100 MB
    rmd_content <- paste0(rmd_content, sprintf('
# Genotype Depth (FORMAT/DP)

```{r load-dp}
# Load all depth data
dp_data <- read_tsv("%s", show_col_types = FALSE, n_max = 1000000)
```

## Site-level Depth Distribution (All Samples)

```{r dp-site-level-hist}
# Site-level depth from qual_fmiss_maf_dp file (if dp column exists)
if ("dp" %%in%% names(qual_data)) {
  dp_99_site <- quantile(qual_data$dp, 0.99, na.rm = TRUE)
  
  qual_data %%>%%
    filter(dp <= dp_99_site) %%>%%
    ggplot(aes(x = dp)) +
    geom_histogram(bins = 50, fill = "#ff7f00", color = "white") +
    geom_vline(aes(xintercept = median(dp, na.rm = TRUE)),
               linetype = "dashed", color = "red", size = 1) +
    annotate("text", 
             x = median(qual_data$dp, na.rm = TRUE),
             y = Inf,
             label = sprintf("Median = %%.1f", median(qual_data$dp, na.rm = TRUE)),
             hjust = -0.1, vjust = 2, color = "red", fontface = "bold") +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Site-level Depth Distribution Across All Samples",
      subtitle = sprintf("Truncated at 99th percentile (%%.1f)", dp_99_site),
      x = "Depth (DP)",
      y = "Number of Variant Sites"
    )
} else {
  cat("*Note: Site-level depth data not available in this dataset.*\\n")
}
```

## Depth Summary Statistics

```{r dp-summary-table}
dp_summary <- dp_data %%>%%
  group_by(sample) %%>%%
  summarise(
    Mean = mean(genotype_depth, na.rm = TRUE),
    Median = median(genotype_depth, na.rm = TRUE),
    SD = sd(genotype_depth, na.rm = TRUE),
    Min = min(genotype_depth, na.rm = TRUE),
    Q25 = quantile(genotype_depth, 0.25, na.rm = TRUE),
    Q75 = quantile(genotype_depth, 0.75, na.rm = TRUE),
    Max = max(genotype_depth, na.rm = TRUE),
    N_Genotypes = n()
  ) %%>%%
  arrange(desc(Median))

kable(
  dp_summary,
  digits = 2,
  format.args = list(big.mark = ","),
  caption = "Summary Statistics of Genotype Depth by Sample"
)
```

## Depth Distribution by Sample

```{r dp-by-sample, fig.height=max(6, length(unique(dp_data$sample)) * 0.5)}
dp_data %%>%%
  filter(genotype_depth <= quantile(genotype_depth, 0.99, na.rm = TRUE)) %%>%%
  ggplot(aes(x = genotype_depth, y = sample)) +
  geom_violin(fill = "#ff7f00", alpha = 0.6) +
  labs(
    title = "Genotype Depth Distribution by Sample",
    subtitle = "Truncated at 99th percentile",
    x = "Genotype Depth (FORMAT/DP)",
    y = NULL
  )
```

', dp_file[1]))
  } else {
    rmd_content <- paste0(rmd_content, '
# Genotype Depth (FORMAT/DP)

*Note: Depth file is too large (>100 MB) for detailed analysis in this report.*

')
  }
}

# No session info needed

# Write RMarkdown file
rmd_file <- file.path(input_dir, paste0(file_prefix, "_report.Rmd"))
writeLines(rmd_content, rmd_file)

# Render the report
output_file <- file.path(
  input_dir,
  paste0(file_prefix, "_report.", output_format)
)

cat(sprintf("Generating %s report...\n", toupper(output_format)))
cat(sprintf("Input directory: %s\n", input_dir))
cat(sprintf("Output file: %s\n", output_file))

tryCatch({
  rmarkdown::render(
    input = rmd_file,
    output_format = ifelse(output_format == "html", "html_document", "pdf_document"),
    output_file = basename(output_file),
    quiet = FALSE
  )
  
  cat(sprintf("\n✓ Report successfully generated: %s\n", output_file))
  
  # Optionally remove the intermediate Rmd file
  # file.remove(rmd_file)
  
}, error = function(e) {
  cat(sprintf("\n✗ Error generating report: %s\n", e$message))
  quit(status = 1)
})
