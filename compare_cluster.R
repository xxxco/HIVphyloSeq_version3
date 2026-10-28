#!/usr/bin/env Rscript
# ============================================================
# compare_cluster_full.R  (REVISED 2 - multi-metric)
#
# Design UNCHANGED: subregion clustering vs full-length clustering,
# within each method (ClusterPicker vs full-length CP reference;
# HIV-TRACE vs full-length HIV-TRACE reference).
#
# METRIC CHANGE (this revision):
#   V-measure is dominated by the large mutually-singleton background
#   that is normal in HIV clustering, so it saturates near 1 and hides
#   real differences between regions/thresholds. We therefore also
#   report:
#     * ARI  - Adjusted Rand Index (chance-corrected; not singleton-inflated)
#     * pairwise precision / recall / F1 on co-clustered ("linked") pairs
#         - pw_recall    = sensitivity: fraction of full-length linked pairs
#                          the subregion also links (structure recovered)
#         - pw_precision = PPV: fraction of subregion linked pairs that are
#                          also linked in full-length (spurious links)
#   V-measure is reported three ways:
#     * V_measure_bothclustered = HEADLINE. Scored ONLY over sequences that are
#         in a cluster (>=2 members in the shared universe) in BOTH the subregion
#         and the full-length reference, so neither side's unclustered sequences
#         form a large group that biases the value. n_both_clustered = its size.
#     * V_measure_refclustered = scored over sequences in a full-length cluster
#         (>=2 members), whether or not the subregion clustered them, so it also
#         scores full-length clusters the subregion missed entirely.
#     * V_measure_full = legacy full-universe value (singleton-inflated); kept
#         only for contrast, NOT for interpretation.
#
# Earlier fixes retained:
#   1. parse_trace() reads Nodes$cluster as a flat array (compact_json)
#      AND the older Nodes$cluster$values layout; fails loudly (NULL)
#      rather than collapsing to one cluster.
#   2. Metrics computed over the FULL universe (clustered + singletons),
#      singletons given unique labels on both sides.
#   3. The sequence "universe" is read from the input alignment FASTAs.
# ============================================================

options(bspm.sudo = TRUE)

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(clevr)
  library(readr)
})

# ─── CONFIGURATION ───────────────────────────────────────────────────

desired_GD <- seq(from = 0.005, to = 0.045, by = 0.005)
BS_TAGS <- c("bs70", "bs90", "bs99")
BS_NUMBER_FOLDER <- "bs_1000"

PROJECT_ROOT <- "/hpc/home/ql165/projects/HIV_V3"
Year         <- "2005-2024_US_CA"

CP_BASE    <- file.path(PROJECT_ROOT, "Cluster", "ClusterPicker", Year)
TRACE_BASE <- file.path(PROJECT_ROOT, "Cluster", "HIV_trace", Year)
# source of the sequence universe: the clustering inputs (duplicate patients
# removed; Cluster/prepare_cluster_inputs.R), so it matches what was clustered
ALN_BASE   <- file.path(PROJECT_ROOT, "Cluster", "input", Year)

CP_LOG_FILE <- "input_clusterPicks_log.txt"

OUTPUT_DIR <- file.path(PROJECT_ROOT, "Cluster", "V_measure_Comparisons",
                        "V_measure_Dual_Ref", Year)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Full-length reference ---
# Default matches the manuscript's reference configuration (GD=3.0%, BS=90%
# for CP; GD=3.0% for HIV-TRACE). The reference GD can be overridden via a
# 4th CLI arg (see arg-parsing block below) to support the reviewer-requested
# sensitivity analysis: re-scoring every subregion against full-length
# partitions defined at OTHER GD thresholds, not just 3.0%. Reference BS
# stays fixed at bs90 in all cases -- the reviewer's request was about the
# GD threshold, not bootstrap support.
FULL_LENGTH_REF_REGION_NAME <- "fulllength"
FULL_LENGTH_REF_GD_NUMERIC <- 0.030
FULL_LENGTH_REF_BS_TAG <- "bs90"
# FULL_LENGTH_REF_GD_TAG, FULL_LENGTH_CP_REF_LOG_PATH, and
# FULL_LENGTH_TRACE_REF_JSON_PATH are built further down, AFTER argument
# parsing, so a 4th CLI arg can override FULL_LENGTH_REF_GD_NUMERIC first.

# ─── HELPERS ─────────────────────────────────────────────────────────

# Read the sequence-ID universe from a FASTA, normalised the SAME way the
# parsers normalise IDs (drop ">", keep first token, strip a leading "Seq_").
read_fasta_ids <- function(path) {
  if (!file.exists(path)) {
    message("  FASTA missing (cannot build sequence universe): ", path)
    return(character(0))
  }
  hdr <- grep("^>", readLines(path, warn = FALSE), value = TRUE)
  ids <- sub("^>", "", hdr)
  ids <- sub("\\s.*$", "", ids)     # first whitespace-delimited token
  sub("^Seq_", "", ids)
}

fasta_for <- function(region) {
  file.path(ALN_BASE, paste0(region, ".fasta"))
}

# Expand a parsed cluster table (clustered sequences only) to the full
# universe: every universe ID gets its cluster label; every ID not in a
# cluster gets its OWN unique singleton label. Returns an integer vector
# named by SequenceID, in universe_ids order.
build_labels <- function(df, universe_ids, cluster_col) {
  lab <- setNames(rep(NA_integer_, length(universe_ids)), universe_ids)
  if (!is.null(df) && nrow(df) > 0) {
    keep <- df$SequenceID %in% universe_ids
    lab[df$SequenceID[keep]] <- as.integer(df[[cluster_col]][keep])
  }
  max_c <- suppressWarnings(max(lab, na.rm = TRUE))
  if (!is.finite(max_c)) max_c <- 0L
  na_idx <- which(is.na(lab))
  if (length(na_idx) > 0) lab[na_idx] <- max_c + seq_along(na_idx)
  lab
}

# Pairwise precision / recall / F1 on co-clustered ("linked") pairs.
# Uses the contingency counts (O(N) memory: only non-empty joint cells are
# enumerated) rather than iterating over all N-choose-2 pairs.
pairwise_prf <- function(pred, ref) {
  comb2 <- function(x) x * (x - 1) / 2
  jt <- dplyr::count(tibble(pred = pred, ref = ref), pred, ref, name = "n")
  tp <- sum(comb2(jt$n))                              # linked in BOTH
  ap <- as.integer(table(pred)); pairs_pred <- sum(comb2(ap))   # linked in pred
  ar <- as.integer(table(ref));  pairs_ref  <- sum(comb2(ar))   # linked in ref
  precision <- if (pairs_pred > 0) tp / pairs_pred else NA_real_
  recall    <- if (pairs_ref  > 0) tp / pairs_ref  else NA_real_
  f1 <- if ((pairs_pred + pairs_ref) > 0) 2 * tp / (pairs_pred + pairs_ref) else NA_real_
  list(pw_precision = precision, pw_recall = recall, pw_f1 = f1,
       linked_pairs_shared = tp, linked_pairs_pred = pairs_pred, linked_pairs_ref = pairs_ref)
}

parse_cp <- function(log_path) {
  if (!file.exists(log_path)) {
    message(paste("ClusterPicker log missing:", log_path))
    return(NULL)
  }
  lines <- readLines(log_path, warn = FALSE)
  rows <- grep("^\\d+\\t", lines, value = TRUE)
  if (!length(rows)) {
    # No cluster rows. Two different cases, which must not be conflated:
    #   (a) ClusterPicker ran and found 0 clusters -- a valid result (common at
    #       low GD with a high bootstrap threshold). Every sequence is then a
    #       singleton, and the metrics are well defined, so return an EMPTY
    #       partition. Returning NULL here would silently drop the whole
    #       region/BS/GD row from the results table.
    #   (b) the log is truncated or unparseable -> NULL, so the row is skipped.
    if (any(grepl("^Found\\t\\d+\\tclusters", lines))) {
      message(paste("ClusterPicker found 0 clusters (all singletons):", log_path))
      return(tibble(SequenceID = character(), CP_cluster = integer()))
    }
    message(paste("No results table in ClusterPicker log (truncated?):", log_path))
    return(NULL)
  }
  bind_rows(lapply(rows, function(r) {
    p <- str_split(r, "\\t", n = 6)[[1]]
    tibble(
      SequenceID = str_remove_all(p[4], "\\[|\\]|\\s") |>
        str_split(",") |> unlist() |> str_replace("^Seq_", ""),
      CP_cluster = as.integer(p[1])
    )
  }))
}

# Handles compact_json (flat array) and older nested ($values) layouts,
# and fails loudly rather than collapsing to one cluster.
parse_trace <- function(json_path) {
  if (!file.exists(json_path)) {
    message(paste("HIV-TRACE JSON missing:", json_path))
    return(NULL)
  }
  j <- tryCatch(fromJSON(json_path, simplifyVector = FALSE),
                error = function(e) {
                  message(paste("Error parsing JSON for", json_path, ":", e$message)); NULL
                })
  if (is.null(j) || is.null(j$trace_results) || is.null(j$trace_results$Nodes)) {
    message(paste("Invalid or empty HIV-TRACE JSON structure for:", json_path))
    return(NULL)
  }
  
  nodes <- j$trace_results$Nodes
  
  cluster_values <- nodes$cluster
  if (is.list(cluster_values) && !is.null(cluster_values$values)) {
    cluster_values <- cluster_values$values
  }
  cluster_ids <- nodes$id
  if (is.list(cluster_ids) && !is.null(cluster_ids$values)) {
    cluster_ids <- cluster_ids$values
  }
  
  if (is.null(cluster_values) || length(cluster_values) == 0 ||
      is.null(cluster_ids)    || length(cluster_ids) == 0) {
    message(paste("Could not extract Nodes$cluster / Nodes$id from:", json_path,
                  "-- check the JSON layout (compact_json?). Skipping (NOT collapsing)."))
    return(NULL)
  }
  if (length(cluster_values) != length(cluster_ids)) {
    message(paste("cluster/id length mismatch in:", json_path, "-- skipping."))
    return(NULL)
  }
  
  tibble(
    SequenceID    = unlist(cluster_ids) |> as.character() |> str_replace("^Seq_", ""),
    TRACE_cluster = unlist(cluster_values) |> as.integer()
  )
}

# Reference-clustered evaluation set: sequences whose full-length reference
# cluster has >= 2 members PRESENT in this universe. Restricting to this set
# removes the trivially-perfect singleton mass that pins the full-universe
# V-measure near 1. (Ref clusters that collapse to a single member within the
# universe are dropped, since a lone member would re-introduce a singleton.)
ref_clustered_ids <- function(ref_df, universe_ids, ref_col) {
  ri <- ref_df[ref_df$SequenceID %in% universe_ids, , drop = FALSE]
  if (nrow(ri) == 0) return(character(0))
  tab <- table(ri[[ref_col]])
  big <- names(tab)[tab >= 2]
  unique(ri$SequenceID[as.character(ri[[ref_col]]) %in% big])
}

# Compute all metrics.
#   V_measure_bothclustered = HEADLINE (clustered by both)
#   V_measure_refclustered = ref-clustered set (also scores misses)
#   V_measure_full         = legacy full-universe value (singleton-inflated)
#   ARI / pairwise         = full universe (already robust to singletons)
calc_metrics <- function(pred_df, ref_df, universe_ids, pred_col, ref_col) {
  na_out <- list(V_measure_bothclustered = NA_real_, n_both_clustered = 0L,
                 V_measure_refclustered = NA_real_, V_measure_full = NA_real_,
                 ARI = NA_real_, pw_precision = NA_real_, pw_recall = NA_real_,
                 pw_f1 = NA_real_, linked_pairs_shared = NA_real_,
                 linked_pairs_pred = NA_real_, linked_pairs_ref = NA_real_,
                 n_ref_clustered = 0L)
  if (length(universe_ids) == 0) {
    message("  Empty sequence universe -> metrics NA.")
    return(na_out)
  }
  
  # Full-universe labels for ARI / pairwise / legacy V-measure.
  pred_labels <- build_labels(pred_df, universe_ids, pred_col)
  ref_labels  <- build_labels(ref_df,  universe_ids, ref_col)
  prf <- pairwise_prf(pred_labels, ref_labels)
  
  # Headline V-measure: sequences clustered (in a cluster of >= 2 within the
  # universe) by BOTH the subregion and the full-length reference, so neither
  # side's unclustered sequences form a large group that biases the value.
  # ref_clustered_ids() works for either side: it keeps members of clusters
  # with >= 2 members in the given set.
  eval_both <- intersect(ref_clustered_ids(ref_df,  universe_ids, ref_col),
                         ref_clustered_ids(pred_df, universe_ids, pred_col))
  vm_both <- NA_real_
  if (length(eval_both) >= 2) {
    vm_both <- clevr::v_measure(build_labels(pred_df, eval_both, pred_col),
                                build_labels(ref_df,  eval_both, ref_col))
  }

  # Secondary V-measure over the reference-clustered set: also scores
  # full-length clusters the subregion missed (pred singletons = misses).
  eval_ids <- ref_clustered_ids(ref_df, universe_ids, ref_col)
  vm_ref <- NA_real_
  if (length(eval_ids) >= 2) {
    pred_e <- build_labels(pred_df, eval_ids, pred_col)  # pred singletons here = real misses
    ref_e  <- build_labels(ref_df,  eval_ids, ref_col)   # no ref singletons by construction
    vm_ref <- clevr::v_measure(pred_e, ref_e)
  }

  list(
    V_measure_bothclustered = vm_both,
    n_both_clustered        = length(eval_both),
    V_measure_refclustered = vm_ref,
    V_measure_full         = clevr::v_measure(pred_labels, ref_labels),
    ARI                    = clevr::adj_rand_index(pred_labels, ref_labels),
    pw_precision           = prf$pw_precision,
    pw_recall              = prf$pw_recall,
    pw_f1                  = prf$pw_f1,
    linked_pairs_shared    = prf$linked_pairs_shared,
    linked_pairs_pred      = prf$linked_pairs_pred,
    linked_pairs_ref       = prf$linked_pairs_ref,
    n_ref_clustered        = length(eval_ids)
  )
}

# ─── MAIN SCRIPT EXECUTION ───────────────────────────────────────────

# Usage:
#   Rscript compare_cluster_full.R regions.txt                        # full sweep
#   Rscript compare_cluster_full.R regions.txt bs90                   # one bootstrap, all GDs
#   Rscript compare_cluster_full.R regions.txt bs90 0.03               # one bootstrap, one GD
#   Rscript compare_cluster_full.R regions.txt "" 0.03                 # all bootstraps, one GD
#   Rscript compare_cluster_full.R regions.txt "" "" 0.02              # all bootstraps, all GDs,
#                                                                       #   reference full-length GD = 2.0%
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript compare_cluster_full.R regions.txt [bs_tag] [gd_fraction] [ref_gd_fraction]")
}

region_file <- args[1]
regions <- readLines(region_file, warn = FALSE)

restrict_bs     <- length(args) >= 2 && nzchar(args[2])
restrict_gd     <- length(args) >= 3 && nzchar(args[3])
restrict_ref_gd <- length(args) >= 4 && nzchar(args[4])

if (restrict_bs) {
  BS_TAGS <- args[2]
  if (!grepl("^bs", BS_TAGS)) {
    warning(sprintf("bs_tag '%s' does not start with 'bs' -- folders use e.g. 'bs90'.", BS_TAGS))
  }
}
if (restrict_gd) {
  desired_GD <- as.numeric(args[3])
  if (is.na(desired_GD)) stop(sprintf("Could not parse GD '%s' as a number.", args[3]))
  if (desired_GD > 0.5)  stop(sprintf("GD '%s' looks like a percent; pass a fraction (0.03).", args[3]))
}
if (restrict_ref_gd) {
  ref_gd_arg <- as.numeric(args[4])
  if (is.na(ref_gd_arg)) stop(sprintf("Could not parse reference GD '%s' as a number.", args[4]))
  if (ref_gd_arg > 0.5)  stop(sprintf("Reference GD '%s' looks like a percent; pass a fraction (0.02).", args[4]))
  FULL_LENGTH_REF_GD_NUMERIC <- ref_gd_arg
}

# Reference path/tag construction -- AFTER arg parsing, so a 4th CLI arg
# (if given) has already overridden FULL_LENGTH_REF_GD_NUMERIC.
FULL_LENGTH_REF_GD_TAG <- sprintf("gd%03d", as.integer(round(FULL_LENGTH_REF_GD_NUMERIC * 1000)))

FULL_LENGTH_CP_REF_LOG_PATH <- file.path(
  CP_BASE,
  paste0(FULL_LENGTH_REF_REGION_NAME, "_CP"),
  BS_NUMBER_FOLDER,
  paste0(FULL_LENGTH_REF_BS_TAG, "_", FULL_LENGTH_REF_GD_TAG),
  CP_LOG_FILE
)

FULL_LENGTH_TRACE_REF_JSON_PATH <- file.path(
  TRACE_BASE,
  FULL_LENGTH_REF_REGION_NAME,
  paste0(FULL_LENGTH_REF_REGION_NAME, "_hivtrace_GD",
         sprintf("%03d", as.integer(round(FULL_LENGTH_REF_GD_NUMERIC * 1000))), ".json")
)

message("\n--- Starting cluster-concordance comparisons (multi-metric) ---\n")
message(sprintf("Project root : %s", PROJECT_ROOT))
message(sprintf("Year         : %s", Year))
message(sprintf("Bootstraps   : %s", paste(BS_TAGS, collapse = ", ")))
message(sprintf("GDs          : %s", paste(sprintf("%.3f", desired_GD), collapse = ", ")))
message(sprintf("Reference GD : %.3f (BS=%s)", FULL_LENGTH_REF_GD_NUMERIC, FULL_LENGTH_REF_BS_TAG))
message("Metrics      : V-measure over clustered-by-both set (headline), ref-clustered set + ARI + pairwise; full V-measure (legacy)")

# --- Load full-length references ---
ref_cp_df <- parse_cp(FULL_LENGTH_CP_REF_LOG_PATH)
if (is.null(ref_cp_df)) stop(sprintf("ERROR: full-length CP reference not loaded from %s", FULL_LENGTH_CP_REF_LOG_PATH))
message(sprintf("Loaded %d full-length CP reference assignments (GD=%.1f%%, BS=%s).",
                nrow(ref_cp_df), FULL_LENGTH_REF_GD_NUMERIC * 100, FULL_LENGTH_REF_BS_TAG))

ref_trace_df <- parse_trace(FULL_LENGTH_TRACE_REF_JSON_PATH)
if (is.null(ref_trace_df)) stop(sprintf("ERROR: full-length HIV-TRACE reference not loaded from %s", FULL_LENGTH_TRACE_REF_JSON_PATH))
message(sprintf("Loaded %d full-length HIV-TRACE reference assignments (%d clusters).",
                nrow(ref_trace_df), length(unique(ref_trace_df$TRACE_cluster))))

# Reference universe = IDs present in the full-length alignment
ref_universe <- read_fasta_ids(fasta_for(FULL_LENGTH_REF_REGION_NAME))

all_results <- list()

for (region in regions) {
  region <- trimws(region)
  if (region == "" || startsWith(region, "#")) next
  message(sprintf("\n=== Processing Region: %s ===", region))
  
  region_universe <- read_fasta_ids(fasta_for(region))
  universe_ids <- intersect(region_universe, ref_universe)
  if (length(universe_ids) == 0) {
    message(sprintf("  No shared sequence universe for %s; skipping region.", region))
    next
  }
  message(sprintf("  Sequence universe (region ∩ full-length): %d", length(universe_ids)))
  
  for (gd_numeric_val in desired_GD) {
    gd_tag_str <- sprintf("gd%03d", as.integer(round(gd_numeric_val * 1000)))
    
    # ── A. ClusterPicker subregion vs full-length CP reference ──
    for (bs_tag in BS_TAGS) {
      cp_log_path <- file.path(CP_BASE, paste0(region, "_CP"), BS_NUMBER_FOLDER,
                               paste0(bs_tag, "_", gd_tag_str), CP_LOG_FILE)
      current_cp_df <- parse_cp(cp_log_path)
      if (is.null(current_cp_df)) {
        message(sprintf("  Skipping CP (BS=%s, GD=%.3f) for %s", bs_tag, gd_numeric_val, region))
      } else {
        m <- calc_metrics(current_cp_df, ref_cp_df, universe_ids, "CP_cluster", "CP_cluster")
        all_results[[length(all_results) + 1]] <- tibble(
          Region = region, Method = "ClusterPicker",
          GD_percent = gd_numeric_val * 100, Bootstrap = bs_tag,
          Reference_Dataset = "Full-Length", Reference_Method = "ClusterPicker",
          Reference_GD = FULL_LENGTH_REF_GD_NUMERIC * 100,
          Reference_Bootstrap = FULL_LENGTH_REF_BS_TAG,
        V_measure_refclustered = round(m$V_measure_refclustered, 6),  # HEADLINE: also scores missed ref clusters
          pw_recall    = round(m$pw_recall, 6),      # sensitivity to full-length structure
          pw_precision = round(m$pw_precision, 6),   # PPV vs full-length
          pw_f1        = round(m$pw_f1, 6),
          V_measure_full = round(m$V_measure_full, 6),  # legacy (singleton-inflated)
          n_ref_clustered     = m$n_ref_clustered,      # seqs the headline V-measure rests on
          linked_pairs_shared = m$linked_pairs_shared,
          linked_pairs_pred   = m$linked_pairs_pred,
          linked_pairs_ref    = m$linked_pairs_ref,
          n_predicted_clusters = length(unique(current_cp_df$CP_cluster)),
          n_reference_clusters = length(unique(ref_cp_df$CP_cluster)),
          n_predicted_seqs = nrow(current_cp_df),
          n_reference_seqs = nrow(ref_cp_df),
          n_universe       = length(universe_ids),
          n_singletons     = length(universe_ids) - sum(current_cp_df$SequenceID %in% universe_ids)
        )
      }
    }
    
    # ── B. HIV-TRACE subregion vs full-length HIV-TRACE reference ──
    trace_json_path <- file.path(TRACE_BASE, region,
                                 paste0(region, "_hivtrace_GD",
                                        sprintf("%03d", as.integer(round(gd_numeric_val * 1000))), ".json"))
    current_trace_df <- parse_trace(trace_json_path)
    if (is.null(current_trace_df)) {
      message(sprintf("  Skipping HIV-TRACE (GD=%.3f) for %s", gd_numeric_val, region))
    } else {
      m <- calc_metrics(current_trace_df, ref_trace_df, universe_ids, "TRACE_cluster", "TRACE_cluster")
      all_results[[length(all_results) + 1]] <- tibble(
        Region = region, Method = "HIV-TRACE",
        GD_percent = gd_numeric_val * 100, Bootstrap = NA_character_,
        Reference_Dataset = "Full-Length", Reference_Method = "HIV-TRACE",
        Reference_GD = FULL_LENGTH_REF_GD_NUMERIC * 100, Reference_Bootstrap = NA_character_,
        V_measure_refclustered = round(m$V_measure_refclustered, 6),  # HEADLINE: also scores missed ref clusters
        pw_recall    = round(m$pw_recall, 6),
        pw_precision = round(m$pw_precision, 6),
        pw_f1        = round(m$pw_f1, 6),
        V_measure_full = round(m$V_measure_full, 6),  # legacy (singleton-inflated)
        n_ref_clustered     = m$n_ref_clustered,
        linked_pairs_shared = m$linked_pairs_shared,
        linked_pairs_pred   = m$linked_pairs_pred,
        linked_pairs_ref    = m$linked_pairs_ref,
        n_predicted_clusters = length(unique(current_trace_df$TRACE_cluster)),
        n_reference_clusters = length(unique(ref_trace_df$TRACE_cluster)),
        n_predicted_seqs = nrow(current_trace_df),
        n_reference_seqs = nrow(ref_trace_df),
        n_universe       = length(universe_ids),
        n_singletons     = length(universe_ids) - sum(current_trace_df$SequenceID %in% universe_ids)
      )
    }
  }
}

# ─── FINAL OUTPUT ────────────────────────────────────────────────────

if (length(all_results) > 0) {
  final_results_df <- bind_rows(all_results) %>%
    arrange(Region, GD_percent, Method, Bootstrap, Reference_Method)
  
  out_name <- "Dual_Ref_Results_full_new.tsv"
  if (restrict_bs || restrict_gd || restrict_ref_gd) {
    tag <- paste0(
      if (restrict_bs) paste0("_", BS_TAGS) else "",
      if (restrict_gd) sprintf("_gd%03d", as.integer(round(desired_GD * 1000))) else "",
      if (restrict_ref_gd) sprintf("_refgd%03d", as.integer(round(FULL_LENGTH_REF_GD_NUMERIC * 1000))) else ""
    )
    out_name <- paste0("Dual_Ref_Results", tag, ".tsv")
  }
  output_tsv_path <- file.path(OUTPUT_DIR, out_name)
  write_tsv(final_results_df, output_tsv_path)
  message(sprintf("\n✓ Wrote metric results to %s", output_tsv_path))
} else {
  message("\n⚠ No results generated.")
}

message("\n--- All comparisons complete ---")