#' Merge overlapping or nearby genomic peaks
#'
#' @description
#' Collapses overlapping or nearby intervals using `GenomicRanges::reduce()`,
#' then returns a mapping of each original `peak_id` to its merged interval ID
#' in the form `"chrom:start-end"`.
#'
#' @param df A data frame with columns:
#'   - `peak_id` (character): unique identifier per peak.
#'   - `chrom` (character): chromosome/contig name.
#'   - `start` (integer/numeric): 1-based start position.
#'   - `end` (integer/numeric): 1-based end position (≥ start).
#' @param max_gap Integer (≥ 0). Maximum distance (in bp) between intervals to
#'   consider them part of the same merged peak. `0` merges only overlapping
#'   ranges; `100` merges ranges ≤100 bp apart, etc.
#'
#' @return
#' A base data frame with columns:
#'   - `peak_id` (character): original peak ID.
#'   - `merged_peak` (character): merged interval ID `"chrom:start-end"`.
#'
#' @details
#' Uses `GenomicRanges` and `IRanges` without attaching the packages. The
#' reduction uses `min.gapwidth = max_gap + 1`, so `max_gap = 0` collapses only
#' overlaps.
#'
#' @examples
#' df <- data.frame(
#'   peak_id = paste0("p", 1:4),
#'   chrom   = c("chr1","chr1","chr1","chr2"),
#'   start   = c(100, 180, 1000, 50),
#'   end     = c(200, 220, 1100, 80),
#'   stringsAsFactors = FALSE
#' )
#' merge_peaks(df, max_gap = 50)
#'
#' @importFrom GenomicRanges GRanges reduce seqnames findOverlaps
#' @importFrom IRanges IRanges start end
#' @importFrom S4Vectors queryHits subjectHits
#' @export
merge_peaks <- function(df, max_gap = 0) {
  # Input check
  required_cols <- c("peak_id", "chrom", "start", "end")
  stopifnot(all(required_cols %in% colnames(df)))
  
  # Use GenomicRanges without loading the package
  gr <- GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    peak_id = df$peak_id
  )
  
  # Merge overlapping or nearby intervals
  reduced <- GenomicRanges::reduce(gr, min.gapwidth = max_gap + 1)
  
  # Create new merged peak identifiers
  reduced$merged_peak <- paste0(
    GenomicRanges::seqnames(reduced), ":",
    IRanges::start(reduced), "-",
    IRanges::end(reduced)
  )
  
  # Map original peaks to merged ones
  overlaps <- GenomicRanges::findOverlaps(gr, reduced)
  
  out_df <- data.frame(
    peak_id = gr$peak_id[S4Vectors::queryHits(overlaps)],
    merged_peak = reduced$merged_peak[S4Vectors::subjectHits(overlaps)],
    stringsAsFactors = FALSE
  )
  
  return(out_df)
}


#' Merge peaks per gene and map originals to gene-level spans
#'
#' @description
#' For each `gene_id` and chromosome, merges all peaks by taking the min start
#' and max end, then joins the merged spans back to the original peaks. Produces
#' a simple ID `merged_peak = "chrom_start_end"`.
#'
#' @param peaks_df A data frame/tibble with columns:
#'   - `peak` (character): original peak ID.
#'   - `chrom` (character): chromosome/contig name.
#'   - `start` (integer/numeric): 1-based start position.
#'   - `end` (integer/numeric): 1-based end position (≥ start).
#'   - `gene_id` (character): gene identifier.
#'
#' @return
#' A tibble/data frame with one row per original peak containing:
#' `peak, merged_peak, gene_id, chrom, start, end, merged_start, merged_end`.
#'
#' @details
#' This is a coarse gene-level merge (min start / max end). If you need
#' gap-aware merging, use `merge_peaks()` first within each gene.
#'
#' @examples
#' library(dplyr)
#' peaks <- tibble::tibble(
#'   peak    = paste0("p", 1:5),
#'   gene_id = c("g1","g1","g2","g2","g2"),
#'   chrom   = c("chr1","chr1","chr2","chr2","chr2"),
#'   start   = c(100, 220, 50, 120, 180),
#'   end     = c(150, 260, 90, 160, 210)
#' )
#' merge_peaks_by_gene(peaks)
#'
#' @importFrom dplyr group_by summarise mutate left_join select
#' @export
merge_peaks_by_gene <- function(peaks_df) {
  required_cols <- c("peak", "chrom", "start", "end", "gene_id")
  if (!all(required_cols %in% names(peaks_df))) {
    stop("Input must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Merge all peaks per gene_id + chrom
  merged_peaks <- peaks_df %>%
    dplyr::group_by(gene_id, chrom) %>%
    dplyr::summarise(
      merged_start = min(start),
      merged_end = max(end),
      .groups = "drop"
    ) %>%
    dplyr::mutate(merged_peak = paste0(chrom, "_", merged_start, "_", merged_end))
  
  # Join back to original to map old → new peak IDs
  peaks_df %>%
    dplyr::left_join(merged_peaks, by = c("gene_id", "chrom")) %>%
    dplyr::select(peak, merged_peak, gene_id, chrom, start, end, merged_start, merged_end)
}

