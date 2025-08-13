# test_peak_utils.R

test_that("merge_peaks() merges overlaps and respects max_gap", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_on_cran()
  
  df <- data.frame(
    peak_id = c("p1","p2","p3","p4"),
    chrom   = c("chr1","chr1","chr1","chr2"),
    start   = c(100,   140,   171,   50),
    end     = c(150,   160,   180,   60),
    stringsAsFactors = FALSE
  )
  
  # max_gap = 0 -> only overlaps are merged: (p1+p2), p3, p4
  out0 <- merge_peaks(df, max_gap = 0)
  expect_true(all(c("peak_id", "merged_peak") %in% names(out0)))
  
  exp0 <- data.frame(
    peak_id     = c("p1","p2","p3","p4"),
    merged_peak = c("chr1:100-160", "chr1:100-160", "chr1:171-180", "chr2:50-60"),
    stringsAsFactors = FALSE
  )
  out0 <- out0[order(out0$peak_id), , drop = FALSE]
  exp0 <- exp0[order(exp0$peak_id), , drop = FALSE]
  expect_identical(out0, exp0)
  
  # Gap between (p1+p2) end=160 and p3 start=171 is 10 bp.
  # max_gap = 9 -> do NOT merge p3; max_gap = 10 -> DO merge all chr1 peaks.
  out9  <- merge_peaks(df, max_gap = 9)
  out10 <- merge_peaks(df, max_gap = 10)
  
  exp9 <- exp0 # same as max_gap = 0 case
  out9 <- out9[order(out9$peak_id), , drop = FALSE]
  exp9 <- exp9[order(exp9$peak_id), , drop = FALSE]
  expect_identical(out9, exp9)
  
  exp10 <- data.frame(
    peak_id     = c("p1","p2","p3","p4"),
    merged_peak = c("chr1:100-180", "chr1:100-180", "chr1:100-180", "chr2:50-60"),
    stringsAsFactors = FALSE
  )
  out10 <- out10[order(out10$peak_id), , drop = FALSE]
  exp10 <- exp10[order(exp10$peak_id), , drop = FALSE]
  expect_identical(out10, exp10)
})

test_that("merge_peaks() validates input columns", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_on_cran()
  
  # Missing required columns
  bad <- data.frame(chrom = "chr1", start = 1, end = 10)
  expect_error(merge_peaks(bad), "all\\(required_cols %in%")
  

})

test_that("merge_peaks() keeps chromosomes separate", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  skip_on_cran()
  
  df <- data.frame(
    peak_id = c("a1","a2","b1","b2"),
    chrom   = c("chr1","chr1","chr2","chr2"),
    start   = c(10,  15,  100,  105),
    end     = c(12,  20,  101,  110),
    stringsAsFactors = FALSE
  )
  
  out <- merge_peaks(df, max_gap = 100)
  # Expect two merged peaks, one per chromosome
  merged_ids <- unique(out$merged_peak)
  expect_equal(length(merged_ids), 2L)
  expect_true(any(grepl("^chr1:", merged_ids)))
  expect_true(any(grepl("^chr2:", merged_ids)))
})

test_that("merge_peaks_by_gene() merges per gene+chrom and maps back", {
  skip_if_not_installed("dplyr")
  skip_on_cran()
  
  peaks <- data.frame(
    peak    = paste0("p", 1:6),
    gene_id = c("g1","g1","g1","g2","g2","g2"),
    chrom   = c("chr1","chr1","chr2","chr2","chr2","chr3"),
    start   = c(100, 220,  50,  120, 180, 10),
    end     = c(150, 260,  90,  160, 210, 40),
    stringsAsFactors = FALSE
  )
  
  out <- merge_peaks_by_gene(peaks)
  
  # Required columns present
  expect_setequal(
    names(out),
    c("peak","merged_peak","gene_id","chrom","start","end","merged_start","merged_end")
  )
  
  # Coarse gene+chrom spans are correct
  # g1 chr1 span: [100, 260]; g1 chr2 span: [50, 90]
  # g2 chr2 span: [120, 210]; g2 chr3 span: [10, 40]
  ref <- within(peaks, {
    merged_start <- c(100,100,50,120,120,10)
    merged_end   <- c(260,260,90,210,210,40)
    merged_peak  <- paste0(chrom, "_", merged_start, "_", merged_end)
  })
  # Align rows by peak id for deterministic compare
  out2 <- out[order(out$peak), , drop = FALSE]
  ref2 <- ref[order(ref$peak), , drop = FALSE]
  # Compare the mapped columns only
  expect_identical(
    as.data.frame(out2[, c("peak","gene_id","chrom","start","end","merged_start","merged_end","merged_peak")]),
    as.data.frame(ref2[, c("peak","gene_id","chrom","start","end","merged_start","merged_end","merged_peak")])
  )
})

test_that("merge_peaks_by_gene() validates input and supports single-peak genes", {
  skip_if_not_installed("dplyr")
  skip_on_cran()
  
  bad <- data.frame(gene_id = "g1", chrom = "chr1", start = 1, end = 10)
  expect_error(merge_peaks_by_gene(bad), "Input must contain columns")
  
  one <- data.frame(
    peak = "p1", gene_id = "gA", chrom = "chr5", start = 1000, end = 2000,
    stringsAsFactors = FALSE
  )
  out <- merge_peaks_by_gene(one)
  expect_equal(nrow(out), 1L)
  expect_equal(out$merged_start, 1000)
  expect_equal(out$merged_end, 2000)
  expect_equal(out$merged_peak, "chr5_1000_2000")
})
