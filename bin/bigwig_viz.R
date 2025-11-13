#!/usr/bin/env Rscript

#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(argparser)
library(future)
library(future.apply)
#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
# Create argument parser
p <- arg_parser("Process peak signal data")

# Add arguments
p <- add_argument(p, "--input_dir", short = "-i", help = "Input directory", type = "character")
p <- add_argument(p, "--control_signal", short = "-ctr", help = "Control tag for signal", type = "character")
p <- add_argument(p, "--condition_signal", short = "-cond", help = "Condition tag for signal", type = "character")
p <- add_argument(p, "--cond1_peaks", short = "-c1", help = "Condition 1 peaks", type = "character")
p <- add_argument(p, "--cond2_peaks", short = "-c2", help = "Condition 2 peaks", type = "character")
p <- add_argument(p, "--common_peaks", short = "-cm", help = "Common peaks", type = "character")
p <- add_argument(p, "--max_width", short = "-max_w", help = "Max peak width", type = "integer", default = 1000)
p <- add_argument(p, "--display_width", short = "-dw", help = "Signal display window", type = "integer", default = 5000)
p <- add_argument(p, "--bin_width", short = "-bw", help = "Signal bin width", type = "integer", default = 100)
p <- add_argument(p, "--cores", short = "-c", help = "Number of cores", type = "integer", default = 1)

# Parse arguments
argv <- parse_args(p)

# Assign to variables
input_dir <- argv$input_dir
control_signal <- argv$control_signal
condition_signal <- argv$condition_signal
cond1_peaks <- argv$cond1_peaks
cond2_peaks <- argv$cond2_peaks
common_peaks <- argv$common_peaks
max_w <- argv$max_width
dw <- argv$display_width
bw <- argv$bin_width
cores <- argv$cores

max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)
plan(multicore, workers = cores)
#-----------------------------------------------------------------------------#
# UTILS
#-----------------------------------------------------------------------------#
filter_peaks <- function(peaks, max_w){
  locs <- width(peaks) <= max_w
  peaks <- peaks[locs]
  return(peaks)
}

get_peak_center <- function(peaks, dw) {
  semi_window <- dw %/% 2
  peak_start <- start(peaks)
  peak_width <- width(peaks) %/% 2
  peak_center <- peak_start + peak_width
  gr <- GRanges(seqnames = seqnames(peaks),
    ranges = IRanges(peak_center - semi_window, peak_center + semi_window),
    strand = "*")
  return(gr)
}
bin_signal <- function(signal, peak_centered, bw) {
  bins <- future_lapply(seq_along(peak_centered), function(i) {
    # Create all bins at once
    vec <- seq(start(peak_centered[i]), end(peak_centered[i]), by = bw)
    n_vec <- length(vec) - 1
    
    if (n_vec == 0) return(numeric(0))
    
    # Create GRanges for all bins at once
    bins_gr <- GRanges(
      seqnames = rep(seqnames(peak_centered[i]), n_vec),
      ranges = IRanges(start = vec[-length(vec)], end = vec[-1]),
      strand = "*"
    )
    
    # Single findOverlaps call
    overlaps <- findOverlaps(signal, bins_gr)
    
    # Initialize result
    loc_signal <- numeric(n_vec)
    
    if (length(overlaps) > 0) {
      # Group by bin and calculate means
      query_idx <- queryHits(overlaps)
      subject_idx <- subjectHits(overlaps)
      scores <- signal$score[query_idx]
      
      # Use tapply for grouping (often faster than aggregate)
      means <- tapply(scores, subject_idx, mean)
      loc_signal[as.integer(names(means))] <- means
    }
    
    return(loc_signal)
  })
  
  return(do.call("rbind", bins))
}
#-----------------------------------------------------------------------------#
# LOAD DATA
#-----------------------------------------------------------------------------#
control_signal <- list.files(input_dir, pattern = control_signal, full.names = TRUE)
control_signal <- grep(".bg", control_signal, value = TRUE)
control_signal <- rtracklayer::import(control_signal, format = "bedGraph")

condition_signal <- list.files(input_dir, pattern = condition_signal, full.names = TRUE)
condition_signal <- grep(".bg", condition_signal, value = TRUE)
condition_signal <- rtracklayer::import(condition_signal, format = "bedGraph")

cond1_peaks <- list.files(input_dir, pattern = cond1_peaks, full.names = TRUE)
cond1_peaks <- grep(".bed", cond1_peaks, value = TRUE)
cond1_peaks <- rtracklayer::import(cond1_peaks)

cond2_peaks <- list.files(input_dir, pattern = cond2_peaks, full.names = TRUE)
cond2_peaks <- grep(".bed", cond2_peaks, value = TRUE)
cond2_peaks <- rtracklayer::import(cond2_peaks)

common_peaks <- list.files(input_dir, pattern = common_peaks, full.names = TRUE)
common_peaks <- grep(".bed", common_peaks, value = TRUE)
common_peaks <- rtracklayer::import(common_peaks)

#-----------------------------------------------------------------------------#
# WRANGLE
#-----------------------------------------------------------------------------#
cond1_peaks_filtered <- filter_peaks(cond1_peaks, max_w)
cond1_peaks_centered <- get_peak_center(cond1_peaks_filtered, dw)
cond1_peaks_matrix <- bin_signal(control_signal,cond1_peaks_centered, bw)