#!/usr/bin/env Rscript

#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Input directory"),
  make_option(c("-ctr", "--control_signal"), type = "character", help = "Control tag for signal"),
  make_option(c("-cond", "--condition_signal"), type = "character", help = "Condition tag for signal"),
  make_option(c("-c1", "--cond1_peaks"), type = "character", help = "Condition 1 peaks"),
  make_option(c("-c2", "--cond2_peaks"), type = "character", help = "Condition 2 peaks"),
  make_option(c("-cm", "--common_peaks"), type = "character", help = "Common peaks")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_dir <- opt$input_dir
control_signal <- opt$control_signal
condition_signal <- opt$condition_signal
cond1_peaks <- opt$cond1_peaks
cond2_peaks <- opt$cond2_peaks
common_peaks <- opt$common_peaks

max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)

#-----------------------------------------------------------------------------#
# LOAD DATA
#-----------------------------------------------------------------------------#
control_signal <- list.files(input_dir, pattern = control_signal, full.names = TRUE)
control_signal <- grep(".bg", control_signal, value = TRUE)
control_signal <- rtracklayer::import(control_signal)

condition_signal <- list.files(input_dir, pattern = condition_signal, full.names = TRUE)
condition_signal <- grep(".bg", condition_signal, value = TRUE)
condition_signal <- rtracklayer::import(condition_signal)

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
# 