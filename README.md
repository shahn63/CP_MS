# CP_MS
# Conformalized Outlier detection for Mass-sepctrometry data

# Description of content:
* CPMS.R:
  
  The main source code.
  CPMS function calculates the conformal prediction set and marginal conformal prediction intervals.
  'samp.plot=T' builds the plot of the conformality scores/MH distance across the samples and returns detected sample IDs.
  One needs to load this cource code via source("CPMS.R").

* CPMS_band.R:
  
  The main source code.
  CPMS_band function builds the plot of the marginal conformality intervals over m/z values with a given detected samples which is selected from the CPMS function.
  The inputs are the CPMS object and detected ID.
  It returns detected m/z value (peak) IDs.
  One needs to load this cource code via source("CPMS_band.R").

* MSexample.R:
  
  This is an example code to reproduce the analysis of the LC_MRM data presented in the manuscript.
  One needs to load the main source codes 'CPMS.R' and 'CPMS_band.R' via source("CPMS.R") and source("CPMS_band.R")

* README.txt:
  
  This file.

# R programming specifications:

* Session info:
It was tested with the following configureation:

R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
