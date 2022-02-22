# MutDens: Mutational density spatial trend in focal regions

## Introduction
MutDens is an R application to thoroughly investigate mutational density patterns for any specific genomic regions. By scanning the bi-directional vicinity regions of the focal positions of any specific genomic feature, MutDens systematically characterizes the mutational density spatial trends for each of six mutational classes after adjusting for total mutation burden and local nucleotide proportion. MutDens is capable of revealing mutational density peaks, dips, and strand biases around genomic features such as transcription start sites and replication origins.

## Usage
git clone https://github.com/hui-sheen/MutDens/
cd MutDens # the root directory of MutDens
R # open an R session
>library('rmarkdown')
>render('mutDens_1reg.Rmd') # Primary analysis of only 1 sample cohort
>render('mutDens_2cht.Rmd') # Comparison analysis between two sample cohorts
>render('mutDens_2reg.Rmd') # Comparison analysis between two region sets

## Input files

## Options

## Output Description

