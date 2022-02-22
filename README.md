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
In the primary analysis modality, MutDens requires two input files, one of somatic mutations and the other of focal genomic positions.  
### mutation file
A mutation file can be a VCF file or a TSV file. If in a TSV file, it should includes four columns for chromosome, position, ref allele, and alternate allele, sequentially. The TSV file should not contain a header row. Example TSV files and VCF files are provided in folder *mutFiles*.  
### focal position file
A focal position file designates point positions in the genome where we will investigate the mutational density trends in the flanking regions. A focal position file should be tab-delimited, and it should contain at least three columns (chromosome, position, strand) with a header row. If the user only cares about positions on the forward strand, just fill the strand column with all "+"s.  
In the two comparison analysis modalities, MutDens may require two mutation files (between-cohort) or two focal position files (between-region).  The necessary input files are specificed in an option file, along with other adjustable options.  

## Options
All options must be specified in an option file. The primary analysis entails 10 options, whereas the other two comparison analyses entail 11 options.

Key | Meaning | Example value
----|---------|--------------
**baseCont** | Base content normalization choice. | *{Origin, TSS, genome}*
**bsz** | Bin size. | *100*
**fTSS** | Feature is TSS or not? | *{TRUE, FALSE}*
**inspan** | Proximal boundary of farther flanking regions | *1000*
**mutF** | Mutation file name | *mutFiles/GACA-JP.tsv*
**outerspan** | Distal boundary of farther flanking regions | *6000*
**pointsF** | Focal position file name | *data/HG19ISq1.tsv*
**sbs** | Single base substitution class | *{C2A, C2G, C2T, T2A, T2C, T2G, SBS6}*
**span** | Length of flanking region in either direction | *2000*

## Output Description

