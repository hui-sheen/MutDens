# MutDens: Mutational density spatial trend in focal regions

## Introduction
MutDens is an R application to thoroughly investigate mutational density patterns for any specific genomic regions. By scanning the bi-directional vicinity regions of the focal positions of any specific genomic feature, MutDens systematically characterizes the mutational density spatial trends for each of six mutational classes after adjusting for total mutation burden and local nucleotide proportion. MutDens is capable of revealing mutational density peaks, dips, and strand biases around genomic features such as transcription start sites (TSS) and replication origins (Origin).

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
The application bundle comprises focal position files for three genomic features:  TSS, Origin, and Retrotransposon Insertion Polymorphism (RIP) sites. 

file name | feature | reference genome
----------|---------|-----------------
gTSS37.tsv | TSS | GRCh37
gTSS38.tsv | TSS | GRCh38
HG19IS.tsv | Origin | GRCh37
HG19ISq1.tsv | Origin (1st quantile) | GRCh37
HG19RIP.tsv | Retrotransposon Insertion Polymorphism | GRCh37
HG38IS.tsv | Origin | GRCh38
HG38ISq1.tsv | Origin (1st quantile) | GRCh38

In the two comparison analysis modalities, MutDens may require two mutation files (between-cohort) or two focal position files (between-region).  The necessary input files are specificed in an option file, along with other adjustable options.  

## Options
All options must be specified in an option file. The primary analysis entails 10 options, whereas the other two comparison analyses entail 11 options. Template option files for the three analysis modalities are included in folder *optFiles*.

Key | Meaning | Example value
----|---------|--------------
**baseCont** | Base content normalization choice. | *{Origin, TSS, genome}*
**bsz** | Bin size. | *100*
**fTSS** | Feature is TSS or not? | *{TRUE, FALSE}*
**inspan** | Proximal boundary of farther flanking regions | *1000*
**mutF** | Mutation file name | *mutFiles/GACA-JP.tsv*
**mutF1** | Mutation file name 1 | *mutFiles/LICA-CN.tsv*
**mutF2** | Mutation file name 2 | *mutFiles/LICA-FR.tsv*
**outerspan** | Distal boundary of farther flanking regions | *6000*
**pointsF** | Focal position file name | *data/HG19ISq1.tsv*
**points1F** | Focal position file name 1 | *data/gTSS37.tsv*
**points2F** | Focal position file name 2 | *data/HG19RIP.tsv*
**sbs** | Single base substitution class | *{C2A, C2G, C2T, T2A, T2C, T2G, SBS6}*
**span** | Length of flanking region in either direction | *2000*

## Output Description

