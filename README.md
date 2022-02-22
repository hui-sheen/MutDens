# MutDens: Mutational density spatial trend in focal regions

## Introduction
MutDens is an R application to thoroughly investigate mutational density patterns for any specific genomic regions. By scanning the bi-directional vicinity regions of the focal positions of any specific genomic feature, MutDens systematically characterizes the mutational density spatial trends for each of six mutational classes of Single Base Substitutions (SBSs) after adjusting for total mutation burden and local nucleotide proportion. MutDens is capable of revealing mutational density peaks, dips, and strand biases around genomic features such as transcription start sites (TSS) and replication origins (Origin).

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
Generally, MutDens takes no more than a few minutes to complete the example analyses as defined by the three option files under *optFiles* (*options_1reg.txt*, *options_2cht.txt*, and *options_2reg.txt*). In the working directory, an HTML report that includes four sections is generated. In addition, a project folder named after distinctive option values is created, and several text files are saved in the project folder.  
### Mutational class distribution
One pie chart illustrates the distribution of SBSs across all six classes. One horizontal stacked barplot indicates how the two complementary SBS forms balance within a mutational class. Legend label "itself" designates the pyrimidine-originated SBS (those with C/T as the source), whereas "revmut" designates the complementary form that originates from purines (with G/A as the source).

### Mutational density curves
In the primary analysis modality, the figure comprises two mutational density curves across the bi-directional flanking regions. The X-axis shows the length of considered flanking nucleotides, and the points in the curves indicate how much mutational density values are in the running adjacent bins. Y-axis may be labelled as # Mutations per Kilo total mutations per Mega-base (MPKM), indicating the total mutational burden has been normalized against; otherwise, Y-axis is labelled as # Mutations per Mega-base (MPM). In the TSS context (fTSS==TRUE), the orange and blue colors designate the coding and template strands of the pyrimidine-initialized SBS, and they can reversely designate the template and coding strands of the complementary purine-initialized SBS. In the non-TSS context (fTSS==FALSE), the orange and blue colors designate the pyrimidine-initialized SBS and the complementary purine-initialized SBS, respectively. 

In the comparison modalities, the figure comprises four mutational density curves. The differentiation between Orange and Blue colors is the same as in the primary analysis modality. The solid/dotted lines are used to differentiate the two cohorts (between-cohort) or the two region sets (between-region). 

### Peaks/Dips
When the Poisson distribution test resulted in a p-value less than 1E-5, a mutational density peak or dip that may span multiple continuous bins is indicated with a light gray rectangle in the mutational density curve plot. A table is shown for each mutational class, which includes p-values and left/right boundary bin indices.

### Mutational density difference
In the primary analysis modality, Paired Wilcoxon test is invoked to compare mutational density values between two complementary forms across all bins. Three different ranges are considered: left wing, right wing, and whole range (both wings). In the comparison analysis modalities,  Paired Wilcoxon test is conducted between two cohorts (between-cohort) or two region sets (between-region), rather than between the two mutational forms.  





