# MutDens: Mutational density spatial trend in focal regions

## Introduction
MutDens is an R application to thoroughly investigate mutational density patterns for any specific genomic regions. By scanning the bi-directional vicinity regions of the focal positions of any specific genomic feature, MutDens systematically characterizes the mutational density spatial trends for each of six mutational classes of Single Base Substitutions (SBSs) after adjusting for total mutation burden and local nucleotide proportion. MutDens is capable of revealing mutational density peaks, dips, and strand biases around genomic features such as transcription start sites (TSS) and replication origins (Origin).  
![Concept Illustration](/fig/illustration.jpg)


## Usage
	git clone https://github.com/hui-sheen/MutDens/  
	cd MutDens # the root directory of MutDens  	
	R # open an R session  
    	>library('rmarkdown')  
    	>render('mutDens_1reg.Rmd') # Primary analysis of only 1 sample cohort, relying on *optFiles/options_1reg.txt*  
    	>render('mutDens_2cht.Rmd') # Comparison analysis between two sample cohorts, relying on *optFiles/options_2cht.txt*
    	>render('mutDens_2reg.Rmd') # Comparison analysis between two region sets, relying on *optFiles/options_2reg.txt*

## Input files
In the primary analysis modality, MutDens requires two input files, one of somatic mutations and the other of focal genomic positions.  
### mutation file
A mutation file can be a VCF file or a TSV file. If in a TSV file, it should includes four columns for chromosome, position, ref allele, and alternate allele, sequentially. The TSV file should not contain a header row. Example TSV files and VCF files are provided in folder *mutFiles*.  
### focal position file
A focal position file designates point positions in the genome where we will investigate the mutational density trends in the flanking regions. A focal position file should be tab-delimited, and it should contain at least two columns (chromosome, position, strand) with a header row. A column with header "strand" is optional, which values of "+" or "-" are given.    
For humans (hg38 and hg19), the *data/* folder contains coordinate files for five genomic features:  Transcription Start Site (TSS), Transcription End Site (TES), replication oritin Initiation Site(IS), enhancer (eRNA), and Retrotransposon Insertion Polymorphism (RIP) sites. TSS and TES are further distinguished into those of protein-coding genes (proTSS and proTES) and lncRNAs (lncTSS and lncTES), respectively. TSS and TES files are also available for the following other model organisms: canFam3, danRer11, dm6, galGal6, hg38, hg19, mm39, rheMac10, rn7, sacCer3. 
Please examine the files in *data/* for a full list of pre-built genomic features in several supported reference genomes. Users' self-prepared coordinate files for interesting features are recommended to be stored there too, following the same format as the prebuilt ones. Presently, MutDens tackles nine species as enumerated above. In principle, MutDens can be extended to tackle any species for which a ![BSgenome.*](https://bioconductor.org/packages/release/BiocViews.html#___AnnotationData) is released in Bioconductor.  
All prebuilt focal position files conform to this file naming format: **[gFt]_[gn].tsv**, where *gFt* designates the genomic feature and *gn* the code for the species or genome. A few example files are explained in the table below.  
file name | feature | reference genome
----------|---------|-----------------
eRNA_hg38.tsv | enhancer | GRCh38
ISq1_hg19.tsv | Replication Origin (1st quantile) | GRCh37
lncTES_hg38.tsv | TES of lncRNA genes | GRCh38
RIP_hg19.tsv | Retrotransposon Insertion Polymorphism | GRCh37
IS_hg38.tsv | Replication Origin | GRCh38
TSS_hg19.tsv | TSS | GRCh37
proTSS_hg38.tsv | TSS of protein-coding genes | GRCh38

In the two comparison analysis modalities, MutDens may require two mutation files (between-cohort) or two focal position files (between-region).  The necessary input files are specificed in an option file, along with other adjustable options.  

## Of Note: Base Content
By default (*calcSrc*=NULL, see below), MutDens seeks to calculate the vicinity base content (A/T/G/C proportion values) on the fly. The resultant base content file is named as **atgc2000-100_[gFt]_[gn].tsv**. The user can move the resultant atgc* file to data/, thereby next time MutDens can immediately locate this base content file for the *gFt* feature in the *gn* reference genome.  
## Options
All options must be specified in an option file. Template option files for the three analysis modalities are included in folder *optFiles*. Users can modify and rename these template options to reflect their real analysis requirements. An option file does not need to be put in *optFiles*, only that its absolute or relative path needs to be correctly specified in the analysis Rmd file.  
The primary analysis entails 12 options, whereas the other two comparison analyses entail 13 and 14 options. Despite the ostensible multitude, most options can use default values. The most important options involve file names for the mutations and the genomic feature(s), which must be prepared by the user.

Key | Meaning | Example value
----|---------|--------------
**baseCont** | Base content normalization choice. | *{IS, TSS, genome}*
**bsz** | Bin size. | *100*
**calcSrc** | Calculation source for base content in feature vicinity | *NULL*, or a file name
**gFt** | Genomic feature | *{IS, TSS, RIP, ...}*
**gFt1** | Genomic feature 1 | *{IS, TSS, RIP, ...}*
**gFt2** | Genomic feature 2 | *{IS, TSS, RIP, ...}*
**gn** | Genome | *{canFam3, danRer11, dm6, galGal6, hg38, hg19, mm39, rheMac10, rn7, sacCer3}*
**inspan** | Proximal boundary of farther flanking regions | *1000*
**mutF** | Mutation file name | *mutFiles/GACA-JP.tsv*
**mutF1** | Mutation file name 1 | *mutFiles/LICA-CN.tsv*
**mutF2** | Mutation file name 2 | *mutFiles/LICA-FR.tsv*
**outerspan** | Distal boundary of farther flanking regions | *6000*
**pointsF** | Focal position file name | *data/HG19ISq1.tsv*
**points1F** | Focal position file name 1 | *data/gTSS37.tsv*
**points2F** | Focal position file name 2 | *data/HG19RIP.tsv*
**sbs** | Single base substitution class | *{C2A, C2G, C2T, T2A, T2C, T2G, SBS6}*
**shapeModel** | Probability model to test for peak/dip | *{pois,nbinom}*
**span** | Length of flanking region in either direction | *2000*
**triCont** | Expand to trinucleotide context? | *{TRUE, FALSE}*

## Output Description
Generally, MutDens takes no more than a few minutes to complete the example analyses as defined by the three option files under *optFiles* (*options_1reg.txt*, *options_2cht.txt*, and *options_2reg.txt*). In the working directory, an HTML report that includes four sections is generated. In addition, a project folder named after distinctive option values is created, and several text files are saved in the project folder.  
### Mutational class distribution
One pie chart illustrates the distribution of SBSs across all six classes. One horizontal stacked barplot indicates how the two complementary SBS forms balance within a mutational class. Legend label "itself" designates the pyrimidine-originated SBS (those with C/T as the source), whereas "revmut" designates the complementary form that originates from purines (with G/A as the source).  
![Example class distribution figures](/fig/classDistrib.JPG)

### Mutational density curves
In the primary analysis modality, the figure comprises two mutational density curves across the bi-directional flanking regions. The X-axis shows the length of considered flanking nucleotides, and the points in the curves indicate how much mutational density values are in the running adjacent bins. Y-axis may be labelled as # Mutations per Kilo total mutations per Mega-base (MPKM), indicating the total mutational burden has been normalized against; otherwise, Y-axis is labelled as # Mutations per Mega-base (MPM). In the TSS context (fTSS==TRUE), the orange and blue colors designate the coding and template strands of the pyrimidine-initialized SBS, and they can reversely designate the template and coding strands of the complementary purine-initialized SBS. In the non-TSS context (fTSS==FALSE), the orange and blue colors designate the pyrimidine-initialized SBS and the complementary purine-initialized SBS, respectively.  

![Two mutational density series](/fig/Curves2.JPG)

In the comparison modalities, the figure comprises four mutational density curves. The differentiation between Orange and Blue colors is the same as in the primary analysis modality. The solid/dotted lines are used to differentiate the two cohorts (between-cohort) or the two region sets (between-region). 

![Four mutational density series](/fig/Curves4.JPG)

### Peaks/Dips
When the Poisson or Negative Binomial distribution test resulted in a p-value less than 1E-5, a mutational density peak or dip that may span multiple continuous bins is indicated with a light gray rectangle in the mutational density curve plot. A table is shown for each mutational class, which includes p-values and left/right boundary bin indices.  
![Central Dip indicated with gray rectangle](/fig/Dip.JPG)

### Statistical test results
#### Wilcoxon test to compare two mutation density series
In the primary analysis modality, Paired Wilcoxon test is invoked to compare mutational density values between two complementary forms across all bins. Three different ranges are considered: left wing, right wing, and whole range (both wings). In the comparison analysis modalities,  Paired Wilcoxon test is conducted between two cohorts (between-cohort) or two region sets (between-region), rather than between the two mutational forms.  
#### WAVK test to decern non-random spatial trends
The WAVK test is implemented by virtue of the ![funtimes](https://cran.r-project.org/web/packages/funtimes/index.html) package. This test is conducted on one mutation density series only, and the results are shown for the two density series separately. Users can use the WAVK test results as a corroboration for any significant peaks/dips detected in the section above.  




