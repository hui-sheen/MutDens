#' Load the corresponding BSgenome obeject in light of gn string
#' @param gn0 String for genome, eg. hg38, mm10
#' @return The corresponding BSgenome object
gnStr2Obj<- function(gn0) {
	gnPkg<- switch(gn0,
		hg38='BSgenome.Hsapiens.UCSC.hg38',
		hg19='BSgenome.Hsapiens.UCSC.hg19',
		mm39='BSgenome.Mmusculus.UCSC.mm39',
		danRer10='BSgenome.Drerio.UCSC.danRer10',
		canFam3='BSgenome.Cfamiliaris.UCSC.canFam3',
		galGal6='BSgenome.Ggallus.UCSC.galGal6',
		rheMac10='BSgenome.Mmulatta.UCSC.rheMac10',
		rn7='BSgenome.Rnorvegicus.UCSC.rn7',
		danRer11='BSgenome.Drerio.UCSC.danRer11',
		dm6='BSgenome.Dmelanogaster.UCSC.dm6',
		sacCer3='BSgenome.Scerevisiae.UCSC.sacCer3',
	)
	if (!require(gnPkg,character.only=T,quietly=T)) BiocManager::install(gnPkg)
	library(gnPkg,character.only=T,quietly=T)
	cmdStr<- paste0('gn<-',gnPkg)
	eval(parse(text=cmdStr))
	gn
}
#' Regularize position table
#' Chromosomes were ensured to be canonical ones. Outputted chr has prefix string of chr
#' @param bed file name for focal positions (or a DF that was just read in)
#' @return DF of 3 cols: chr,pos,str
regPosTbl<- function(bed,gn0) {
	chrs<- switch(gn0,
		hg38=c(1:22,c('X','Y')),
		hg19=c(1:22,c('X','Y')),
		mm39=c(1:19,c('X','Y')),
		rheMac10=c(1:20,c('X','Y')),
		rn7=c(1:20,c('X','Y')),
		galGal6=c(1:33,c('W','Z')),
		canFam3=c(1:38,'X'),
		dm6=c('X','Y','2L','2R','3L','3R','4'),
		sacCer3=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI'),
		danRer11=c(1:25)
	)
	if (is.null(dim(bed))) {
		tbl<- read.delim(bed,header=T,as.is=T)[,1:3]
	} else {
		tbl<- bed
	}
	colnames(tbl)<- c('chr','pos','str')
	if (any(!grepl('chr',tbl$chr))) {
		tbl<- subset(tbl,chr%in%chrs)
		tbl[,1] <- paste0('chr',tbl[,1])
	}	else {
		padChr<- paste0('chr',chrs)
		tbl<- subset(tbl,chr%in%padChr)
	}
	tbl
}

getBinsG <- function(origin,bsz,span) {
  library(GenomicRanges)
  nBins <- span%/%bsz
  G <- GRangesList('0'=origin)
  dn.i <- resize(shift_bystr(origin,-1),width=bsz,fix='end',ignore.strand=FALSE) #1st downstream window uses points0 as an endpoint
  up.i <- resize(origin,width=bsz,fix='start',ignore.strand=FALSE)
  for (i in 1:nBins) {
    BSZ <- ifelse(i!=nBins,bsz,(bsz+span%%bsz))
    dn.i <- shift_bystr(dn.i,BSZ)
    up.i <- shift_bystr(up.i,-BSZ)
    G <- c(GRangesList(up.i),G,GRangesList(dn.i))
    names(G)[1] <- -i
    names(G)[length(G)] <- i
  }
  G <- G[setdiff(names(G),'0')]
  G
}

#' Calculate numbers and proportions of four nucs of a nuc segment
#'
#' @param seg Segment (Vector) of A/T/G/C letters
#' @return Vector of eight values, first 4 proportions, & second 4 frequency numbers
atgcStat <- function(seg) {
	seg<- toupper(seg) #Coercing all letters to upper cases
  unknown <- toupper(seg)=='N'
  seg <- seg[!unknown]
  An <- sum(seg=='A')
  Tn <- sum(seg=='T')
  Gn <- sum(seg=='G')
  Cn <- sum(seg=='C')
  vec1 <- c(An,Tn,Gn,Cn)
  vec2 <- vec1/length(seg)
  vec <- c(vec1,vec2)
  names(vec) <- paste0(rep(c('A','T','G','C'),2),rep(c('num','freq'),each=4))
  vec  
}
#' Summarize all GRanges into A/T/C/G percentages for overall G (one conceptual Bin).
#' @param gn BSgenome object for specific species
#' @return Vector of 8 values, first 4 numbers and second 4 percentages
sumGrInG <- function(gr,gn) {
  st0 <- start(gr)
  st0[st0<=1] <- 1
	ed0<- end(gr)
	chrNames<- as.character(seqnames(gr))
	allowable<- seqlengths(gn)[chrNames] 
	st<- ifelse(st0>allowable-1,allowable,st0)
	ed<- ifelse(ed0>allowable,allowable,ed0)
  Seqs <- getSeq(gn,seqnames(gr),st,ed,as.character=T)
  Nucs <- strsplit(Seqs,'') #Converting nuc strings into segments (letter vectors)
  grRes <- sapply(Nucs,atgcStat) #Matrix of 8 rows for A/T/G/C stats, with cols for segments
  ATGC <- rowSums(grRes[1:4,]) #Total number of each nuc types in all segments
  atgc <- ATGC/sum(ATGC)
  c(ATGC,round(atgc,3))
}

#' Generate conceptual Grand Bins and loop over Bins to obtain A/T/C/G statistics (8-valued vectors) in each Bin
#' @param gr Origin GRanges
#' @param gn BSgenome object for specific species
#' @return First component $prop has exact 8 cols required as a MutDens baseCont table
loopOverBins <- function(gr,span,bin,gn) {
  G <- getBinsG(gr,bin,span)
  nBins <- span%/%bin
  FourFreq <- matrix(0,nr=nBins*2,nc=8)  
  for (i in 1:(2*nBins)) {
    FourFreq[i,] <- sumGrInG(G[[i]],gn)
  }
  colnames(FourFreq) <- paste0(rep(c('A','T','G','C'),2),rep(c('num',''),each=4)) #Omitting suffix freq for percentage columns
  rownames(FourFreq) <- names(G)
  FourProp<- FourFreq[,5:8]
	list(prop=FourProp,freq=FourFreq)
}
#' baseCont() differs from atgc_binG() in that tbl specifies a point with or without strand.
#' @param gn BSgenome object for specific species
baseCont <- function(tbl,span,bsz,gn,prefix='atgcComp/hg38_dextIS6k') {
  nBins <- span%/%bsz
  if (is.element('str',colnames(tbl))) {
    gr <- GRanges(
      seqnames=tbl$chr,
      ranges=IRanges(start=tbl$pos,width=1), #field named changed to pos
      strand=tbl$str #field named changed to str
    ) } else {
    gr <- GRanges(
      seqnames=tbl$chr,
      ranges=IRanges(start=tbl$pos,width=1), #field named changed to pos
      strand='+'
    )      
  }
	gnObj<- gnStr2Obj(gn)
  Res0 <- loopOverBins(gr,span,bsz,gnObj) 
  Res <- Res0$prop
  write.table(Res,paste0(prefix,'.tsv'),col.names=NA,sep='\t',quote=F)
	tiff(paste0(prefix,'.tif'))
		cols=c('lightblue','darkblue','pink','magenta')
		pches=c(16,4,16,4)
		par(mar=c(5,10,6,2),xpd=T)
		plot(as.numeric(rownames(Res))*bsz,Res[,'A'],type='p',bty='L',
			col=cols[1],pch=pches[1],ylim=c(0.19,0.31),las=1, cex.lab=2,cex.axis=1.5,
			xlab='Distance from origin',ylab='')
		mtext(text='Base content in bins',side=2,line=4,cex=2)
		points(as.numeric(rownames(Res))*bsz,Res[,'T'],col=cols[2],pch=pches[2])#cex=2
		points(as.numeric(rownames(Res))*bsz,Res[,'G'],col=cols[3],pch=pches[3]) #,cex=2
		points(as.numeric(rownames(Res))*bsz,Res[,'C'],col=cols[4],pch=pches[4])# ,cex=2
		legend('top',legend=c('A','T','G','C'),horiz=T,pch=pches,col=cols,cex=2,pt.cex=2,
			bty='n',inset=c(0,-0.15))
	dev.off()
  Res0$prop # Same as the written file
}
