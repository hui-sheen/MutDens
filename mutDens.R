### Assign (NA,NA) to nbparam & two rows of (NA,NA) to nbparam & nbparam2 when options$shapeModel is pois
source('mutDens_baseCont.R')
#gFeatures<- c('TSS','TES','IS','RIP','PAS','Junc')
gns<- c('hg38','mm39','rn7','galGal6','danRer11','dm6','sacCer3','rheMac10','canFam3')#,'hg19','mm9') #Genome scope
#' Extract mutation sites for a pair of mutation types (e.g. C>T & G>A)
#' @return Point mutation sites specified in three columns: chr,pos,strand (+, C2T; -,G2A)
extMutSite <- function(vcf,sbs='C2T',chrCol=1,posCol=2,refCol=4,altCol=5,smpCol=NULL,head=FALSE) {
	vcf <- read.delim(vcf,header=head,as.is=T,comment.char='#')
  mut <- vcf[,c(chrCol,posCol,refCol,altCol)]
  colnames(mut) <- c('chr','pos','ref','alt')
  chrs <- c('X','Y',1:22)
  mut$chr <- gsub('chr','',mut$chr)
  mut <- subset(mut,chr%in%chrs)
  mut <- unique(mut)
	if (!is.null(smpCol)) {
		if (is.numeric(smpCol)) {
			GT0=vcf[,9+smpCol]
		} else {
			GT0=vcf[,smpCol]
		}
		GT <- sapply(strsplit(GT0,':'),function(x) x[1])
		mut <- subset(mut,GT%in%c('1/1','1/0','0/1')) #assumes first-ever element in sample column is 1/1 or the alike.
	}
  ix2 <- subset2sbsRows(mut,sbs)
  if ( any(ix2$ix.p) & any(ix2$ix.n) ) {
    mut.p <- data.frame(unique(mut[ix2$ix.p,c('chr','pos')]),strand='+')
    mut.n <- data.frame(unique(mut[ix2$ix.n,c('chr','pos')]),strand='-')
    mut <- rbind(mut.p,mut.n)
    mut$strand <- as.character(mut$strand)
  } else {
    mut <- matrix(NA,nr=0,nc=3) #Return null mut matrix if not seeing any mutation of the required type.
  }
  mut
}

#' subset2sbsRows(): From mut data, identify two subsetting indix sets for sbs (ix.p) and its complementary sbs (ix.n).
#' @param mut four-columned matrix, headers set as chr, pos, ref, alt.
#' @return list of two components for mutation form 1 (+) and form 2 (-).
#' OUTPUT ix.p: row indices for interested mutations that are on forward strand.
#' OUTPUT ix.n: row indices for interested mutations that are on reverse strand.
subset2sbsRows <- function(mut,sbs='C2T') {
  sbs <- toupper(unlist(strsplit(sbs,'2')))
  src <- sbs[1]
  tgt <- sbs[2]
  maptable <- c('A','T','G','C')
  names(maptable) <- c('T','A','C','G')
  rev.src <- maptable[src]
  rev.tgt <- maptable[tgt]
  ix.p <- (toupper(mut$ref)==src&toupper(mut$alt)==tgt)
  ix.n <- (toupper(mut$ref)==rev.src&toupper(mut$alt)==rev.tgt)
  list(ix.p=ix.p,ix.n=ix.n)
}

#' Formalize name C2T to paired names C>T & G>A
#' @param sbs Informal representation of SBS (e.g. C2T)
#' @return Two character strings for paired SBS forms (e.g., C>T & G>A)
sbs2mut2 <- function(sbs) {
  sbs <- toupper(unlist(strsplit(sbs,'2')))
  src <- sbs[1]
  tgt <- sbs[2]
  maptable <- c('A','T','G','C')
  names(maptable) <- c('T','A','C','G')
  rev.src <- maptable[src]
  rev.tgt <- maptable[tgt]
	mut <- paste(src,tgt,sep='>')
	revmut <- paste(rev.src,rev.tgt,sep='>')
	c(mut,revmut)
}

#' Plot two parallel lines for two coupled binwise mutational density series
#' @param mutDens two-row matrix of mutational densities for template/coding strands or two complementary forms respectively
#' @param bandIx Boundary bin indices for peak/dip bands. Matrix of two columns (left & right). Setting to NULL means no peak/dip or not drawing
plot2mutDens <- function(mutDens,totMut=NULL,
	lnames=c('C>T','G>A'),bsz=100,bandIx=NULL,title='Mutation Density around TSS',cols=c('orange','blue'),legend=T) 
{
  if (!is.null(totMut) ) {
    mutDens <- mutDens/(totMut/1000)
		ylab='MPKM'
  } else {
		ylab='Mutations per Mb'
	}
	par(mar=c(10,11,4,2),xpd=TRUE)
  xcoords <- as.numeric(colnames(mutDens))*bsz
	ymin <- min(mutDens)
	ymax <- max(mutDens)+(max(mutDens)-min(mutDens))/10
  plot(xcoords,mutDens[2,],ylim=c(ymin,ymax),bty='l',type='n',xlab='',ylab='',cex.axis=1.5,las=1,cex.main=3,main=title)
  mtext('Relative distance to origin (bp)',side=1,line=3,cex=1.5)
  mtext(ylab,side=2,line=4.5,cex=1.5)
  if (!is.null(bandIx)) {
    L1.bounds <- bandIx[!is.na(bandIx[,1]),1:2]
    L2.bounds <- bandIx[!is.na(bandIx[,3]),3:4]
    bounds <- rbind(L1.bounds,L2.bounds)
    for (i in 1:nrow(bounds)) { #Add gray band for peak/dip
      rect(xcoords[bounds[i,1]],ymin,xcoords[bounds[i,2]],ymax,col=rgb(211/255,211/255,211/255,alpha=0.5),border=NA)
    }
  }
  for (i in 1:nrow(mutDens)) {
    lines(xcoords,mutDens[i,],type='l',col=cols[i],lwd=4)
  }
	if (legend)
  	legend('topleft',col=cols,legend=lnames,lwd=4,cex=1,bty='n',inset=c(-0.2,-0.4))
}

#' Plot four parallel lines in two sets of orange & blue pairs for two cohorts or two regions
#' @param mutDens1 two-row matrix of parallel mutational densities for cohort/region 1
#' @param mutDens1 two-row matrix of parallel mutational densities for cohort/region 2
plot4mutDens <- function(mutDens1,mutDens2,totMut1=NULL,totMut2=NULL,
	lnames=c('C>T (C1)','G>A (C1)','C>T (C2)','G>A (C2)'),bsz=100,
	title='Mutation Density around TSS',cols=c('orange','blue'),ltys=c('solid','dotted'),legend=T) 
{
	if (!is.null(totMut1)&!is.null(totMut2) ) {
		mutDens1 <- mutDens1/(totMut1/1000)
		mutDens2 <- mutDens2/(totMut2/1000) 
		ylab='MPKM'
	} else {
		ylab='Mutations per Mb'
	}
  par(mar=c(10,11,4,2),xpd=TRUE)
  xcoords <- as.numeric(colnames(mutDens1))*bsz
	mutDens <- rbind(mutDens1,mutDens2)
  ymin <- min(mutDens)
  ymax <- max(mutDens)+(max(mutDens)-min(mutDens))/10
  plot(xcoords,mutDens[2,],ylim=c(ymin,ymax),bty='l',type='n',xlab='',ylab='',cex.axis=1.5,las=1,cex.main=3,main=title)
  mtext('Relative distance to origin (bp)',side=1,line=3,cex=1.5)
  mtext(ylab,side=2,line=4,cex=1.5)
	cols <- rep(cols,2)
	ltys <- rep(ltys,each=2)
  for (i in 1:nrow(mutDens)) {
    lines(xcoords,mutDens[i,],type='l',col=cols[i],lty=ltys[i],lwd=4)
  }
  if (legend)
    legend('topleft',col=cols,lty=ltys,legend=lnames,lwd=4,cex=1,bty='n',inset=c(-0.2,-0.4))
}

#' Generate a number of equal-sized bins in bi-directional vicinity of points0.
#' 
#' Depends on GenomicRanges
#' @param points0 zero-points to extend downstream & upstream. Two or Three cols of chr,pos,&strand. Optional 4th col for names
#' @param bsz Bin size
#' @param span Length of adjacent span in one direction to detect peak/dip
getBinsG2 <- function(points0,bsz,span) {
	library(GenomicRanges)
	if (!any(grepl('^str',colnames(points0)))) {
		pre3<- data.frame(points0[,1:2],strand='+')
		if (ncol(points0)>=3) {
			points0<- data.frame(pre3,points0[,3:ncol(points0)])
		} else {
			points0<- pre3
		}
	}
	colnames(points0)[1:3] <- c('chr','pos','strand') 	 	
	origin <- GRanges(
		seqname=points0$chr,
		ranges=IRanges(start=points0$pos,width=1),#name=points0$regName),
		strand=points0$strand
	)
	nBins <- span%/%bsz
	G <- GRangesList('0'=origin)
	#Initial dn/up.i blocks were located at the opposite realms, ready to move to the expected directions
	# First stride of dn.i block covers the 1st position (origin GR, as above)
  dn.i <- resize(shift_bystr(origin,-1),width=bsz,fix='end',ignore.strand=FALSE)
	up.i <- resize(origin,width=bsz,fix='start',ignore.strand=FALSE)
  for (i in 1:nBins) {
    BSZ <- ifelse(i!=nBins,bsz,(bsz+span%%bsz))
    dn.i <- shift_bystr(dn.i,BSZ)
		up.i <- shift_bystr(up.i,-BSZ)
    G <- c(GRangesList(up.i),G,GRangesList(dn.i))
		names(G)[1] <- -i
    names(G)[length(G)] <- i
  }
  if (!any(strand(origin)=='-')) {
    for (i in (-nBins):(-1)) {
			strand(G[[as.character(i)]]) <- '-'
		} #This is for IS case - setting all upstream bins to minus strand
  }
	G <- G[setdiff(names(G),'0')]
	G
}
 
#' Get base content data, either static four prortions or matrix of bins of four proportions  
#' Accommodating select genomes: hg38/19, mm10/9, rn6, galGal5,danRer10, dm6, sacCer3, rheMac8,canFam3
#' Pre-built GCcont matrices for select gFeatures (mostly for humans): TSS (other gns), TES (other gns), IS, RIP, PAS, Junc
#' The default output is a 4-value vector of rep(0.25,4)
#' Invokes baseCont() and series of subfunctions from separate file baseCont.R
#' @param gFt One of controlled key words for select genomic features, for which base content matrix has been prebuilt [gFt static/genome for genome-wide static ones]
#' @param points0 Points0 will induce definition of calcSrc
#' @param calcSrc Calculation source; either NULL, focal position table, or formatted base cont (4-tuple vector or binwise matrix)
getBaseCont<- function(span=2000,bin=100,gn='mm10',gFt='custom',calcSrc=NULL, 
	gns=gns,prefix=paste0('atgc',span,'-',bin,'_',gFt,'_',gn)) {
	baseCont.res<- static4<- switch(gn,
		hg38=c(A=0.591/2,T=0.591/2,G=0.409/2,C=0.409/2),
		hg19=c(A=0.591/2,T=0.591/2,G=0.409/2,C=0.409/2),
		rep(0.25,4)
	)
	libName<- paste0('data/atgc2000-100_',gFt,'_',gn,'.tsv')
	errMsg_unknownGn<- 'For unknown genome, non-null readily available baseCont vec/matrix is expected and passed by calcSrc!'
	wrnMsg_staticVct<- 'Genome-wide static base content (possibly four 0.25s) is used!'
	wrnMsg_readyCont<- 'Existant base cont in a file name is supplied by user!'
	wrnMsg_calcFly<- 'Bin-wise base content will be calculated on the fly!'
	if (gn%in%gns) {
		if ( span==2000&bin==100 & ( file.exists(libName) ) ) { #(gFt%in%gFeatures) 
      baseCont.res<- read.delim(libName,row.names=1) #load baseCont matrix for gFt for hg genome (or TES/TSS for all gns)
		} else if ( !is.null(calcSrc) & !tolower(gFt)%in%c('static','genome') & !is.null(dim(calcSrc)) ){
			warning(wrnMsg_calcFly) #Calc on the fly will be evoked;	
	  	calcSrc<- regPosTbl(calcSrc,gn) #Ensure chrs are regularized for smoothy getSeq()
    	baseCont.res<- baseCont(calcSrc,span,bin,gn,prefix=prefix) #paste0('atgc',span,'-',bin,'_',gFt,'_',gn,'.tsv'))#$prop
		} else {
			warning(wrnMsg_staticVct)
		}
	} else {
		warning('#TODO: what is direction on rare genome???')
	}
	baseCont.res
}

#' Count mutations falling into left/right adjacent bins
#' Two vectors of Count-like values are returned for a pair of features: coding/template (TSS) or itself/revmut (IS)
#' Output two similar series of summary values: mutCnt (Count) & mutDens (Density). mutDens adjusts for base content in each bin or in whole genome, and mutDens is normalized to Mutation per Mega-Base.
#' @param points0 DF read in from focal position file (will be converted to GRanges G in one commandline)
#' @param gFt Key word for genomic feature - IS,TSS,TES,etc. points to matrix file for -K bins & K bins; custom means calculating baseCont on the fly for any focal vicinities
#' @param calcSrc can be NULL or 4-tuple or points0 matrix (input to getBinsG())
mutDensX <- function(mut,points0,sbs='C2T',bsz=100,span=2000,gFt='custom',gn='hg38',
	calcSrc=NULL) {
	library(GenomicRanges)
	G <- getBinsG2(points0,as.numeric(bsz),as.numeric(span)) 
  mut <- GRanges(
    seqnames=mut$chr,
    ranges=IRanges(mut$pos,width=1,names=rownames(mut)),
    strand=mut$strand
  )# GRanges for mutations
	mutCnt_both <- countOverlaps(G,mut,ignore.strand=T)
	mut_strMatch <- countOverlaps(G,mut,ignore.strand=F) #strMatch: stranded match #C2T on right-to-origin & G2A on left-to-origin
	# Two rows of count values are extracted differently for TSS and IS (default)
	fTSS<- ifelse(gFt%in%c('TSS','proTSS','lncTSS'),TRUE,FALSE)
	if (fTSS) { #This is TSS case
    mutCnt.temp <- mutCnt_both - mut_strMatch #template strand
    mutCnt <- rbind(Coding=mut_strMatch,Template=mutCnt.temp)
	} else { #This is IS or default case
		offset <- mutCnt_both-mut_strMatch #C2T on left-to-origin & G2A on right-to-origin
		nT <- length(G); nL <- nT%/%2
		mutCnt.it <- c(offset[1:nL],mut_strMatch[(nL+1):nT]) #C2T
		mutCnt.rev <- c(mut_strMatch[1:nL],offset[(nL+1):nT])#G2A
		mutCnt <- rbind(itself=mutCnt.it,revmut=mutCnt.rev)	
	}
	if (is.null(calcSrc)) calcSrc<- points0
	baseCont.res<- getBaseCont(span,bsz,gn,gFt=gFt,calcSrc=calcSrc,gns=gns) #gFeatures & gns are defined on very top as global variable
	base <- substr(sbs,1,1)
	if (is.null(dim(baseCont.res))) {
		mutable <- length(G[[1]])*bsz*baseCont.res[base]
		mutDens <- mutCnt/(mutable/1000000)
	} else {
		mutable1 <- length(G[[1]])*bsz*baseCont.res[names(G),base]
		base2 <- substr(sbs2mut2(sbs)[2],1,1) # C converted to G; T converted to A
		mutable2 <- length(G[[1]])*bsz*baseCont.res[names(G),base2]
		mutable <- rbind(mutable1,mutable2)
		mutDens <- mutCnt/(mutable/1000000) #keep output as 2-rowed vector still
	}	
	list(mutDens=mutDens,mutCnt=mutCnt) 
}

#' Shift GRanges x-bp in the downstream direction, considering strand of each GRange
shift_bystr <- function(gr,x) {
  library('GenomicRanges',quietly=T,warn.conflicts=F)
  str <- as.vector(strand(gr))
	rev.ix <- str=='-'
	step = rep(x,length(gr))
	step[rev.ix] = -x
	gr <- shift(gr,step)
}
#' Statistically test the difference between two rows of mutational densities.
#' totMut uses default 1000, meaning technically no further normalization to MPKM is conducted.
#' @param mutDens two-rowed vector of mutational density
#' @param totMut total counts for itself & revmut from VCF file (for mpKm)
#' @return a 4-by-3 matrix cocerning three range modalities for three tests & meanDiff.
mutDensCmpr <- function(mutDens,totMut=c(1000,1000)) {
  mutDens <- mutDens/(totMut/1000) # Normalize to MPKM
  rtv <- matrix(NA,nr=4,nc=3)
  rownames(rtv) <- c('meanDiff','Wilcoxon','WAVKtrend1','WAVKtrend2') #diff is [,1]-[,2] #WAVKtrend tests L1 & L2 separately
  colnames(rtv) <- c('whole','left','right')
  if (any(mutDens!=0)) {
    mutDens.lst <- list(whole=mutDens,left=mutDens[,1:ncol(mutDens)%/%2],right=mutDens[,(ncol(mutDens)%/%2+1):ncol(mutDens)])
    for (item in names(mutDens.lst)) {
      mutDens.i <- mutDens.lst[[item]]
      rtv['Wilcoxon',item] <- wilcox.test(mutDens.i[1,],mutDens.i[2,],paired=TRUE,alternative='two.sided')$p.val
      rtv[c('WAVKtrend1','WAVKtrend2'),item] <- tryCatch(  {
				apply(mutDens.i,1,
					function(x) notrend_test(x,test='WAVK',factor.length = 'adaptive.selection',q=0.9)$p.value ) },
				error=function(e) {rtv[c('WAVKtrend1','WAVKtrend2'),item]}  )
      rtv['meanDiff',item] <- mean(mutDens.i[1,]-mutDens.i[2,])
    }
  }
  return(rtv)
}

#' Statistically test the difference between two rows of mutational densities.
#' totMut uses default 1000, meaning technically no further normalization to MPKM is conducted.
#' @param mutDens two-rowed vector of mutational density
#' @param totMut total counts for itself & revmut from VCF file (for mpKm)
#' @return a 3-by-3 matrix cocerning two test methods & three range modalities.
mutDensCmpr0 <- function(mutDens,totMut=c(1000,1000)) {
  mutDens <- mutDens/(totMut/1000) # Normalize to MPKM
  rtv <- matrix(NA,nr=3,nc=3)
  rownames(rtv) <- c('Wilcoxon','KS','meanDiff') #diff is [,1]-[,2]
  colnames(rtv) <- c('whole','left','right')
  if (any(mutDens!=0)) {
    rtv <- matrix(NA,nr=3,nc=3)
    rownames(rtv) <- c('Wilcoxon','KS','meanDiff') #diff is [,1]-[,2]
    colnames(rtv) <- c('whole','left','right')
    mutDens.lst <- list(whole=mutDens,left=mutDens[,1:ncol(mutDens)%/%2],right=mutDens[,(ncol(mutDens)%/%2+1):ncol(mutDens)])
    for (item in names(mutDens.lst)) {
      mutDens.i <- mutDens.lst[[item]]
      rtv['Wilcoxon',item] <- wilcox.test(mutDens.i[1,],mutDens.i[2,],paired=TRUE,alternative='two.sided')$p.val
      rtv['KS',item] <- ks.test(mutDens.i[1,],mutDens.i[2,],alternative='two.sided')$p.val
      rtv['meanDiff',item] <- mean(mutDens.i[1,]-mutDens.i[2,])
    }
  }
  return(rtv)
}

#' Takes in two-rowed bin-wise Cnt values & backgrould lamda (2) and detect peak(s) & dip
#' @param shapeModel background probability model, [pois,nbinom]
#' @return 4-by-6 matrix, rows for Peak, PeakL, PeakR, Dip; columns for two sets of p, lIx, rIx
mutDensShape <- function(mutCnt,lamda2,nbparam2,pTh=1e-5,shapeModel) {
	cnt1 <- mutCnt[1,]
	cnt2 <- mutCnt[2,]
	peakMat1 <- judgePeak(cnt1,lamda2[1],nbparam2[1,],pTh,shapeModel)
	peakMat2 <- judgePeak(cnt2,lamda2[2],nbparam2[2,],pTh,shapeModel)
	colnames(peakMat1) <- paste('L1',colnames(peakMat1),sep='_')
	colnames(peakMat2) <- paste('L2',colnames(peakMat2),sep='_')
	dip1 <- judgeDip(cnt1,lamda2[1],nbparam2[1,],pTh,shapeModel)
	dip2 <- judgeDip(cnt2,lamda2[2],nbparam2[2,],pTh,shapeModel)
	res <- cbind(peakMat1,peakMat2)
	res <- rbind(res,matrix(c(dip1,dip2),nr=1))
	rownames(res)[4] <- 'Dip'
	res			
}

#' estimate average count in local vicinity (10K minus 2K in one direction)
#' Removed two arguments foot & steps
#' Obtain sufficient aligned-Bins as background for Poisson-Lamda background estimation
#' Ignoring innate strand of points0; working on + strand only
getBinsOut <- function(points0,bsz=100,outerspan=7000,inspan=2000) {
  library(GenomicRanges,quietly=T,warn.conflicts=F)
  colnames(points0)[1:2] <- c('chr','pos')#,'strand')
  innerBins <- (outerspan-inspan)%/%bsz
  G <- NULL
	dn.j <- GRanges(
		seqname=points0$chr,
		ranges=IRanges(start=points0$pos+inspan,width=bsz),
		strand='+'
	)
	up.j <- GRanges(
		seqname=points0$chr,
		ranges=IRanges(start=points0$pos-inspan,width=bsz),
		strand='+'
	)#start=points0$pos-inspan, so that foreground left-flank is ovoided
	for (j in 1:innerBins) {#innner loop of innerBins
		dn.j <- shift(dn.j,bsz)
		up.j <- shift(up.j,-bsz)
		G <- c(GRangesList(up.j),G,GRangesList(dn.j))
	}
  G
}
#' Expand out all 16 [Di-Nuc] context types for triNuc with central ignored
#' @return res Vector of 16 triNuc specifications (e.g.A*A or A*C)
contextTypes<- function(totBase=3) {
	if (totBase==3) {
		base4<- c('A','C','G','T')
		lBase16<- rep(base4,each=4)
		rBase16<- rep(base4,4)
		res<- paste0(lBase16,'N',rBase16)
	} else {
		warning('Have not considered contexts other than triNuc! NULL type is returned!!!')
		res<- NULL
	}
	res	
}

#' Extract mutation sites for a pair of mutation types (e.g. C>T & G>A)
#' @return Point mutation sites specified in three columns: chr,pos,strand (+, C2T; -,G2A)
#' @param mut0 Data frame to define the chromosome coordinate of central sbs (of triNuc context). Must include chr & pos
#' @param gn BSgenome object for specific species (most typically, BSgenome.Hsapiens.UCSC.hg38)
#' @param adjBase Number of adjacent/neighboring bases to consider in uni-direction
#' @return Mutation data frame with one more field ($context) characterizing triNuc type
contextTyping<- function( mut0,gn,adjBase=1,chrs=c(1:22,c('X','Y')) ) {
  mut<- mut0
  if (any(!grepl('chr',mut$chr))) {
    mut<- subset(mut,chr%in%chrs)
    mut[,1] <- paste0('chr',mut[,1])
  } else {
    padChr<- paste0('chr',chrs)
    mut<- subset(mut,chr%in%padChr)
  }
  ctrPos<- mut$pos
  st0<- ctrPos-adjBase
  ed0<- ctrPos+adjBase
  allowable<- seqlengths(gn)[mut$chr]
  startPos<- ifelse(st0>allowable-1,allowable,st0)
  endPos<- ifelse(ed0>allowable,allowable,ed0)
  context<- toupper(getSeq(gn,mut$chr,startPos,endPos,as.character=T))
  ctrNuc<- toupper( substr(context,adjBase+1,adjBase+1) )
  suspicious<- ctrNuc!=mut$ref
  substr(context,adjBase+1,adjBase+1)<- 'N' # [TriNuc] Context type, e.g., A*A or A*C
  context.df<- data.frame(chr=gsub('chr','',mut$chr),pos=mut$pos,context=context)
  if (any(suspicious)) {
    warning(sum(suspicious),'out of',length(suspicious),'mutation bases discordant with getSeq() central were discarded!!!')
    context.df<- context.df[!suspicious,]
  }
  mut<- merge(mut0,context.df)
  mut
}

#' Count totals of two sub-classes of mutations for each of six mutational classes
#'
#' @return two-component list of 2-by-6 (or 2-by-96 for triCont is T) counts and actual mutations for six classes
classDistrib_wg <- function(vcf,gn,chrCol=3,posCol=4,refCol=5,altCol=6,triCont=FALSE){ #wg: whole genome
	SBS6 <- c('C2A','C2G','C2T','T2A','T2C','T2G')
	dinucs <- contextTypes(totBase=3) #Expand out all 16 dinucs (occupying 1st & 3rd bases of triNuc)
	if (!triCont) {
		triSBSs<- SBS6
	} else {
		dinucs <- contextTypes(totBase=3)  
		SBSs<- rep(SBS6,each=length(dinucs)) #Each of all 6 sbs' is expanded into 16 repeats
		triSBSs<- paste(SBSs,dinucs,sep='_') # "C2T_A*C"
	}
	sbsCnt2 <- matrix(0,nr=2,nc=length(triSBSs),dimnames=list(c('itself','revmut'),triSBSs)) # two rows of six classes, one row for itself (+) & the other for revmut (-)
	mutS <- vector('list',length=length(triSBSs))
	names(mutS) <- triSBSs
	for (sbs in SBS6) {
		mut<- extMutSite(vcf,sbs,chrCol=chrCol,posCol=posCol,refCol=refCol,altCol=altCol)
		if (!triCont) {
			mutS[[sbs]]<- mut
			if (nrow(mut)>0) {
				sbsCnt2[,sbs] <- as.vector(table(mut$strand)[c('+','-')])
			}
		} else {
			superMut.typed<- contextTyping(mut,gn=gn,adjBase=1)
			superMut.split<- split(superMut.typed[,setdiff(colnames(superMut.typed),'context')],superMut.typed$context)
			for (dinuc in dinucs) { #Loop over 16 types
				trisbs<- paste(sbs,dinuc,sep='_')
				mutS[[trisbs]]<- mut<- superMut.split[[dinuc]]
				if (nrow(mut)>0) {
					sbsCnt2[,trisbs] <- as.vector(table(mut$strand)[c('+','-')])
				}
			}
		}
	}
	list(cnt2=sbsCnt2,mutS=mutS)
}

#' Count totals of two sub-classes of mutations for each of six mutational classes, in the focal region (2k+2k vicinity of points0)
#'
#' points0's strands are uniformly coerced to all +'s, while mutS's strands for it&revmut are retained.
#' @param span (foregraound) uni-directional span (left-flank or right-flank)
#' @param mutS List of six components of mutation classes ($mutS output of classDistrib_wg())
#' @return 2-by-6 (or 2-by-96 for triCont) counts for the focal regions (2k+2k vicinity)
classDistrib_fr <- function(mutS,points0,span=2000) { # _fr focal_region
	library(GenomicRanges)
	points0 <- read.delim(points0,as.is=T,head=T)
	colnames(points0)[1:2] <- c('chr','pos')
  origin <- GRanges(
    seqname=points0$chr,
    ranges=IRanges(start=points0$pos,width=1),
    strand='+' # Regardless of TSS or IS
  )
	span2 <- resize(origin,width=span,fix='end',ignore.strand=TRUE)
	span2 <- resize(span2,width=2*span,fix='start',ignore.strand=TRUE) #Do not care about detailed 1-bp endpoint
	SBSs <- names(mutS) #c('C2A','C2G','C2T','T2A','T2C','T2G')
	sbsCnt2 <- matrix(0,nr=2,nc=length(SBSs),dimnames=list(c('itself','revmut'),SBSs))
	for (sbs in SBSs) {
		mut <- mutS[[sbs]]
		mut <- GRanges(
			seqnames=mut$chr,
			ranges=IRanges(mut$pos,width=1,names=rownames(mut)),
			strand=mut$strand
		)
		cnt_both <- sum(countOverlaps(mut,span2,ignore.strand=TRUE))
		cnt_strMatch <- sum(countOverlaps(mut,span2,ignore.strand=FALSE))
		sbsCnt2[,sbs] <- c(cnt_strMatch, cnt_both-cnt_strMatch)
	}
	sbsCnt2
}

# plot_sym(): Plot symmetry situations for six mutation classes (six horizontal bars to signify symmetry situation).
#' @param sym 2-by-6 matrix of mutational counts, $cnts2 output of classDistrib_wg()
plot_sym <- function(sym,sample='') {
	colnames(sym) <- gsub('2','>',colnames(sym))
  par(mar=c(6,4,4,10),cex.axis=1.5,cex.main=1.5)
  sym.normed <- t(t(sym)/colSums(sym))
  coord=barplot(sym.normed[,ncol(sym.normed):1],horiz=T,las=1,main=paste(sample,'\nSymmetry'),names.arg=rep('',ncol(sym)),
		legend.text=T, args.legend=list(x='right',pt.cex=2,cex=1.5,bty='n',inset=c(-0.4,0)))
	cols <- rev(rainbow(ncol(sym)))
	for (i in 1:length(cols)) {
		axis(2,at=coord[i],labels=rev(colnames(sym))[i],col.axis=cols[i],las=1)
	}	
  abline(v=0.5,lty='dashed',col='red',lwd=2)
		
}

plot_pieANDsym<- function(cnt2rows,title='Whole-genome') {
	pie(colSums(cnt2rows),labels=gsub('2','>',colnames(cnt2rows)),col=rainbow(ncol(cnt2rows)),main=paste(title,'proportion',sep='\n'))
	plot_sym(cnt2rows,title)
}

#' Estimate parameters for Negative Bayes distribution
#' Drafted by Chung-I on 3/17/22
nbParamEstimate <- function(cnts) {
	library(MASS)
	fit <- fitdistr(cnts,'Negative Binomial',method='SANN')
	size<-fit$estimate[1]
	mu<-fit$estimate[2]
	c(size,mu)
}
#' Estimate parameters (lamda; size & mu) for Poisson and Negative Binomial distribution
#' @param G GRanges of outer 100-bp bins of all points
#' @return Three components: lamda2 of two scalars, nbparam2 two rows of two scalars each, and cnts of two rows of background cnts
paramLocal <- function(mut,G,shapeModel) {
  library(GenomicRanges)
  mut <- GRanges(
    seqnames=mut$chr,
    ranges=IRanges(mut$pos,width=1,names=rownames(mut)),
    strand=mut$strand
  )
  mutCnt_both <- countOverlaps(G,mut,ignore.strand=T)
  mut_strMatch <- countOverlaps(G,mut,ignore.strand=F)
  cnts <- list(itself=mut_strMatch,revmut=mutCnt_both-mut_strMatch)
	if (shapeModel=='nbinom') {
		size_mu_1 <- nbParamEstimate(cnts$itself)
		size_mu_2 <- nbParamEstimate(cnts$revmut) 	
	} else {
		size_mu_1<- size_mu_2<- c(NA,NA)
	}
  lamda2 <- c(itself=mean(cnts$itself),revmut=mean(cnts$revmut))
	nbparam2 <- rbind(size_mu_1,size_mu_2)
	rownames(nbparam2) <- c('itself','revmut')
	colnames(nbparam2) <- c('size','mu')
  list(lamda2=lamda2,nbparam2=nbparam2,cnts=cnts)
}

#' @param shapeModel background probability model, [pois,nbinom]
#' @return Vector of three scalars, 1st for minP, 2nd & 3rd for left/right dip boundary indices
judgeDip <- function(mutCnt1,lamda,nbparam,pTh=1e-5,shapeModel='pois') {
  nBin <- length(mutCnt1)
  midL <- nBin%/%2; midR <- midL+1
	if (shapeModel=='pois') {
  	dip.p <- ppois(mutCnt1,lamda,lower.tail=TRUE) #test for Dip rather than Peak
	} else {
		dip.p <- pnbinom(mutCnt1,size=nbparam[1],mu=nbparam[2],lower.tail=TRUE)
	}
  dip.yes <- dip.p<=pTh
  if (dip.yes[midL] | dip.yes[midR]) {
    ixItv <- boundPeak(midL,dip.yes)
    p <- min(dip.p[ixItv])
    res <- c(p,ixItv)
  } else {
    res <- rep(NA,3)
  }
  res
}
#' @param shapeModel background probability model, [pois,nbinom]
#' @return 3-by-3 matrix with rows for Peak, PeakL, &PeakR, and columns for minP and left/right peak boundary indices
#' Identify focal peak bin and merge adjacent ones into a final peak
#' Operate on one party of the mutation pair (so it is called mutCnt1)
judgePeak <- function(mutCnt1,lamda,nbparam,pTh=1e-5,shapeModel='pois') {
	nBin <- length(mutCnt1)
	midL <- nBin%/%2; midR <- midL+1
	if (shapeModel=='pois') {
		peak.p <- ppois(mutCnt1,lamda,lower.tail=FALSE)
	} else {
		peak.p <- pnbinom(mutCnt1,size=nbparam[1],mu=nbparam[2],lower.tail=FALSE)
	}
	peak.yes <- peak.p<=pTh
	peakMat <- matrix(NA,nr=3,nc=3) #Ultimate output of 3-by-3 matrix
	rownames(peakMat) <- c('Peak','PeakL','PeakR')
	colnames(peakMat) <- c('p','ixL','ixR')
	if (peak.yes[midL] | peak.yes[midR]) {
		ixItv <- boundPeak(midL,peak.yes)
		peakP <- min(peak.p[midL:midR])
		peakMat['Peak',] <- c(peakP,ixItv)
	} else {
		peakL.p <- peak.p[1:midL]
		peakL.yes <- peak.yes[1:midL]
    peakR.p <- peak.p[midR:nBin]
    peakR.yes <- peak.yes[midR:nBin]
		if (any(peakL.yes)) { #Left-flank has a peak root
			lIx <- which.min(peakL.p)
			lItv <- boundPeak(lIx,peakL.yes)
			peakMat['PeakL',] <-c(min(peakL.p),lItv) 
		}
    if (any(peakR.yes)) { #Right-flank has a peak root
      rIx <- which.min(peakR.p)
      rItv <- boundPeak(rIx,peakR.yes)+midL
      peakMat['PeakR',] <-c(min(peakR.p),rItv)
    }
	}
	peakMat
}
#' Find left & right boundary indices for the affirmed central position
#' @param cIx central indix or, more generally, one index of the binVec of significance
#' @param binVec Binary vector indicating if each cell is significant (p<p.th)
#' @return left/right indices for boundary
boundPeak <- function(cIx,binVec) {
	lVec <- binVec[1:cIx]; 
	rVec <- binVec[(cIx+1):length(binVec)]
	ixL <- max(which(!lVec))+1
  ixR <- min(which(!rVec)+cIx)-1
	ixL <- ifelse(is.finite(ixL),ixL,1) #finite occurs when whole left-flank are significant
	ixR <- ifelse(is.finite(ixR),ixR,length(binVec)) #finite occurs when whole right-flank are significant
	ixItv <- c(ixL,ixR)
	ixItv # index interval
}

flatten6<- function(cnt96) { #cnt96: Two rows of 96 triNuc classes' counts
	fct6cls<- gsub('_.*','',colnames(cnt96))
	cnt96.split<- split(data.frame(t(cnt96)),fct6cls) #split works well with explicit df
	cnt6<- sapply(cnt96.split,colSums)
}

