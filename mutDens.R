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

#' subset2sbsRows(): From mut data, identify two subsetting indix sets for sbs (ix.p) or complementary sbs (ix.m).
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
#' @param points0 zero-points to extend downstream & upstream. Three cols of chr,pos,&strand. Optional 4th col for names
#' @param bsz Bin size
#' @param span Length of adjacent span in one direction to detect peak/dip
getBinsG2 <- function(points0,bsz,span) {
	library(GenomicRanges)
	if (ncol(points0)==3) {
		regName=paste0('region',1:nrow(points0))
		points0$regName=regName
	} 
	colnames(points0)[1:4] <- c('chr','pos','strand','regName')	 	
	origin <- GRanges(
		seqname=points0$chr,
		ranges=IRanges(start=points0$pos,width=1,name=points0$regName),
		strand=points0$strand
	)
	nBins <- span%/%bsz
	G <- GRangesList('0'=origin)
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
 
#' Count mutations falling into left/right adjacent bins
#' Two vectors of Count-like values are returned for a pair of features: coding/template (TSS) or itself/revmut (IS)
#' Output two similar series of summary values: mutCnt (Count) & mutDens (Density). mutDens adjusts for base content in each bin or in whole genome, and mutDens is normalized to Mutation per Mega-Base.
#' @param G GrangeList where each component contains aligned i-th bin in the adjacent span of origins. Three cols of chr,pos,&strand. Optional 4th col for names
#' @param baseCont one of three controlled keywords to designate base content info in vicinity of origins. IS & TSS points to matrix file for -K bins & K bins.
#' @param fTSS Focal feature (positions) denotes TSS or not
mutDensX <- function(mut,G,sbs='C2T',baseCont=c('Origin','TSS','genome')[1],bsz=100,fTSS=FALSE) {
	library(GenomicRanges)
  mut <- GRanges(
    seqnames=mut$chr,
    ranges=IRanges(mut$pos,width=1,names=rownames(mut)),
    strand=mut$strand
  )# GRanges for mutations
	mutCnt_both <- countOverlaps(G,mut,ignore.strand=T)
	mut_strMatch <- countOverlaps(G,mut,ignore.strand=F) #strMatch: stranded match #C2T on right-to-origin & G2A on left-to-origin
	# Two rows of count values are extracted differently for TSS and IS (default)
	if (fTSS) { #This is TSS case
    mutCnt.temp <- mutCnt_both - mut_strMatch
    mutCnt <- rbind(Coding=mut_strMatch,Template=mutCnt.temp)
	} else { #This is IS or default case
		offset <- mutCnt_both-mut_strMatch #C2T on left-to-origin & G2A on right-to-origin
		nT <- length(G); nL <- nT%/%2
		mutCnt.it <- c(offset[1:nL],mut_strMatch[(nL+1):nT]) #C2T
		mutCnt.rev <- c(mut_strMatch[1:nL],offset[(nL+1):nT])#G2A
		mutCnt <- rbind(itself=mutCnt.it,revmut=mutCnt.rev)	
	}
	if (baseCont=='Origin') {
		baseCont <- read.delim('data/atgc2000-100_Origin.tsv',row.names=1)
	} else if (baseCont=='TSS') {
		baseCont <- read.delim('data/atgc2000-100_gTSS.tsv',row.names=1)
	} else {
		baseCont <- c(A=0.591/2,T=0.591/2,G=0.409/2,C=0.409/2) #per 10.1186/s13104-019-4137-z
	}
	base <- substr(sbs,1,1)
	if (is.null(dim(baseCont))) {
		mutable <- length(G[[1]])*bsz*baseCont[base]
		mutDens <- mutCnt/(mutable/1000000)
	} else {
		mutable1 <- length(G[[1]])*bsz*baseCont[names(G),base]
		base2 <- substr(sbs2mut2(sbs)[2],1,1) # C converted to G; T converted to A
		mutable2 <- length(G[[1]])*bsz*baseCont[names(G),base2]
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
#' @param mutDens two-rowed vector of mutational density
#' @param totMut total counts for itself & revmut from VCF file (for mpKm)
#' @return a 3-by-3 matrix cocerning two test methods & three range modalities.
mutDensCmpr <- function(mutDens,totMut=c(1000,1000)) {
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
#' @return 4-by-6 matrix, rows for Dip, Peak, PeakL, PeakR; columns for two sets of p, lIx, rIx
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
#' Obtain sufficient aligned-Bins as background for Poisson-Lamda background estimation
#' Assumes strand of points0 be all +
#' @param foot Offset original point by this distance
#' @param steps Reset original point this many times. steps was fixed to 1 to de facto mute it.
getBinsOut <- function(points0,bsz=100,outerspan=6000,inspan=1000,foot=10,steps=1) { # steps=10
  library(GenomicRanges)
	colnames(points0)[1:3] <- c('chr','pos','strand')
	points.i <- points0
	points.i[,2] <- points.i[,2]+foot # so as to remain unshifted for steps=1
	innerBins <- (outerspan-inspan)%/%bsz
	G <- NULL
	for (i in 1:steps) {#outer loop of steps; actually useless.
		points.i[,2] <- points.i[,2]-foot 
	  dn.j <- GRanges(
  	  seqname=points.i$chr,
    	ranges=IRanges(start=points.i$pos+inspan,width=bsz),
    	strand='+'
  	)
    up.j <- GRanges(
      seqname=points.i$chr,
      ranges=IRanges(start=points.i$pos-inspan,width=bsz),
      strand='+'
    )
		for (j in 1:innerBins) {#innner loop of innerBins
			dn.j <- shift(dn.j,bsz)
			up.j <- shift(up.j,-bsz)
			G <- c(GRangesList(up.j),G,GRangesList(dn.j))
		}
	}
	G	
}

#' Count totals of two sub-classes of mutations for each of six mutational classes
#' @return two-component list of 2-by-6 counts and actual mutations for six classes
classDistrib_wg <- function(vcf,chrCol=3,posCol=4,refCol=5,altCol=6){
	SBSs <- c('C2A','C2G','C2T','T2A','T2C','T2G')
	sbsCnt2 <- matrix(0,nr=2,nc=length(SBSs),dimnames=list(c('itself','revmut'),SBSs)) # two rows of six classes, one row for itself (+) & the other for revmut (-)
	mutS <- vector('list',length=length(SBSs))
	names(mutS) <- SBSs
	for (sbs in SBSs) {
		mut <- mutS[[sbs]] <- extMutSite(vcf,sbs,chrCol=chrCol,posCol=posCol,refCol=refCol,altCol=altCol)
		if (nrow(mut)>0) {
			sbsCnt2[,sbs] <- as.vector(table(mut$strand)[c('+','-')])
		}
	}
	list(cnt2=sbsCnt2,mutS=mutS)
}

#' Count totals of two sub-classes of mutations for each of six mutational classes, in the focal region (2k+2k vicinity of points0)
#' @return 2-by-6 counts for the focal regions (2k+2k vicinity)
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
	SBSs <- c('C2A','C2G','C2T','T2A','T2C','T2G')
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

# plot_sym(): Given 2-by-6 matrix of mutational counts, plot six horizontal bars to signify symmetry situation.
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
 
#' Estimate average value Lamda for Poisson distribution
#' @param G GRanges of outer 100-bp bins of all points
lamdaLocal0 <- function(mut,G) {
	library(GenomicRanges)
  mut <- GRanges(
    seqnames=mut$chr,
    ranges=IRanges(mut$pos,width=1,names=rownames(mut)),
    strand=mut$strand
  )
	mutCnt_both <- countOverlaps(G,mut,ignore.strand=T)
	mut_strMatch <- countOverlaps(G,mut,ignore.strand=F)
	cnts <- list(itself=mut_strMatch,revmut=mutCnt_both-mut_strMatch)
	lamda2 <- c(itself=mean(cnts$itself),revmut=mean(cnts$revmut))
	list(lamda2=lamda2,cnts=cnts)
}
nbParamEstimate <- function(cnts) {
	library(MASS)
	fit <- fitdistr(cnts,'Negative Binomial')
	size<-fit$estimate[1]
	mu<-fit$estimate[2]
	c(size,mu)
}
#' Estimate parameters (lamda; size & mu) for Poisson and Negative Binomial distribution
#' @param G GRanges of outer 100-bp bins of all points
paramLocal <- function(mut,G) {
  library(GenomicRanges)
  mut <- GRanges(
    seqnames=mut$chr,
    ranges=IRanges(mut$pos,width=1,names=rownames(mut)),
    strand=mut$strand
  )
  mutCnt_both <- countOverlaps(G,mut,ignore.strand=T)
  mut_strMatch <- countOverlaps(G,mut,ignore.strand=F)
  cnts <- list(itself=mut_strMatch,revmut=mutCnt_both-mut_strMatch)
	size_mu_1 <- nbParamEstimate(cnts$itself)
	size_mu_2 <- nbParamEstimate(cnts$revmut) 	
  lamda2 <- c(itself=mean(cnts$itself),revmut=mean(cnts$revmut))
	nbparam2 <- rbind(size_mu_1,size_mu_2)
	rownames(nbparam2) <- c('itself','revmut')
	colnames(nbparam2) <- c('size','mu')
  list(lamda2=lamda2,nbparam2=nbparam2,cnts=cnts)
}
#' @param shapeModel background probability model, [pois,nbinom]
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
	peakMat <- matrix(NA,nr=3,nc=3)
	rownames(peakMat) <- c('Peak','PeakL','PeakR')
	colnames(peakMat) <- c('p','ixL','ixR')
	if (peak.yes[midL] | peak.yes[midR]) {
		ixItv <- boundPeak(midL,peak.yes)
		peakP <- min(peak.p[midL:midR]) # Output p & index interval
		peakMat['Peak',] <- c(peakP,ixItv)
	} else {
		peakL.p <- peak.p[1:midL]
		peakL.yes <- peak.yes[1:midL]
    peakR.p <- peak.p[midR:nBin]
    peakR.yes <- peak.yes[midR:nBin]
		if (any(peakL.yes)) {
			lIx <- which.min(peakL.p)
			lItv <- boundPeak(lIx,peakL.yes)
			peakMat['PeakL',] <-c(min(peakL.p),lItv) 
		}
    if (any(peakR.yes)) {
      rIx <- which.min(peakR.p)
      rItv <- boundPeak(rIx,peakR.yes)+midL
      peakMat['PeakR',] <-c(min(peakR.p),rItv)
    }
	}
	peakMat
}
#' Find left & right boundary indices for the affirmed central position
#' @param cIx central indix
#' @param binVec Binary vector indicating if each cell is significant (p<p.th)
#' @return left/right indices for boundary
boundPeak <- function(cIx,binVec) {
	lVec <- binVec[1:cIx]; 
	rVec <- binVec[(cIx+1):length(binVec)]
	ixL <- max(which(!lVec))+1
  ixR <- min(which(!rVec)+cIx)-1
	ixL <- ifelse(is.finite(ixL),ixL,1)
	ixR <- ifelse(is.finite(ixR),ixR,length(binVec))
	ixItv <- c(ixL,ixR)
	ixItv # index interval
}


