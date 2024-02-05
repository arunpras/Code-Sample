##########################################################################################
# Written by Arun Manoharan
# This script plots the GWAS association results
# This script is called from gwasPlots.R and should be in the same directory
#
# p is a vector of p-values in genomic order and chromosome is a vector of chromosomes numbers
# 1:25 for 1:22, X,Y,XY.  You need to edit the first line of the function if you have more
# chromosome indicators or the numbering is different.
##########################################################################################


# You probably also want a better color scheme.  One way is to use the RColorBrewer package

#library(RColorBrewer)
#palette(brewer.pal(n=5,"Set2"))

###############################################################
# Function: generate manhatten plot
###############################################################
palette(c("tan4","darkblue"));
stripe<-function(p,chromosome,...){
    colv <- rep(c("lightblue","lavender"),100);
    chl <- as.numeric(names(table(chromosome)));
    chromcode <- c(1:22,"X","XY","Y","M")
    chromcode <- chromcode[chl];

    # convert p-values to friendly form
    logp <- -log(p,10)

    # N is the number of SNPS, ymax is its log
    N<-length(logp)
    ymax<-log(N,10)+4

    # find the first index in chromosome fo each chromosome number
    chromstart <- which(c(1,diff(chromosome))>=1)

    # last index in chromosome fo each chromosome number (not really, its the next chromosome)
    chromend   <-c(chromstart[-1],N)

    x <- (1:N) + chromosome * (chromend[1]/6)
    y <- pmin(ymax,logp)

    plot(x,y,cex=0.3+1.5*y/ymax, col=chromosome, pch=ifelse(y==ymax,24,19),
           bg=chromosome,xlab="Chromosome",ylab=quote(-log[10](p-value)),
           xaxt="n",ylim=c(0,ymax),axes=F,cex.lab=1.5,...)

     
    centers<- (x[chromstart]+x[chromend]-(chromend[1]/6))/2
    centers[length(centers)]<-x[N]

    stopifnot(length(centers)==length(chromcode))
    axis(1,at=centers,label=chromcode,cex.axis=1.1,font.axis=2);
    axis(2,cex.axis=1.2,font.axis=2);
    }

###############################################################
# Function: QQ plot generation
###############################################################
plot.qq <- function(z, neglog10=TRUE, ...){
	m <- length(z)
	if(neglog10)	plot( x=-log10((1:m)/m), y = -log10(sort(z)), ... )
	else  plot( x=(1:m)/m, y = sort(z), ... )
}
points.qq <- function(z, neglog10=TRUE, ...){
	m <- length(z)
	if(neglog10)	points( x=-log10((1:m)/m), y = -log10(sort(z)), ... )
	else  points( x=(1:m)/m, y = sort(z), ... )
}
###############################################################


