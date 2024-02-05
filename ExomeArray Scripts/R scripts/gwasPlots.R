#!/usr/bin/env Rscript
##########################################################################################
# Written by Arun Manoharan
# This script plots the GWAS association results
# The input is the plink logistic regression command output
# The script outputs the manhatten and QQ plots
# usage: 
# ./gwasPlots.R -i <project>.EUR.logistic.assoc.logistic -o <project>.EUR.gwasplot
##########################################################################################

#get input arguments
require(getopt)
opt <- getopt(matrix(c(
	'input','i',1,"character",
	'output','o',1,"character"
	),ncol=4,byrow=TRUE)
    );

# call gwtestplots_ms.r script
source("gwtestplots_ms.r");
pcom <- paste("cat ",opt$input," | grep -e TEST -e ADD",sep="");
res <- read.table(pipe(pcom),head=TRUE,as.is=TRUE);

manf <- paste(opt$output,".manhatten.png",sep="");
qqf <- paste(opt$output,".qq.png",sep="");

print(manf);

#plot manhatten
png(manf,width=1200,height=700);
par(mar=c(5.1,5.1,3.1,2.1));
pvdat=res[,c("P","CHR")];
pvdat=na.omit(pvdat);
stripe(pvdat[,1], pvdat[,2],cex.main=1.5);
abline(h=-log10(5e-8),col="gray32");
dev.off();

# plot QQ
png(qqf,height=700,width=700);
plot.qq(pvdat[,1],neglog10=TRUE,xlab="Expected -log10(p-value)",
    ylab="Observed -log10(p-value)",cex.lab=1.5,cex.axis=2);
abline(0,1,col="red",lwd=2);
dev.off() ;


##########################################################################################

