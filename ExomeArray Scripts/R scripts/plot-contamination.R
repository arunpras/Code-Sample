#!/usr/bin/env Rscript
##########################################################################################
# Written by Arun Manoharan
# This script plots the sample related pairs
# The input is the plink genome command output
# usage: 
# ./plot-contamination.R -i <project>.test.rel.cnt -o <project>.EUR.contamination
##########################################################################################


#get input arguments
require(getopt)
opt <- getopt(matrix(c('input','i',1,"character",'output','o',1,"character"),ncol=4,byrow=T));
infile <- read.table(opt$input,as.is=T);
outfile <- paste(opt$output,".png",sep="");
outfile1 <- paste(opt$output,"-hist.png",sep="");

# Plot related pairs
png(outfile)
plot(infile$V2,xlab="samples",ylab="No of related samples",main="plot of relatedness", col="blue",cex=0.4)
abline(h=30,lty=2,col="red")
dev.off()
png(outfile1)
hist(infile$V2[infile$V2<100],xlab="No of related samples", ylab="No of samples",main="",col="lightblue")
dev.off()
##########################################################################################
