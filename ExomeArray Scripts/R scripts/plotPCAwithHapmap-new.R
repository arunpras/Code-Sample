#!/usr/bin/env Rscript
##########################################################################################
# Written by Arun Manoharan
# This script plots the PCA output
# The input is the Eigenstrat output
# This is important to see if the case and control population are homogeneous
# usage: 
# ./plotPCAwithHapmap-new.R  -i <project>.evec -o <projet>.PCA.png
##########################################################################################

#get input arguments
require(getopt)
opt <- getopt(matrix(c('input','i',1,"character",'output','o',1,"character",'datafile','d',1,"character"),ncol=4,byrow=T));

# This reads the population information from hapmap
pop <- read.table("/gnet/is2/research/bioinfo/HumGenet/immuno/immunochip-2011/qc/hapmap/relationships_w_pops_041510.txt",as.is=TRUE);
dat <- read.table(opt$input,skip=1,as.is=T);
names(dat) <- c("ID",paste("EV",1:(ncol(dat)-2),sep="."),"disease");

pop$ID <- paste(pop$V1,pop$V2,sep=":");
dat$pop <- pop$V7[match(dat$ID,pop$ID)];
dat$pop[is.na(dat$pop)] = "Genentech";

# Plot the Eigen vectors output
png(opt$output);
cm=c("violet","blue","blueviolet","cornflowerblue"    ,"yellow","orange","red","magenta","darkgreen",  "brown","cyan","pink");
names(cm)=c("ASW","LWK","MKK","YRI","CHB","CHD","JPT","CEU","TSI","GIH","MEX","Genentech");

with(dat,plot(EV.1,EV.2,col=cm[pop]));
legend("bottomleft",names(cm),col=cm,pch=1);
dev.off();


##########################################################################################
