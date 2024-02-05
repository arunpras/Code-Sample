#!/usr/bin/env Rscript

##########################################################################################
# Written by Arun Manoharan
# This R script identifies the significant Eigen vectors
# The input is the Eigenstrat program output and the output is a list of Eigen vectors
# that are significant
# The script uses generalized linear models using a binomial regression
# usage: 
# ./signf.evec.R -i hapmap.SLEexome.fam -d hapmap.SLEexome.evec -o hapmap.SLEexome.signf.evec
##########################################################################################

require(getopt)
opt <- getopt(matrix(c('input','i',1,"character",'output','o',1,"character",'datafile','d',1,"character"),ncol=4,byrow=T));

fam <- read.table(opt$input,as.is=T);
evec <- read.table(opt$datafile,head=T,as.is=T,comment.char="");
names(evec) <- paste("EV",1:ncol(evec),sep=".");
names(evec)[ncol(evec)] <- "case.control";
evec$disease <- NA;
evec$disease[evec$case.control == "Case"] <- 1;
evec$disease[evec$case.control == "Control"] <- 0;

fmla <- paste("EV",1:(ncol(evec)-2),sep=".");
fmla <- paste(fmla,collapse="+");
fmla <- paste("disease ~ ",fmla);
fmla <- as.formula(fmla);

gm <- glm(fmla,data=evec,family=binomial());
gms <- summary(gm);
sng <- coefficients(gms)[-1,4];
sng <- names(sng[sng<0.05]);
famids <- unlist(lapply(rownames(evec),function(x)unlist(strsplit(x,":"))[[1]]));
indvids <- unlist(lapply(rownames(evec),function(x)unlist(strsplit(x,":"))[[2]]));

covdat <- cbind(famids,indvids,evec[,sng]);
names(covdat)[1:2] <- c("FID","IID");
covdatfile <- paste(opt$output,"plink.covar",sep=".");
sampkeepfile <- paste(opt$output,"keep",sep=".");

write.table(covdat,file=covdatfile,row.names=F,quote=F);
write.table(covdat[,1:2],file=sampkeepfile,row.names=F,col.names=F,quote=F);
##########################################################################################
