#!/usr/bin/env Rscript
##########################################################################################
# Written by Arun Manoharan
# This R script identifies the samples and SNPs passing the missingness cutoffs
# The input is the plink missingness command output and the output is a list of samples 
# and SNPs passing the input cutoffs
# It also plots the SNP and sample missingness
# usage: 
# ./missingness.R -i <project>.missing.imiss -l <project>.missing.lmiss -s 0.05 -p 0.1 -o <project>.ms
##########################################################################################

#get input arguments
require(getopt)
opt <- getopt(c(
	'imiss','i',1,"character",
	'lmiss','l',1,"character",
	'snpcut','s',1,"double",
	'sampcut','p',1,"double",
	'out','o',1,"character"
	)
    );

print(opt$imiss);
print(opt$lmiss);
imiss <- read.table(opt$imiss,head=TRUE);
lmiss <- read.table(opt$lmiss,head=TRUE);

print(dim(imiss));

imiss$color[imiss$F_MISS > opt$sampcut] <- "red"
imiss$color[imiss$F_MISS < opt$sampcut] <- "blue"
lmiss$color[lmiss$F_MISS > opt$snpcut] <- "red"
lmiss$color[lmiss$F_MISS < opt$snpcut] <- "blue"

sampfail<-length(imiss$color[imiss$F_MISS > opt$sampcut])
snpfail<-length(lmiss$color[lmiss$F_MISS > opt$snpcut])

#plot sample missingness
png(paste(opt$out,".samp.png",sep=""));
with(imiss,plot(sort(F_MISS),xlab="Sample number",ylab="Missing genotype rate",main=paste("Missingess in samples",sampfail, "failed",sep=" ")));
abline(h=opt$sampcut);
dev.off();

# plot SNP missingness
png(paste(opt$out,".snp.png",sep=""));
with(lmiss,plot(sort(F_MISS),xlab="Snp number",ylab="Missing genotype rate",main=paste("Missingess in snps",snpfail, "failed",sep=" ")));
abline(h=opt$snpcut);
dev.off();

#identify samples passing cutoff
kept.samp <- subset(imiss,F_MISS < opt$sampcut);
samp.file <- paste(opt$out,".samp.keep.txt",sep="");
write.table(kept.samp,file=samp.file,quote=F,row.names=F,col.names=F);

#identify samples failing cutoff
rm.samp <- subset(imiss,F_MISS >= opt$sampcut);
rmsamp.file <- paste(opt$out,".samp.rm.txt",sep="");
write.table(rm.samp,file=rmsamp.file,quote=F,row.names=F);

print("wrote samp file");

#identify SNPs passing cutoff
kept.snp <- subset(lmiss,F_MISS < opt$snpcut);
snp.file <- paste(opt$out,".snp.keep.txt",sep="");
write.table(kept.snp,file=snp.file,quote=F,row.names=F,col.names=F);
print("wrote snp file");

summary.file <- paste(opt$out,".missing.summary.txt",sep="");

cat(nrow(imiss) - nrow(kept.samp)," samples removed.\n",file=summary.file,append=T);
cat(nrow(lmiss) - nrow(kept.snp)," snps removed.\n",file=summary.file,append=T);
print("wrote summary file");
##########################################################################################
