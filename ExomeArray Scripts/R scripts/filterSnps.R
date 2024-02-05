#!/usr/bin/env Rscript
##########################################################################################
# Written by Arun Manoharan
# This script filters bad SNPs and generates list of clean SNPs
# The script takes in as input: plink hardy command output, plink freq command output and 
#    plink test-missing command output
# The script also outputs the list of clean SNPs
# usage: 
# ./filterSnps.R -h <project>).EUR.hwe.hwe -f <project>.EUR.related-removed.maf.frq -t <project>.EUR.test_missing.missing -m 0.00 -o <project>.EUR.clean.snps -d 1e-5 -c 1e-5 -n <project>.missing.lmiss -r 0.05
##########################################################################################

#get input arguments
require(getopt)
opt <- getopt(c(
	'hardy','h',1,"character",
	'freq','f',1,"character",
	'testmiss','t',1,"character",
	'out','o',1,"character"
	)
    );

print(opt);

# cutoffs to be used
hardy.cutoff <- 1e-4;
freq.cutoff  <- 0.01;
testmiss.cutoff <- 1e-4;

#filter SNPs failing HWE
hpipe <- paste("cat ",opt$hardy," | grep -e TEST -e UNAFF",sep="");
hardy <- read.table(pipe(hpipe),head=TRUE,as.is=TRUE);
hardy.filt <- subset(hardy,TEST=="UNAFF" & P > hardy.cutoff);
hardy.fail <-subset(hardy,TEST=="UNAFF" & P <= hardy.cutoff);
paste("Hardy fail:",length(hardy.fail$V1))
dim(hardy.fail)

#filter SNPs failing MAF
freq <- read.table(opt$freq,head=TRUE,as.is=TRUE);
freq.filt <- subset(freq,MAF > freq.cutoff);
freq.fail <- subset(freq,MAF <= freq.cutoff);
paste("MAF fail:",length(freq.fail$V1))
dim(freq.fail)

#filter SNPs failing differential missingness
testmiss <- read.table(opt$testmiss,head=TRUE,as.is=TRUE);
testmiss.filt <- subset(testmiss,P > testmiss.cutoff);
testmiss.fail <- subset(testmiss, P <= testmiss.cutoff);
paste("diff missingness fail:",length(testmiss.fail$V1))
dim(testmiss.fail)

#write out the clean SNP list
clean.snps  <- intersect(intersect(hardy.filt$SNP,freq.filt$SNP),testmiss.filt$SNP);
cat(clean.snps,file=opt$out,sep="\n");

##########################################################################################



