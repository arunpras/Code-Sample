#!/usr/bin/env Rscript

##########################################################################################
# Written by Arun Manoharan
# This script plots the heterozygosity and filters the individuals with 
# >3SD heterozygosity
# The input is the plink het command output and the output is a list of samples 
# passing the  cutoffs
# It also plots the heterozygosity
# usage: 
# ./heterozygosity.R  -h <project>.EUR.het -s 4 -o <project>.EUR
##########################################################################################


#get input arguments
require(getopt)
opt <- getopt(c(
	'het','h',1,"character",
	'cutoff','s',1,"character",
	'out','o',1,"character"
	)
    );

print(opt);

het <- read.table(opt$het, head=TRUE)

het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM.";
Hlo = with(het, mean(HET_RATE)-sd(HET_RATE)*opt$cutoff)
Hhi = with(het, mean(HET_RATE)+sd(HET_RATE)*opt$cutoff)
het$color[het$HET_RATE<Hlo | het$HET_RATE>Hhi] <- "red"
het$color[het$HET_RATE>Hlo & het$HET_RATE<Hhi] <- "blue"

sampfail<-length(het$color[het$HET_RATE<Hlo | het$HET_RATE>Hhi])

#plot the heterozygosity histogram
png(filename=paste(opt$out,".het.hist.png",sep=""))
hist(het$HET_RATE, main = "Histogram of Heterozygosity", xlab = "Heterozygosity")
abline(v = Hlo, lty = 2, col = "red")
abline(v = Hhi, lty = 2, col = "red")
abline(v = mean(het$HET_RATE), lty = 2, col = "black")
dev.off()

#plot the heterozygosity frequency
png(filename=paste(opt$out,".het.plot.png",sep=""))
plot(het$HET_RATE,het$F, main=paste( "Heterozygosity",sampfail, "failed",sep=" "), xlab="Heterozygosity", ylab="F",col=het$color)
abline(v = Hlo, lty = 2, col = "red")
abline(v = Hhi, lty = 2, col = "red")
abline(v = mean(het$HET_RATE), lty = 2, col = "black")
dev.off()

# write list of heterozygosity passing and failing samples
write.table(het, file=paste(opt$out,".het.out.txt",sep=""), row.name=F, quote=F, col.name=T)
drop = het[which(het$HET_RATE<Hlo | het$HET_RATE>Hhi), ]
write.table(paste(drop$FID, drop$IID, sep = " "), file=paste(opt$out,".het.samp.rm.txt",sep = ""), row.name = F, quote = F, col.name = F)

##########################################################################################
