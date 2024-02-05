#!/usr/bin/env Rscript
##########################################################################################
# Written by Arun Manoharan
# This script plots the sample relatedness
# The plot identifies siblings, parent-offspring,twins and duplicates
# The input is the plink genome command output
# The script also outputs the sample pairs which are related
# usage:
# ./ibd.R -i <project>.EUR.genome.rel -o <project>.EUR.ibdplot -v 0.4 -m <project>.missing.imiss
##########################################################################################

#get input arguments
require(getopt)
opt <- getopt(c(
'ibd','i',1,"character",
'out','o',1,"character",
'vline','v',1,"double",
'hline','h',1,"double",
'imiss','m',1,"character"
)
);

print(opt);
rel <- read.table(opt$ibd,head=TRUE);

#plot relatedness
png(paste(opt$out,".png",sep=""));
with(rel,plot(Z0,Z1,xlim=c(0,1),ylim=c(0,1),xlab="k0",ylab="k1"));
# this plots the reference points
points(c(0,0,0.25,0.5,0.75),c(0,1,0.5,0.5,0.25),pch=16,col=c("magenta","cyan","red","yellow","green"),cex=0.75);

# plot legends
legend(0.175,1,c("Reported unrelated"),col=c("black"),pch=c(1),bty="n");
legend(0.60,1,c("monozy twins/duplicate","parent-offspring","full sibs","half sibs, grandparent, aunt","first cousin"),title="Expected values",pch=16,col=c("magenta","cyan","red","yellow","green"),bty="n");

if(exists("vline",opt))
abline(v=opt$vline,lty=2,col="red")
if(exists("hline",opt))
abline(h=opt$hline,lty=2,col="red")
dev.off();

relpairs <- subset(rel,Z0 < opt$vline);


imiss <- read.table(opt$imiss,head=TRUE);
relpairs$id1 <- with(relpairs,paste(FID1,IID1,sep=" "));
relpairs$id2 <- with(relpairs,paste(FID2,IID2,sep=" "));
imiss$id <- with(imiss,paste(FID,IID,sep=" "));

relpairs$m1 <- imiss$F_MISS[match(relpairs$id1,imiss$id)];
relpairs$m2 <- imiss$F_MISS[match(relpairs$id2,imiss$id)];

rmid <- character();

if(dim(relpairs)[1]>0)
{
    relpairs$rm1 <- 0;
    relpairs$rm2 <- 0;
    
    for(i in 1:nrow(relpairs))
    {
        if(relpairs$rm1[i]==0 & relpairs$rm2[i]==0){
            mid <- relpairs$id1[i];
            if(relpairs$m1[i] < relpairs$m2[i])
            mid <- relpairs$id2[i];
            relpairs$rm1[relpairs$id1==mid] <- 1;
            relpairs$rm2[relpairs$id2==mid] <- 1;
            rmid <- c(rmid,mid);
        }
    }
    
}



#write output files
write.table(relpairs,file = paste(opt$out,".related.genome",sep=""),quote=F,row.names=F);
cat(rmid,file = paste(opt$out,".related.removed",sep=""),sep="\n");



##########################################################################################

