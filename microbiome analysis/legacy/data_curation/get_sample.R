#!/usr/bin/Rscript
a<-read.table("a.txt",sep="\t")
b<-read.table("b.txt",sep="\t")
b$V2<-""
b$V2<-a$V2[match(b$V1,a$V1)]
write.table(b$V2,"c.txt",quote=F,row.names=F,col.names=F)
