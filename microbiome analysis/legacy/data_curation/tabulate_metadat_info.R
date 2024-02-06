#######################################################################################################
# Written by Arun Manoharan
# This script pulls the stats for the american gut project

#######################################################################################################
setwd("~/Desktop/other_projects/qiime2/datasets/baseline/")
a<-read.table("10317_20170816-124134.txt",sep="\t",head=T,quote="")
dim(a)
x<-as.data.frame(c(""))
dim(table(a[,193]))

write.table(x,"summary.txt",sep="\t",col.names=F,row.names=F,quote=F)
for(i in 2:dim(a)[2])
{
  if(i!=27 && i!=56 && i!=100 && i !=193)
  {
    b<-as.data.frame(table(a[,i]))
    names(b)[1]<-"values"
    b$criteria<-colnames(a)[i]
    b<-b[order(b$Freq),]
    print(i)
    print(dim(b))
    write.table(b,"summary.txt",sep="\t",col.names=T,row.names=F,quote=F,append=T)
    write.table(x,"summary.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    
  }
}