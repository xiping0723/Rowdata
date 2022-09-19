library("limma")

setwd("C:\\Users\\11649\\Desktop\\DBT\\列线图\\nomg\\01diff")              
gene="DBT"                                                          
normalNum=72                                                        
tumorNum=535                                                  

rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)           
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)


uniq=rbind(ID=colnames(data),data)
write.table(uniq,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)       


Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))
single=cbind(ID=colnames(data),expression=data[gene,],Type)
colnames(single)=c("ID",gene,"Type")
write.table(single,file="singleGene.txt",sep="\t",quote=F,row.names=F)
