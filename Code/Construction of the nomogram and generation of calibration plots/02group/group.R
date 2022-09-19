install.packages("caret")
install.packages("survival") 


library(survival)
library(caret)
dir="C:\\Users\\11649\\Desktop\\DBT\\列线图\\nomg\\02group"
setwd(dir)
TCGA<-read.table("input.txt",header=T,sep="\t",check.names = F,stringsAsFactors = F)
set.seed(300)
data<-createDataPartition(y=TCGA$ID,p=0.50,list=F)
trian_data<-TCGA[data, ]
test_data<-TCGA[-data,] 
write.table(trian_data, "trian_data.txt",row.names = F,quote = F,sep = "\t")
write.table(test_data, "test_data.txt",row.names = F,quote = F,sep = "\t")


