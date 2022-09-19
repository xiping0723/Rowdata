dir="C:\\Users\\11649\\Desktop\\DBT\\列线图\\nomg\\06risk test"
setwd(dir)


library(survival)


TCGA<-read.table("test_data.txt",header=T,sep="\t")

TCGA$Age<-factor(TCGA$Age,labels=c("<=60",">60"))
TCGA$T.satge<-factor(TCGA$T.satge,labels=c("T1&T2","T3&T4"))
TCGA$M.stage<-factor(TCGA$M.stage,labels=c("M0","M1"))
TCGA$N.stage<-factor(TCGA$N.stage,labels=c("N0","N1"))
TCGA$Pathologic.grade<-factor(TCGA$Pathologic.grade,labels=c("Stage I&Stage II","Stage III&Stage IV"))
TCGA$Histologic.grade<-factor(TCGA$Histologic.grade,labels=c("G1&G2","G3&G4"))
TCGA$DBT<-factor(TCGA$DBT,labels=c("high","low"))


cox2 <- coxph(Surv(survival_time,status) ~ Age+T.satge+N.stage+M.stage+Pathologic.grade
              +Histologic.grade+DBT,data=TCGA)

risk_score<-predict(cox2,type="risk",newdata=TCGA)
risk_level<-as.vector(ifelse(risk_score>median(risk_score),"High","Low"))
write.table(cbind(id=rownames(cbind(TCGA[,1:2],risk_score,risk_level)),cbind(TCGA[,1:2],risk_score,risk_level)),"risk_score.txt",sep="\t",quote=F,row.names=F)




library(survival)
library(timeROC)
TCGA1<-read.table("risk_score.txt",header=T,sep="\t")
predict_1_year<- 1*365
predict_3_year<- 3*365 
predict_5_year<- 5*365 

ROC<-timeROC(T=TCGA1$survival_time,delta=TCGA1$status,
             marker=TCGA1$risk_score,cause=1,
             weighting="marginal",
             times=c(predict_1_year,predict_3_year,predict_5_year),ROC=TRUE)

pdf("ROC.pdf")
plot(ROC,time=predict_1_year,title=FALSE,lwd=3)
plot(ROC,time=predict_3_year,col="blue",add=TRUE,title=FALSE,lwd=3)
plot(ROC,time=predict_5_year,col="green",add=TRUE,title=FALSE,lwd=3)
legend("bottomright",
       c(paste("AUC of 1 year survival: ",round(ROC$AUC[1],3)),
         paste("AUC of 3 year survival: ",round(ROC$AUC[2],3)),
         paste("AUC of 5 year survival: ",round(ROC$AUC[3],3))),col=c("red","blue","green"),lwd=3)
dev.off()
