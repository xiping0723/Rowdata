dir="C:\\Users\\11649\\Desktop\\DBT\\列线图\\nomg\\03mult"
setwd(dir)


library(survival)


TCGA<-read.table("trian_data.txt",header=T,sep="\t") 


TCGA$Age<-factor(TCGA$Age,labels=c("<=60",">60"))
TCGA$T.satge<-factor(TCGA$T.satge,labels=c("T1&T2","T3&T4"))
TCGA$M.stage<-factor(TCGA$M.stage,labels=c("M0","M1"))
TCGA$N.stage<-factor(TCGA$N.stage,labels=c("N0","N1"))
TCGA$Pathologic.grade<-factor(TCGA$Pathologic.grade,labels=c("Stage I&Stage II","Stage III&Stage IV"))
TCGA$Histologic.grade<-factor(TCGA$Histologic.grade,labels=c("G1&G2","G3&G4"))
TCGA$DBT<-factor(TCGA$DBT,labels=c("high","low"))



ddist <- datadist(TCGA)
options(datadist='ddist')


cox2 <- cph(Surv(survival_time,status) ~Age+T.satge+N.stage+M.stage+Pathologic.grade
            +Histologic.grade+DBT,surv=T,x=T, y=T,data=TCGA)


surv <- Survival(cox2)
sur_1_year<-function(x)surv(1*365*1,lp=x)
sur_3_year<-function(x)surv(1*365*3,lp=x)
sur_5_year<-function(x)surv(1*365*5,lp=x)

nom_sur<- nomogram(cox2,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,
                   funlabel=c('1-Year Survival','3-Year Survival','5-Year survival'),
                   maxscale=100,
                   fun.at=c('0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))

pdf("nom.pdf",13,8)
plot(nom_sur)
dev.off()

nom_sur

fmla1 <- as.formula(Surv(survival_time,status) ~Age+T.satge+N.stage+M.stage+Pathologic.grade
                    +Histologic.grade+DBT)
cox3 <- coxph(fmla1,data=TCGA)
summary(cox3)



cox4 <- cph(Surv(survival_time,status) ~ Age+T.satge+N.stage+M.stage+Pathologic.grade
            +Histologic.grade+DBT,surv=T,x=T, y=T,time.inc = 1*365*1,data=TCGA)
cal <- calibrate(cox4, cmethod="KM", method="boot", u=1*365*1, m= 40, B=1000)

pdf("calibrate1.pdf")
plot(cal,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of 1-Year Survival",ylab="Actual 1-Year Survival",col="blue",sub=F)
mtext("")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")
dev.off()


cox5 <- cph(Surv(survival_time,status) ~ Age+T.satge+N.stage+M.stage+Pathologic.grade
            +Histologic.grade+DBT,surv=T,x=T, y=T,time.inc = 1*365*3,data=TCGA)
cal <- calibrate(cox5, cmethod="KM", method="boot", u=1*365*3, m= 40, B=1000)

pdf("calibrate3.pdf")
plot(cal,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of 3-Year Survival",ylab="Actual 3-Year Survival",col="blue",sub=F)
mtext("")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")
dev.off()



cox6 <- cph(Surv(survival_time,status) ~ Age+T.satge+N.stage+M.stage+Pathologic.grade
            +Histologic.grade+DBT,surv=T,x=T, y=T,time.inc = 1*365*5,data=TCGA)
cal <- calibrate(cox6, cmethod="KM", method="boot", u=1*365*5, m= 40, B=1000)

pdf("calibrate5.pdf")
plot(cal,lwd=2,lty=1,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram-Predicted Probability of 5-Year Survival",ylab="Actual 5-Year Survival",col="blue",sub=F)
mtext("")
box(lwd = 0.5)
abline(0,1,lty = 3,lwd = 2,col = "black")
dev.off()  
