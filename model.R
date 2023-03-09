######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("glmnet")
#install.packages("survival")
#install.packages('survminer')


#���ð�
library(glmnet)
library(survival)
library(survminer)
inputFile="uniSigExp.txt"        #��������������ı��������ļ�
setwd("C:\\Users\\c'w'l\\Desktop\\new\\8.model")       #���ù���Ŀ¼
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #��ȡ�����ļ�

#COXģ�͹���
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox, direction="both")
multiCoxSum=summary(multiCox)

#���ģ�������Ϣ
outMultiTab=data.frame()
outMultiTab=cbind(
		          coef=multiCoxSum$coefficients[,"coef"],
		          HR=multiCoxSum$conf.int[,"exp(coef)"],
		          HR.95L=multiCoxSum$conf.int[,"lower .95"],
		          HR.95H=multiCoxSum$conf.int[,"upper .95"],
		          pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
write.table(outMultiTab, file="multiCox.txt", sep="\t", row.names=F, quote=F)

#��������ļ�
score=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)
outCol=c("futime", "fustat", coxGene)
risk=as.vector(ifelse(score>median(score), "high", "low"))
outTab=cbind(rt[,outCol], riskScore=as.vector(score), risk)
write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056