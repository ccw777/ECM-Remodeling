######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("survival")
#install.packages("survminer")


#���ð�
library(survival)
library(survminer)
setwd("C:\\Users\\c'w'l\\Desktop\\ECMR\\new\\3.ECMRcluster")     #���ù���Ŀ¼

#�����������ߵĺ���
bioSurvival=function(inputFile=null, outFile=null){
  #��ȡ�����ļ�
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #�Ƚϸߵͷ�����������죬�õ�������pֵ
  diff=survdiff(Surv(futime, fustat) ~ EAGcluster, data=rt)
  pValue=1-pchisq(diff$chisq, df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ EAGcluster, data = rt)
  #print(surv_median(fit))
  
  #������������
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     surv.median.line = "hv",
                     legend.title="EAGcluster",
                     legend.labs=c("EAGcluster1", "EAGcluster2"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("#e64b35", "#4dbbd5"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile, onefile=FALSE, width=6.5, height=5.5)
  print(surPlot)
  dev.off()
}

#���ú�����������������
bioSurvival(inputFile="risk.txt", outFile="survival.pdf")


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056