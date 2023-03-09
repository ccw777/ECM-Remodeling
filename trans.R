######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)             #���ð�
inputFile="symbol.txt"     #�����ļ�
setwd("D:\\biowolf\\ICI\\07.FPKM2TPM")     #���ù���Ŀ¼

#��ȡ�����ļ�,���������ļ�����
outTab=data.frame()
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)
data=t(data)

#FPKMת��ΪTPM
fpkmToTpm=function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(data, 2, fpkmToTpm)

#���ת�����
tpmOut=rbind(ID=colnames(tpm), tpm)
write.table(tpmOut, file="TCGA.TPM.txt", sep="\t", col.names=F, quote=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056