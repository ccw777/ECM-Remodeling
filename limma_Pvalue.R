######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("impute")

logFoldChange=0.5
adjustP=0.05

library(limma)
library("impute")
setwd("C:\\Users\\c'w'l\\Desktop\\new\\4.diffECMRgene")
rt=read.table("ECMRgeneEXP.txt",sep="\t",header=T)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

#impute missing expression data
mat=impute.knn(exp)
rt=mat$data

rt=avereps(rt)     #基因对应多个探针取均值
#normalize
pdf(file="rawBox.pdf")
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off()
rt=normalizeBetweenArrays(as.matrix(rt))
pdf(file="normalBox.pdf")
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

#rt=log2(rt)             #取log值

#differential
#class <- c("con","con","treat","con","treat","treat")
class <- c(rep("con",347),rep("treat",392))    #需要修改
design <- model.matrix(~0+factor(class))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

#write expression level of diff gene
hmExp=rt[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

#volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$P.Value))
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$P.Value), allDiff$logFC, xlab="-log10(P.Value)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=subset(allDiff, P.Value<adjustP & logFC>logFoldChange)
points(-log10(diffSub$P.Value), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=subset(allDiff, P.Value<adjustP & logFC<(-logFoldChange))
points(-log10(diffSub$P.Value), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()

######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056