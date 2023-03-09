######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)


library(estimate)         #引用包
inputFile="merge.txt"     #输入文件名称
setwd("D:\\biowolf\\ICI\\13.estimate")     #设置工作目录

#获取交集基因
filterCommonGenes(input.f=inputFile, 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

#计算免疫得分
estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

#输出每个样品的打分
scores=read.table("estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores=scores[,1:3]
out=rbind(ID=colnames(scores), scores)
write.table(out, file="scores.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
