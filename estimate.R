######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)


library(estimate)         #���ð�
inputFile="merge.txt"     #�����ļ�����
setwd("D:\\biowolf\\ICI\\13.estimate")     #���ù���Ŀ¼

#��ȡ��������
filterCommonGenes(input.f=inputFile, 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

#�������ߵ÷�
estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

#���ÿ����Ʒ�Ĵ��
scores=read.table("estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
scores=scores[,1:3]
out=rbind(ID=colnames(scores), scores)
write.table(out, file="scores.txt", sep="\t", quote=F, col.names=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056