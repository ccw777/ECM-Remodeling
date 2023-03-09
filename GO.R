######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05     #p值过滤条件
qvalueFilter=1        #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("D:\\biowolf\\ICI\\28.GO")           #设置工作目录
rt=read.table("featureGenes.exp.txt", header=T, sep="\t", check.names=F)     #读取输入文件
allGenes=unique(as.vector(rt[,1]))
geneType=gsub("(.*?)\\|(.*?)", "\\2", allGenes)
allGenes=gsub("(.*?)\\|.*", "\\1", allGenes)

for(i in levels(factor(geneType))){
	#基因名字转换为基因id
	genes=allGenes[geneType==i]
	entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
	entrezIDs=as.character(entrezIDs)
	gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
	
	#GO富集分析
	kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
	GO=as.data.frame(kk)
	GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
	#保存富集结果
	write.table(GO, file=paste0(i, ".GO.txt"), sep="\t", quote=F, row.names = F)
	
	#定义显示Term数目
	showNum=10
	if(nrow(GO)<30){
		showNum=nrow(GO)
	}
	
	#柱状图
	pdf(file=paste0(i, ".barplot.pdf"), width=9, height=7)
	bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
	print(bar)
	dev.off()
			
	#气泡图
	pdf(file=paste0(i, ".bubble.pdf"), width=9, height=7)
	bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
	print(bub)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
