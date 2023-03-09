######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05     #pֵ��������
qvalueFilter=1        #�������pֵ��������

#������ɫ
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("D:\\biowolf\\ICI\\28.GO")           #���ù���Ŀ¼
rt=read.table("featureGenes.exp.txt", header=T, sep="\t", check.names=F)     #��ȡ�����ļ�
allGenes=unique(as.vector(rt[,1]))
geneType=gsub("(.*?)\\|(.*?)", "\\2", allGenes)
allGenes=gsub("(.*?)\\|.*", "\\1", allGenes)

for(i in levels(factor(geneType))){
	#��������ת��Ϊ����id
	genes=allGenes[geneType==i]
	entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
	entrezIDs=as.character(entrezIDs)
	gene=entrezIDs[entrezIDs!="NA"]        #ȥ������idΪNA�Ļ���
	
	#GO��������
	kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
	GO=as.data.frame(kk)
	GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
	#���渻�����
	write.table(GO, file=paste0(i, ".GO.txt"), sep="\t", quote=F, row.names = F)
	
	#������ʾTerm��Ŀ
	showNum=10
	if(nrow(GO)<30){
		showNum=nrow(GO)
	}
	
	#��״ͼ
	pdf(file=paste0(i, ".barplot.pdf"), width=9, height=7)
	bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
	print(bar)
	dev.off()
			
	#����ͼ
	pdf(file=paste0(i, ".bubble.pdf"), width=9, height=7)
	bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
	print(bub)
	dev.off()
}


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056