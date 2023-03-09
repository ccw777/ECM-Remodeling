######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("pheatmap")


library(pheatmap)       #引用包
setwd("C:\\Users\\c'w'l\\Desktop\\new\\8.model\\riskplot")      #设置工作目录

#定义风险曲线的函数
bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null, heatmapFile=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
	rt=rt[order(rt$riskScore),]      #按照风险打分对样品排序
		
	#绘制风险曲线
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=riskScoreFile, width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)", ylab="Risk score",
		 col=c(rep("green",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
	dev.off()
		
	#绘制生存状态图
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="green"
	pdf(file=survStatFile, width=7, height=4)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
		
	#绘制风险热图
	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(type=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=heatmapFile, width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation, 
		     cluster_cols = F,
		     cluster_rows = F,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("green",5), "white", rep("red",5)))(50),
		     fontsize_col=3,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}

#调用函数，绘制风险曲线
bioRiskPlot(inputFile="risk.txt",
            riskScoreFile="riskScore.pdf",
            survStatFile="survStat.pdf",
            heatmapFile="heatmap.pdf")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
