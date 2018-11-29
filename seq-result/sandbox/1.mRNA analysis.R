# 0.install packages
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("RColorBrewer")
biocLite("impute")
biocLite("limma")
biocLite("pheatmap")
biocLite("ggplot2")
biocLite("wateRmelon")
limmaUsersGuide(view=TRUE)

# 1.set the working directory
setwd("C:\\Users\\hp-hp\\GSE\\GSE5609")
library(impute)
library(wateRmelon)
library(affy)
library(affyPLM)

# 2.use the geneid matrix to do analysis
probe_exp<-read.table("GSE5609.txt",header=T,sep="\t",row.names=1)
mat=as.matrix(probe_exp)
mat=impute.knn(mat)
matData=mat$data
matData=matData+0.00001

# 3.rawdata exhibition
matData=matData[rowMeans(matData)>0.005,]
pdf(file="rawBox.pdf")
boxplot(matData,col = "blue",xaxt = "n",outline = F)
dev.off()

# 4.normalization
matData = betaqn(matData)
pdf(file="normalBox.pdf")
boxplot(matData,col = "red",xaxt = "n",outline = F)
dev.off()
write.table(matData,file="norm.txt",sep="\t",quote=F)

# 5.don't use the log.pl
library(limma)
eset<-read.table("norm.txt",sep='\t',header=T,row.names=1)

# 6.change the sample number
condition=factor(c(rep("Nor",2),rep("Aggre",2)))

# 7.calculate foldchange
design<-model.matrix(~-1+condition)
colnames(design)<-c("Nor","Aggre")
contranst.matrix<-makeContrasts(contrasts="Aggre-Nor",levels=design)
fit<-lmFit(eset,design)
fit1<-contrasts.fit(fit,contranst.matrix)
fit2<-eBayes(fit1)
dif<-topTable(fit2,coef="Aggre-Nor",n=nrow(fit2),adjust="BH")
genesymbol<-rownames(dif)
dif<-cbind(genesymbol,dif)
write.table(dif,file="probeid.Foldchange.txt",sep='\t',quote=F,row.names=F)

# 8.plot volcano
pdf(file="Volcano.pdf")
yMax=max(-log10(dif$adj.P.Val))
xMax=max(abs(dif$logFC))
plot(dif$logFC,-log10(dif$adj.P.Val), xlab="log2(FC)",ylab="-log10(adj.P.Val)",main="Volcano", xlim=c(-xMax,xMax),ylim=c(0,yMax),yaxs="i",pch=20, cex=0.4,col="grey")
diffSub1=subset(dif, dif$adj.P.Val<0.05 & dif$logFC>2)
diffSub2=subset(dif, dif$adj.P.Val<0.05 & dif$logFC<(-2))
points(diffSub1$logFC,-log10(diffSub1$adj.P.Val), pch=20, col="red",cex=0.4)
points(diffSub2$logFC,-log10(diffSub2$adj.P.Val), pch=20, col="green",cex=0.4)
cut = -log10(0.05)
abline(h=cut,lty=2,lwd=1,col="blue")
abline(v=c(-log2(4),log2(4)),lty=2,lwd=1,col="blue")
dev.off()


# 9.sreening diff probeid
dif2<-topTable(fit2,coef="Aggre-Nor",n=nrow(fit2),lfc=log2(2),adjust="BH")
dif2<-dif2[dif2[,"adj.P.Val"]<0.05,]
genesymbol<-rownames(dif2)
dif2<-cbind(genesymbol,dif2)
write.table(dif2,file="diff.probeid.FDR0.05.txt",sep='\t',quote=F,row.names=F)

# 10.exp for diff probeid
dif2$genesymbol<-rownames(dif2)
loc<-match(dif2$genesymbol,rownames(eset))
DEG_exp<-eset[loc,]
genesymbol<-rownames(DEG_exp)
DEG_exp<-cbind(genesymbol,DEG_exp)
write.table(DEG_exp,file="probeid.FDR0.05.exprs.txt",sep='\t',quote=F,row.names=F)

# 11.Out of R studio!!! Use the probeid to geneid.pl and probeid to genesymbol to convert

# 12. plot pheatmap
library(pheatmap)
DEG_exp<-read.table("genesymbol.FDR0.05.exprs.txt",
                    sep='\t',header=T,row.names=1)
pdf(file="pheatmap.pdf",width=20,height=10)
pheatmap(DEG_exp,
         color=colorRampPalette(c("green","black","red"))(100),
         fontsize_row=30,
         fontsize_col=30,scale="row",
         border_color=NA,cluster_col = FALSE)
dev.off()
