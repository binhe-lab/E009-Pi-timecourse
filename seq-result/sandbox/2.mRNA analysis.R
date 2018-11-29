# 1.install packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("affyPLM")
#biocLite("RColorBrewer")
#biocLite("impute")
#biocLite("limma")
#biocLite("pheatmap")
#biocLite("ggplot2")
#biocLite("wateRmelon")
install.packages("scales")
# 2.set the working directory
setwd("C:\\Users\\hp-hp\\GSE\\GSE5609")
library(impute)
library("wateRmelon")

# 3.use the geneid matrix to do analysis
probe_exp<-read.table("GSE5609.txt",header=T,sep="\t",row.names=1)
mat=as.matrix(probe_exp)
mat=impute.knn(mat)
matData=mat$data
matData=matData+0.00001

# 4.normalization
matData=matData[rowMeans(matData)>0.005,]
pdf(file="rawBox.pdf")
boxplot(matData,col = "blue",xaxt = "n",outline = F)
dev.off()
matData = betaqn(matData)
pdf(file="normalBox.pdf")
boxplot(matData,col = "red",xaxt = "n",outline = F)
dev.off()
write.table(matData,file="norm.txt",sep="\t",quote=F)

#### 5. Use the log.pl
library(limma)
eset<-read.table("log_norm.txt",sep='\t',header=T,row.names=1)

# 6.change the sample number
condition=factor(c(rep("healthy",5),rep("DCIS",9)))
design<-model.matrix(~-1+condition)
colnames(design)<-c("DCIS","healthy")
contranst.matrix<-makeContrasts(contrasts="DCIS-healthy",levels=design)
fit<-lmFit(eset,design)
fit1<-contrasts.fit(fit,contranst.matrix)
fit2<-eBayes(fit1)
dif<-topTable(fit2,coef="DCIS-healthy",n=nrow(fit2),adjust="BH")
genesymbol<-rownames(dif)
dif<-cbind(genesymbol,dif)
write.table(dif,file="probeid.Foldchange.txt",sep='\t',quote=F,row.names=F)

# Plot volcano
pdf(file="Volcano.pdf")
yMax=max(-log10(dif$adj.P.Val))
xMax=max(abs(dif$logFC))
plot(dif$logFC,-log10(dif$adj.P.Val), xlab="log2(FC)",ylab="-log10(adj.P.Val)",main="Volcano", xlim=c(-xMax,xMax),ylim=c(0,yMax),yaxs="i",pch=20, cex=0.4,col="grey")
diffSub1=subset(dif, dif$adj.P.Val<0.05 & dif$logFC>1)
diffSub2=subset(dif, dif$adj.P.Val<0.05 & dif$logFC<(-1))
points(diffSub1$logFC,-log10(diffSub1$adj.P.Val), pch=20, col="red",cex=0.4)
points(diffSub2$logFC,-log10(diffSub2$adj.P.Val), pch=20, col="green",cex=0.4)
cut = -log10(0.05)
abline(h=cut,lty=2,lwd=1,col="blue")
abline(v=c(-log2(2),log2(2)),lty=2,lwd=1,col="blue")
dev.off()


# 7.set the foldchange value in the bracket
dif2<-topTable(fit2,coef="DCIS-healthy",n=nrow(fit2),lfc=log2(2),adjust="BH")

# 8.set the adj.P.Val in the bracket
dif2<-dif2[dif2[,"adj.P.Val"]<0.05,]
genesymbol<-rownames(dif2)
dif2<-cbind(genesymbol,dif2)
dif2<-dif2[order(dif2$logFC),]
write.table(dif2,file="diff.probeid.FDR0.05.txt",sep='\t',quote=F,row.names=F)

dif2$genesymbol<-rownames(dif2)
loc<-match(dif2$genesymbol,rownames(eset))
DEG_exp<-eset[loc,]
genesymbol<-rownames(DEG_exp)
DEG_exp<-cbind(genesymbol,DEG_exp)
write.table(DEG_exp,file="probeid.FDR0.05.exprs.txt",sep='\t',quote=F,row.names=F)

#### 9.Out of R studio!!! Use the probeid to geneid.pl and probeid to genesymbol to convert
DEG_exp<-read.table("mRNA.FDR0.05.exprs.txt",sep='\t',header=T,row.names=1)
flag1 <- 53 #narrow down to which pathway
pathway <- read.table("KEGG_human.txt",header=F,sep="\t")
genelist <- as.character(pathway[flag1,2])
b = as.vector(genelist)
c=unlist(strsplit(b,split=";"))
genelist<-c
genelist
pathway_gene_exp<-DEG_exp[intersect(row.names(DEG_exp),genelist),]
row.names(DEG_exp)
intersect(row.names(DEG_exp),genelist)
designNC =  c( rep("healthy",4), rep("DCIS",4))
group_info <- data.frame(row.names=names(DEG_exp),groups=designNC)

top10_down_up <- DEG_exp[c(1:10,(nrow(DEG_exp)-9):nrow(DEG_exp)),]
write.table(top10_down_up,file="top10_down_up_mRNA.txt",quote=F,sep = "\t",row.names=T,col.names = NA)


DEG_exp1<-read.table("top10_down_up_mRNA.txt",sep='\t',header=T,row.names=1)
library(pheatmap)

#pdf(file="pheatmap_mRNA_pathway.pdf",width=20,height=10)
#pheatmap(pathway_gene_exp,color=colorRampPalette(c("green","black","red"))(100),fontsize_row=10
         #,fontsize_col=10,scale="row",border_color=NA,cluster_col = FALSE,annotation_col=group_info)
#dev.off()

pdf(file="pheatmap_top_10.pdf",width=20,height=10)
pheatmap(DEG_exp1,color=colorRampPalette(c("green","black","red"))(100),fontsize_row=10
         ,fontsize_col=10,scale="row",border_color=NA,cluster_col = FALSE,annotation_col=group_info)
dev.off()
