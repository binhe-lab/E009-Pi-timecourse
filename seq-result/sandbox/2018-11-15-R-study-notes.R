##how to install package
##install.packages("packagename")


## 5595 objects, 41 variables

sample <- fread("Ex009_experiment_set_up_20171019.csv")
## 40 objects, 6 variables

## limma
source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
help(package=limma) ## userguide in R


setwd("C:\\Users\\Jinye Liang\\Documents\\R\\GSE\\E009-Pi-timecourse")
library(impute)
library(wateRmelon)
library(affy)
library(affyPLM)

probe_exp<-read.table("Ex009.txt",header=T,sep="\t",row.names=1) ##header:whether we have first row as our name for variables, row.names:whether we have the first column for our row names
mat=as.matrix(probe_exp) ##as.matrix() Convert a data.table to a matrix
mat=impute.knn(mat) ## impute.knn a function to impute missing expression data, using nearest neighbor averaging
matData=mat$data
matData=matData+0.00001

matData=matData[rowMeans(matData)>0.005,]
pdf(file="rawBox.pdf")
boxplot(matData,col = "blue",xaxt = "n",outline = F)
dev.off()

matData = betaqn(matData) ##betaqn	quantile normalizes betas, the real process to normalize the data
pdf(file="normalBox.pdf")
boxplot(matData,col = "red",xaxt = "n",outline = F)
dev.off()
write.table(matData,file="norm.txt",sep="\t",quote=F)

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
