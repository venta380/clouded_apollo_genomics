library(SNPRelate)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)


#vcf.fn<-"Auto_only_4fold.recode.vcf"
#snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")
snpgdsSummary("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.ids <- data.frame(sample.id)
pop_code <- read.table("pop.group2", header = F)
df3 = merge(sample.ids, pop_code, by.x=c("sample.id"), by.y=c("V1"))


snpset <- snpgdsLDpruning(genofile, ld.threshold=0.4, num.thread=20)
names(snpset)
head(snpset$chr1)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=20)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

sample.id = df3$sample.id
pop = factor(df3$V2)
PCA_stuff=data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
#tab <- data.frame(sample.id = ccm_pca$sample.id, pop = factor(pop_code$V1)[match(ccm_pca$sample.id, sample.id)],colour=factor(pop_code$V2)[match(ccm_pca$sample.id, sample.id)] , EV1 = ccm_pca$eigenvect[,1], EV2 = ccm_pca$eigenvect[,2], stringsAsFactors = FALSE)
tab = merge(PCA_stuff, df3, by.x=c("sample.id"), by.y=c("sample.id"))

write.csv(tab, file = "PCA_final_data.csv")


plot(tab$EV2~tab$EV1, col=as.character(tab$V3), xlab="Eigenvector 2", ylab="Eigenvector 1", cex.lab=1.5, cex = 2, pch=19)
legend("bottomleft", cex = 1.5, pch=19, col=c('#D35400','#F4D03F','#2C3E50','#27AE60'), legend = c('Blekinge','Nordens Ark', 'Roslagen', 'Västernorrland'))


lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.character(tab$V3), labels=lbls, pch=19)


#idb
pop_code <- scan("pop.txt", what=character())
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
ibd <- snpgdsIBDMoM(genofile, sample.id=sample.id, snp.id=snpset.id, maf=0.05, missing.rate=0.05, num.thread=2)

ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)
write.csv(ibd.coeff, file = "Relate.csv")


#dendro
ibs <- snpgdsIBS(genofile, num.thread=20)
pop.idx <- order(pop_code)
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)

plot(x, y, col=race, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")
legend("topleft", legend=levels(race), pch="o", text.col=1:nlevels(race))


set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))
plot(rv2$dendrogram, leaflab="none", ylab="Sequence Dissimilarity")
legend("bottomleft", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)

