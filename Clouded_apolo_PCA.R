#install.packages(c("SNPRelate", "ggplot2", "devtools", "dplyr", "reshape2"))
library(SNPRelate)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)


#sed -i 's/CAJQZP01//' SNPS_out_PCA.recode.vcf
#vcftools --vcf 4D.recode.vcf.recode.vcf --chr /scratch/vt20265/clouded_apollo/Genome/auto_for_PCA --recode  --out Auto_only_4fold




vcf.fn<-"Auto_only_4fold.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")
snpgdsSummary("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.ids <- data.frame(sample.id)
pop_code <- read.table("pop.group", header = F)
df3 = merge(sample.ids, pop_code, by.x=c("sample.id"), by.y=c("V1"))



#ccm_pca<-snpgdsPCA(genofile)
#names(ccm_pca)

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
tab <- read.csv(file="PCA_final_data.csv", header=TRUE, sep=",")
#tab=tab[-c(21,22,47,48,13,17,18), ]
head(tab)

pdf("rplot2.pdf") 
plot(tab$EV2~tab$EV1, col=as.character(tab$V3), xlab="Eigenvector 2", ylab="Eigenvector 1",pch=19, cex.lab=1.5)
legend("bottomleft", pch=c(19), col=c('#D35400','#F4D03F','#2C3E50','#27AE60'), legend = c('Blekinge','NordensArk', 'Uppland', 'Västernorrland'))
#with(tab, text(tab$EV2~tab$EV1, labels = as.character(tab$sample.id)), pos = 4)
dev.off() 


#pdf("rplot2_ids.pdf") 
#plot(tab$EV2~tab$EV1, col=as.character(tab$V3), xlab="Eigenvector 2", ylab="Eigenvector 1",pch=19, cex.lab=1.5)
#with(tab, text(tab$EV2~tab$EV1, labels = as.character(tab$sample.id)), pos = 4)
#dev.off() 
#
#
pdf("rplot_2_Uppland.pdf") 
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.character(tab$V3), labels=lbls, pch=19)
dev.off() 


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

pdf("Heatmap.pdf") 
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))
dev.off() 

loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)

pdf("MDS.pdf") 
plot(x, y, col=race, xlab = "", ylab = "",
    main = "Multidimensional Scaling Analysis (IBS)")
legend("topleft", legend=levels(race), pch="o", text.col=1:nlevels(race))
dev.off() 

set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

pdf("Dendro.pdf") 
plot(rv2$dendrogram, leaflab="none", ylab="Sequence Dissimilarity")
legend("bottomleft", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)

dev.off() 


####
#LEA

library(RColorBrewer)

library(LEA)

jBrewColors <- brewer.pal(n = 10, name = "Paired")
V3<-c("1025","1026","1027","730","731","732","733","737","756","757","871","883","885","962","963","965","ESB1","ESB10","ESB3","ESB4","ESB5","ESB8","ESB9","HH1","HH10","HH2","HH3","HH4","HH5","HH6","HH7","HH8","HH9","HI023","HI033","NJ1","NJ116","NJ203","PL1","PL10","PL2","PL4","PL5","PL6","PL7","PL8","PL9","Plex_BLZ_25_M","Plex_BLZ_4_M","Plex_BLZ_81_F","Plex_BMU_10_M","Plex_BMU_11_F","Plex_BMU_20_F","Plex_CRC_23_M","Plex_CRC_25_F","Plex_CRC_81_F","Plex_FLs_MIA11_F","Plex_FLs_MIA122_M","Plex_FLs_MIA126_F","Plex_FLs_MIA16_M","Plex_FLs_MIA2514_M","Plex_FLs_MIA40_F","Plex_FLs_MIA454_F","Plex_PRI_120_F","Plex_PRI_131_M","Plex_PRI_81_F","T14","T9","mex1527","mex536","mex915","mex919","mex986","stm146","stm163")



West<-data.frame('sample.id'=c("HH1","HH2","HH3","HH4","HH5","HH6","HH7","HH8","HH9","HH10","PL1","PL2","PL4","PL5","PL6","PL7","PL8","PL9","PL10","ESB1","ESB3","ESB4","ESB5","ESB8","ESB9","ESB10"))
East<-data.frame('sample.id'=c("stm163","stm146","T9","T14","NJ203","NJ116","NJ1","HI023","HI033","mex986","mex919","mex915","mex536","mex1527"))
Texas_Migrent<-data.frame('sample.id'=c("730","732","756","757","871","885","883"))
Texas_resident<-data.frame('sample.id'=c("731","733","737","1025","1026","1027","962","963","965"))
BLZ<-data.frame('sample.id'=c("Plex_BLZ_25_M","Plex_BLZ_4_M","Plex_BLZ_81_F"))
BMU<-data.frame('sample.id'=c("Plex_BMU_10_M","Plex_BMU_11_F","Plex_BMU_20_F"))
CRC<-data.frame('sample.id'=c("Plex_CRC_23_M","Plex_CRC_25_F","Plex_CRC_81_F"))
FLs<-data.frame('sample.id'=c("Plex_FLs_MIA11_F","Plex_FLs_MIA122_M","Plex_FLs_MIA126_F","Plex_FLs_MIA16_M","Plex_FLs_MIA2514_M","Plex_FLs_MIA40_F","Plex_FLs_MIA454_F"))
PRI<-data.frame('sample.id'=c("Plex_PRI_120_F","Plex_PRI_131_M","Plex_PRI_81_F"))





obj.snmf =  load.snmfProject("SNPS_out_mig_res.snmfProject")

jBrewColors <- c('#7F7F7F','#4472C4')
#obj.snmf = snmf("temp_calls_SNPs_only_PCA.recode.geno", K = 1:6, ploidy = 2, entropy = T,alpha = 100, project = "new")

pdf(file= './error.pdf' ,onefile=T,paper='A4', )
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
dev.off() 


###

qmatrix1 = Q(obj.snmf, K = 2)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3= V3)

qmat2 <- qmatrix_1[order(qmatrix1[,1]),]

qmat_2 = data.frame(V2 = qmat2[,2], V1 = qmat2[,1])

West_2=merge(qmat2, West, by.x=c("V3"), by.y=c("sample.id"))
East_2=merge(qmat2, East, by.x=c("V3"), by.y=c("sample.id"))
Texas_Migrent_2=merge(qmat2, Texas_Migrent, by.x=c("V3"), by.y=c("sample.id"))
Texas_resident_2=merge(qmat2, Texas_resident, by.x=c("V3"), by.y=c("sample.id"))
BLZ_2=merge(qmat2, BLZ, by.x=c("V3"), by.y=c("sample.id"))
BMU_2=merge(qmat2, BMU, by.x=c("V3"), by.y=c("sample.id"))
CRC_2=merge(qmat2, CRC, by.x=c("V3"), by.y=c("sample.id"))
FLs_2=merge(qmat2, FLs, by.x=c("V3"), by.y=c("sample.id"))
PRI_2=merge(qmat2, PRI, by.x=c("V3"), by.y=c("sample.id"))


West_2=West_2[order(West_2[,2]),]
East_2=East_2[order(East_2[,2]),]
Texas_Migrent_2=Texas_Migrent_2[order(Texas_Migrent_2[,2]),]
Texas_resident_2=Texas_resident_2[order(Texas_resident_2[,2]),]
BLZ_2=BLZ_2[order(BLZ_2[,2]),]
BMU_2=BMU_2[order(BMU_2[,2]),]
CRC_2=CRC_2[order(CRC_2[,2]),]
FLs_2=FLs_2[order(FLs_2[,2]),]
PRI_2=PRI_2[order(PRI_2[,2]),]

West_2=data.frame(V1=West_2$V1, V2=West_2$V2)
East_2=data.frame(V1=East_2$V1, V2=East_2$V2)
Texas_Migrent_2=data.frame(V1=Texas_Migrent_2$V1, V2=Texas_Migrent_2$V2)
Texas_resident_2=data.frame(V1=Texas_resident_2$V1, V2=Texas_resident_2$V2)
BLZ_2=data.frame(V1=BLZ_2$V1, V2=BLZ_2$V2)
BMU_2=data.frame(V1=BMU_2$V1, V2=BMU_2$V2)
CRC_2=data.frame(V1=CRC_2$V1, V2=CRC_2$V2)
FLs_2=data.frame(V1=FLs_2$V1, V2=FLs_2$V2)
PRI_2=data.frame(V1=PRI_2$V1, V2=PRI_2$V2)

pdf(file= './K_2.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,9), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(0.3467,1.0))
plot1=barplot(t(West_2), col = jBrewColors, border = NA, ylab = "K=2", xlab = "West")
par(pin=c(0.1867,1.0))
plot1=barplot(t(East_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "East")
par(pin=c(0.0933,1.0))
plot1=barplot(t(Texas_Migrent_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "Texas_Migrent")
par(pin=c(0.1200,1.0))
plot1=barplot(t(Texas_resident_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "Texas_resident")
par(pin=c(0.0400,1.0))
plot1=barplot(t(BLZ_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "BLZ")
par(pin=c(0.0400,1.0))
plot1=barplot(t(BMU_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "BMU")
par(pin=c(0.0400,1.0))
plot1=barplot(t(CRC_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "CRC")
par(pin=c(0.0933,1.0))
plot1=barplot(t(FLs_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "FLs")
par(pin=c(0.0400,1.0))
plot1=barplot(t(PRI_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "PRI")
dev.off() 


#tab <- read.csv(file="PCA_final_data.csv", header=TRUE, sep=",")
#tab <- data.frame(sample.id=tab$sample.id, colour=tab$V3)
#tab = merge(tab,qmat2, by.y=c("V3"), by.x=c("sample.id"), sort=F)
#tab = tab[order(tab$V1), ]
#pdf(file= './K_2_col.pdf' ,onefile=T,paper='A4', )
#barplot(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), col=tab$colour)
#dev.off() 

###
jBrewColors <- c('#47DDED','#4472C4','#7F7F7F')

qmatrix1 = Q(obj.snmf, K = 3)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3 = qmatrix1[,3], V4= V3)

qmat2 <- qmatrix_1[with(qmatrix_1, order(qmatrix_1[,3]+qmatrix_1[,2])),]
qmat2 <- qmat2[with(qmat2, order(qmat2[,1]+qmat2[,2])),]

qmat_2 = data.frame(V3 = qmat2[,3],V2 = qmat2[,2], V1 = qmat2[,1])


West_2=merge(qmat2, West, by.x=c("V4"), by.y=c("sample.id"))
East_2=merge(qmat2, East, by.x=c("V4"), by.y=c("sample.id"))
Texas_Migrent_2=merge(qmat2, Texas_Migrent, by.x=c("V4"), by.y=c("sample.id"))
Texas_resident_2=merge(qmat2, Texas_resident, by.x=c("V4"), by.y=c("sample.id"))
BLZ_2=merge(qmat2, BLZ, by.x=c("V4"), by.y=c("sample.id"))
BMU_2=merge(qmat2, BMU, by.x=c("V4"), by.y=c("sample.id"))
CRC_2=merge(qmat2, CRC, by.x=c("V4"), by.y=c("sample.id"))
FLs_2=merge(qmat2, FLs, by.x=c("V4"), by.y=c("sample.id"))
PRI_2=merge(qmat2, PRI, by.x=c("V4"), by.y=c("sample.id"))


West_2=West_2[with(West_2, order(West_2[,4]+West_2[,3])),]
East_2=East_2[with(East_2, order(East_2[,4]+East_2[,3])),]
Texas_Migrent_2=Texas_Migrent_2[with(Texas_Migrent_2, order(Texas_Migrent_2[,4]+Texas_Migrent_2[,3])),]
Texas_resident_2=Texas_resident_2[with(Texas_resident_2, order(Texas_resident_2[,4]+Texas_resident_2[,3])),]
BLZ_2=BLZ_2[with(BLZ_2, order(BLZ_2[,4]+BLZ_2[,3])),]
BMU_2=BMU_2[with(BMU_2, order(BMU_2[,4]+BMU_2[,3])),]
CRC_2=CRC_2[with(CRC_2, order(CRC_2[,4]+CRC_2[,3])),]
FLs_2=FLs_2[with(FLs_2, order(FLs_2[,4]+FLs_2[,3])),]
PRI_2=PRI_2[with(PRI_2, order(PRI_2[,4]+PRI_2[,3])),]

West_2=West_2[with(West_2, order(West_2[,2]+West_2[,3])),]
East_2=East_2[with(East_2, order(East_2[,2]+East_2[,3])),]
Texas_Migrent_2=Texas_Migrent_2[with(Texas_Migrent_2, order(Texas_Migrent_2[,2]+Texas_Migrent_2[,3])),]
Texas_resident_2=Texas_resident_2[with(Texas_resident_2, order(Texas_resident_2[,2]+Texas_resident_2[,3])),]
BLZ_2=BLZ_2[with(BLZ_2, order(BLZ_2[,2]+BLZ_2[,3])),]
BMU_2=BMU_2[with(BMU_2, order(BMU_2[,2]+BMU_2[,3])),]
CRC_2=CRC_2[with(CRC_2, order(CRC_2[,2]+CRC_2[,3])),]
FLs_2=FLs_2[with(FLs_2, order(FLs_2[,2]+FLs_2[,3])),]
PRI_2=PRI_2[with(PRI_2, order(PRI_2[,2]+PRI_2[,3])),]

West_2=data.frame(V1=West_2$V1, V2=West_2$V2, V3=West_2$V3)
East_2=data.frame(V1=East_2$V1, V2=East_2$V2, V3=East_2$V3)
Texas_Migrent_2=data.frame(V1=Texas_Migrent_2$V1, V2=Texas_Migrent_2$V2, V3=Texas_Migrent_2$V3)
Texas_resident_2=data.frame(V1=Texas_resident_2$V1, V2=Texas_resident_2$V2, V3=Texas_resident_2$V3)
BLZ_2=data.frame(V1=BLZ_2$V1, V2=BLZ_2$V2, V3=BLZ_2$V3)
BMU_2=data.frame(V1=BMU_2$V1, V2=BMU_2$V2, V3=BMU_2$V3)
CRC_2=data.frame(V1=CRC_2$V1, V2=CRC_2$V2, V3=CRC_2$V3)
FLs_2=data.frame(V1=FLs_2$V1, V2=FLs_2$V2, V3=FLs_2$V3)
PRI_2=data.frame(V1=PRI_2$V1, V2=PRI_2$V2, V3=PRI_2$V3)



pdf(file= './K_3.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,9), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(0.3467,1.0))
plot1=barplot(t(West_2), col = jBrewColors, border = NA, ylab = "K=3", xlab = "West")
par(pin=c(0.1867,1.0))
plot1=barplot(t(East_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "East")
par(pin=c(0.0933,1.0))
plot1=barplot(t(Texas_Migrent_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "Texas_Migrent")
par(pin=c(0.1200,1.0))
plot1=barplot(t(Texas_resident_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "Texas_resident")
par(pin=c(0.0400,1.0))
plot1=barplot(t(BLZ_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "BLZ")
par(pin=c(0.0400,1.0))
plot1=barplot(t(BMU_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "BMU")
par(pin=c(0.0400,1.0))
plot1=barplot(t(CRC_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "CRC")
par(pin=c(0.0933,1.0))
plot1=barplot(t(FLs_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "FLs")
par(pin=c(0.0400,1.0))
plot1=barplot(t(PRI_2), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "PRI")
dev.off() 

tab <- read.csv(file="PCA_final_data.csv", header=TRUE, sep=",")
tab <- data.frame(sample.id=tab$sample.id, colour=tab$V3)
tab = merge(tab,qmat2, by.y=c("V4"), by.x=c("sample.id"), sort=F)
tab = tab[with(tab, order(tab$V3+tab$V2)),]
tab = tab[with(tab, order(tab$V1+tab$V2)),]

pdf(file= './K_3_col.pdf' ,onefile=T,paper='A4', )
barplot(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), col=tab$colour)
dev.off() 

###

qmatrix1 = Q(obj.snmf, K = 4)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3 = qmatrix1[,3], V4 = qmatrix1[,4],V5= V3)

qmat2 <- qmatrix_1[with(qmatrix_1, order(qmatrix_1[,3]+qmatrix_1[,2])),]
qmat2 <- qmat2[with(qmat2, order(qmat2[,1]+qmat2[,2])),]

qmat_2 = data.frame(V4 = qmat2[,4], V3 = qmat2[,3],V2 = qmat2[,2], V1 = qmat2[,1])

pdf(file= './K_4.pdf' ,onefile=T,paper='A4', )
plot1=barplot(t(qmat_2), col = jBrewColors, border = NA, ylab = "K=2")
dev.off() 
