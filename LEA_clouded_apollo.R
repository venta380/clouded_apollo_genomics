library(SNPRelate)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(LEA)



obj.snmf = snmf("Auto_only_4fold.geno", K = 1:6, ploidy = 2, entropy = T,alpha = 100, project = "new", CPU = 20)

obj.snmf =  load.snmfProject("Auto_only_4fold.snmfProject")



pdf(file= './CV.pdf' ,onefile=T,paper='A4', )
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
dev.off() 


V3<-c("S121","S122","S123","S124","S125","S126","S127","S128","S129","S130","S131","S132","S133","S134","S135","S136","S137","S138","S139","S140","S141","S142"   ,"S143"   ,"S144"   ,"S145"   ,"S146"   ,"S147"   ,"S148"   ,"S149"   ,"S151"   ,"S152"   ,"S153"   ,"S154"   ,"S155"   ,"S156"   ,"S157"   ,"S158"   ,"S159")


Uppland<-data.frame('sample.id'=c("S121","S122","S123","S124","S125","S126","S127","S128","S129","S130"))
Blekinge<-data.frame('sample.id'=c("S131","S132","S133","S134","S135","S136","S137","S138","S139","S140"))
NordensArk<-data.frame('sample.id'=c("S141","S142","S143","S144","S145","S146","S147","S148","S149"))
Västernorrland<-data.frame('sample.id'=c("S151","S152","S153","S154","S155","S156","S157","S158","S159"))





#plot K=2
jBrewColors <- c('#2C3E50','#27AE60')

qmatrix1 = Q(obj.snmf, K = 2)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3= V3)

qmat2 <- qmatrix_1[order(qmatrix1[,1]),]

qmat_2 = data.frame(V2 = qmat2[,2], V1 = qmat2[,1])

Uppland<-data.frame('sample.id'=c("S121","S122","S123","S124","S125","S126","S127","S128","S129","S130"))
Blekinge<-data.frame('sample.id'=c("S131","S132","S133","S134","S135","S136","S137","S138","S139","S140"))
NordensArk<-data.frame('sample.id'=c("S141","S142","S143","S144","S145","S146","S147","S148","S149"))
Västernorrland<-data.frame('sample.id'=c("S151","S152","S153","S154","S155","S156","S157","S158","S159"))




Uppland=merge(qmat2, Uppland, by.x=c("V3"), by.y=c("sample.id"))
Blekinge=merge(qmat2, Blekinge, by.x=c("V3"), by.y=c("sample.id"))
NordensArk=merge(qmat2, NordensArk, by.x=c("V3"), by.y=c("sample.id"))
Västernorrland=merge(qmat2, Västernorrland, by.x=c("V3"), by.y=c("sample.id"))

Uppland=data.frame(V1=Uppland$V1, V2=Uppland$V2)
Blekinge=data.frame(V1=Blekinge$V1, V2=Blekinge$V2)
NordensArk=data.frame(V1=NordensArk$V1, V2=NordensArk$V2)
Västernorrland=data.frame(V1=Västernorrland$V1, V2=Västernorrland$V2)


pdf(file= './K_2.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,4), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(nrow(Uppland)/nrow(qmat2),1.0))
plot1=barplot(t(Uppland), col = jBrewColors, yaxt='n', border = NA, ylab = "K=2", xlab = "Uppland")
par(pin=c(nrow(Blekinge)/nrow(qmat2),1.0))
plot1=barplot(t(Blekinge), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Blekinge")
par(pin=c(nrow(NordensArk)/nrow(qmat2),1.0))
plot1=barplot(t(NordensArk), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "NordensArk")
par(pin=c(nrow(Västernorrland)/nrow(qmat2),1.0))
plot1=barplot(t(Västernorrland), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Västernorrland")

dev.off() 







#plot K=3
jBrewColors <- c('#2C3E50','#D35400', '#27AE60' )

qmatrix1 = Q(obj.snmf, K = 3)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3 = qmatrix1[,3], V4= V3)

qmat2 <- qmatrix_1[order(qmatrix1[,1]),]

qmat_2 = data.frame(V3 = qmat2[,3], V2 = qmat2[,2], V1 = qmat2[,1])


Uppland<-data.frame('sample.id'=c("S121","S122","S123","S124","S125","S126","S127","S128","S129","S130"))
Blekinge<-data.frame('sample.id'=c("S131","S132","S133","S134","S135","S136","S137","S138","S139","S140"))
NordensArk<-data.frame('sample.id'=c("S141","S142","S143","S144","S145","S146","S147","S148","S149"))
Västernorrland<-data.frame('sample.id'=c("S151","S152","S153","S154","S155","S156","S157","S158","S159"))

Uppland=merge(qmat2, Uppland, by.x=c("V4"), by.y=c("sample.id"))
Blekinge=merge(qmat2, Blekinge, by.x=c("V4"), by.y=c("sample.id"))
NordensArk=merge(qmat2, NordensArk, by.x=c("V4"), by.y=c("sample.id"))
Västernorrland=merge(qmat2, Västernorrland, by.x=c("V4"), by.y=c("sample.id"))

Uppland=data.frame(V1=Uppland$V1, V2=Uppland$V2, V3=Uppland$V3)
Blekinge=data.frame(V1=Blekinge$V1, V2=Blekinge$V2, V3=Blekinge$V3)
NordensArk=data.frame(V1=NordensArk$V1, V2=NordensArk$V2, V3=NordensArk$V3)
Västernorrland=data.frame(V1=Västernorrland$V1, V2=Västernorrland$V2, V3=Västernorrland$V3)


pdf(file= './K_3.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,4), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(nrow(Uppland)/nrow(qmat2),1.0))
plot1=barplot(t(Uppland), col = jBrewColors, yaxt='n', border = NA, ylab = "K=3", xlab = "Uppland")
par(pin=c(nrow(Blekinge)/nrow(qmat2),1.0))
plot1=barplot(t(Blekinge), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Blekinge")
par(pin=c(nrow(NordensArk)/nrow(qmat2),1.0))
plot1=barplot(t(NordensArk), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "NordensArk")
par(pin=c(nrow(Västernorrland)/nrow(qmat2),1.0))
plot1=barplot(t(Västernorrland), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Västernorrland")

dev.off() 





#plot K=4
jBrewColors <- c('#2C3E50','#27AE60', '#D35400', '#F4D03F')

qmatrix1 = Q(obj.snmf, K = 4)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3 = qmatrix1[,3], V4 = qmatrix1[,4],V5= V3)

qmat2 <- qmatrix_1[order(qmatrix1[,1]),]

qmat_2 = data.frame(V4 = qmat2[,4], V3 = qmat2[,3], V2 = qmat2[,2], V1 = qmat2[,1])



Uppland<-data.frame('sample.id'=c("S121","S122","S123","S124","S125","S126","S127","S128","S129","S130"))
Blekinge<-data.frame('sample.id'=c("S131","S132","S133","S134","S135","S136","S137","S138","S139","S140"))
NordensArk<-data.frame('sample.id'=c("S141","S142","S143","S144","S145","S146","S147","S148","S149"))
Västernorrland<-data.frame('sample.id'=c("S151","S152","S153","S154","S155","S156","S157","S158","S159"))

Uppland=merge(qmat2, Uppland, by.x=c("V5"), by.y=c("sample.id"))
Blekinge=merge(qmat2, Blekinge, by.x=c("V5"), by.y=c("sample.id"))
NordensArk=merge(qmat2, NordensArk, by.x=c("V5"), by.y=c("sample.id"))
Västernorrland=merge(qmat2, Västernorrland, by.x=c("V5"), by.y=c("sample.id"))

Uppland=data.frame(V1=Uppland$V1, V2=Uppland$V2, V3=Uppland$V3, V4=Uppland$V4)
Blekinge=data.frame(V1=Blekinge$V1, V2=Blekinge$V2, V3=Blekinge$V3, V4=Blekinge$V4)
NordensArk=data.frame(V1=NordensArk$V1, V2=NordensArk$V2, V3=NordensArk$V3, V4=NordensArk$V4)
Västernorrland=data.frame(V1=Västernorrland$V1, V2=Västernorrland$V2, V3=Västernorrland$V3, V4=Västernorrland$V4)











pdf(file= './K_4.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,4), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(nrow(Uppland)/nrow(qmat2),1.0))
plot1=barplot(t(Uppland), col = jBrewColors, yaxt='n', border = NA, ylab = "K=4", xlab = "Uppland")
par(pin=c(nrow(Blekinge)/nrow(qmat2),1.0))
plot1=barplot(t(Blekinge), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Blekinge")
par(pin=c(nrow(NordensArk)/nrow(qmat2),1.0))
plot1=barplot(t(NordensArk), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "NordensArk")
par(pin=c(nrow(Västernorrland)/nrow(qmat2),1.0))
plot1=barplot(t(Västernorrland), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Västernorrland")

dev.off() 




#plot K=5
jBrewColors <- c('#2C3E50','#27AE60', '#D35400', '#F4D03F', '#1E90FF')

qmatrix1 = Q(obj.snmf, K = 5)
qmatrix_1 = data.frame(V1 = qmatrix1[,1], V2 = qmatrix1[,2], V3 = qmatrix1[,3], V4 = qmatrix1[,4],V5 = qmatrix1[,5] , V6= V3)

qmat2 <- qmatrix_1[order(qmatrix1[,1]),]

qmat_2 = data.frame(V5 = qmat2[,5], V4 = qmat2[,4], V3 = qmat2[,3], V2 = qmat2[,2], V1 = qmat2[,1])



Uppland<-data.frame('sample.id'=c("S121","S122","S123","S124","S125","S126","S127","S128","S129","S130"))
Blekinge<-data.frame('sample.id'=c("S131","S132","S133","S134","S135","S136","S137","S138","S139","S140"))
NordensArk<-data.frame('sample.id'=c("S141","S142","S143","S144","S145","S146","S147","S148","S149"))
Västernorrland<-data.frame('sample.id'=c("S151","S152","S153","S154","S155","S156","S157","S158","S159"))

Uppland=merge(qmat2, Uppland, by.x=c("V6"), by.y=c("sample.id"))
Blekinge=merge(qmat2, Blekinge, by.x=c("V6"), by.y=c("sample.id"))
NordensArk=merge(qmat2, NordensArk, by.x=c("V6"), by.y=c("sample.id"))
Västernorrland=merge(qmat2, Västernorrland, by.x=c("V6"), by.y=c("sample.id"))

Uppland=data.frame(V1=Uppland$V1, V2=Uppland$V2, V3=Uppland$V3, V4=Uppland$V4, V5=Uppland$V5)
Blekinge=data.frame(V1=Blekinge$V1, V2=Blekinge$V2, V3=Blekinge$V3, V4=Blekinge$V4, V5=Blekinge$V5)
NordensArk=data.frame(V1=NordensArk$V1, V2=NordensArk$V2, V3=NordensArk$V3, V4=NordensArk$V4, V5=NordensArk$V5)
Västernorrland=data.frame(V1=Västernorrland$V1, V2=Västernorrland$V2, V3=Västernorrland$V3, V4=Västernorrland$V4, V5=Västernorrland$V5)


pdf(file= './K_5.pdf' ,onefile=T,paper='A4', )
par(mfrow=c(1,4), mai = c(1, 0.1, 0.1, 0.1))
#layout(hights=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), widths=c(0.3467, 0.1867, 0.0933, 0.1200, 0.0400, 0.0400, 0.0400, 0.0933, 0.0400))
par(pin=c(nrow(Uppland)/nrow(qmat2),1.0))
plot1=barplot(t(Uppland), col = jBrewColors, yaxt='n', border = NA, ylab = "K=5", xlab = "Uppland")
par(pin=c(nrow(Blekinge)/nrow(qmat2),1.0))
plot1=barplot(t(Blekinge), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Blekinge")
par(pin=c(nrow(NordensArk)/nrow(qmat2),1.0))
plot1=barplot(t(NordensArk), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "NordensArk")
par(pin=c(nrow(Västernorrland)/nrow(qmat2),1.0))
plot1=barplot(t(Västernorrland), col = jBrewColors, yaxt='n', border = NA, ylab = "", xlab = "Västernorrland")

dev.off() 

