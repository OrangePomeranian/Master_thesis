library(dplyr)
library(SNPassoc)
library(qqman)

# Editing files -----------------------------------------------------
chr = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.01.ped", header = F) #zmiana
map = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.01.map", header = F) #zmiana
casecontrol = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/casecontrol.csv",sep = ',', header = F)
#casecontrol = data.frame(casecontrol$X1)
dim(chr)
chr_temp = chr[,7:ncol(chr)]
for (i in 1:nrow(chr_temp))
{
  for (j in 1:ncol(chr_temp))
  {
    if (chr_temp[i,j] == 'TRUE')
      chr_temp[i,j] = 'T'
  }
}
temp = matrix(NA,nrow(chr_temp),ncol(chr_temp)/2)
for (i in 1:ncol(temp))
{
  temp[,i] = paste(chr_temp[,2*i-1], sep = "", chr_temp[,2*i])
}
temp[1:5,1:10]
for (i in 1:nrow(temp))
{
  for (j in 1:ncol(temp))
  {
    if (temp[i,j] == '00')
      temp[i,j] = NA
  }
}
temp.chr = cbind(casecontrol,temp)
colnames(temp.chr) = c("casecontrol",map[,2])
# write.csv(temp.chr,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.01.csv', row.names=TRUE)
sum(temp.chr[,2:ncol(temp.chr)]=='00')
sum(temp=='TRUE')

### ________________________________________________________________________________
chr = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.ped", header = F) #zmiana
map = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.map", header = F) #zmiana
casecontrol = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/casecontrol.csv",sep = ',', header = F)
#casecontrol = data.frame(casecontrol$X1)
dim(chr)
chr_temp = chr[,7:ncol(chr)]
for (i in 1:nrow(chr_temp))
{
  for (j in 1:ncol(chr_temp))
  {
    if (chr_temp[i,j] == 'TRUE')
      chr_temp[i,j] = 'T'
  }
}
temp = matrix(NA,nrow(chr_temp),ncol(chr_temp)/2)
for (i in 1:ncol(temp))
{
  temp[,i] = paste(chr_temp[,2*i-1], sep = "", chr_temp[,2*i])
}
temp[1:5,1:10]
for (i in 1:nrow(temp))
{
  for (j in 1:ncol(temp))
  {
    if (temp[i,j] == '00')
      temp[i,j] = NA
  }
}
temp.chr = cbind(casecontrol,temp)
colnames(temp.chr) = c("casecontrol",map[,2])
# write.csv(temp.chr,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.05.csv', row.names=TRUE)
sum(temp.chr[,2:ncol(temp.chr)]=='00')
sum(temp=='TRUE')

# Significant -------------------------------------------------------------
library(dplyr)
library(SNPassoc)

output_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.05.csv", header = TRUE, sep =',')
output_maf_0.01 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.01.csv", header = TRUE, sep =',')
inf_maf_0.01 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.01.map", header = FALSE)
inf_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.map", header = FALSE)
output_maf_0.05 = output_maf_0.05[,2:3745]
output_maf_0.01 = output_maf_0.01[,2:4052]
inf_maf_0.01 = inf_maf_0.01[,2:4]
inf_maf_0.05 = inf_maf_0.05[,2:4]

#### dla 0.05
myDat.HapMap_1<-setupSNP(output_maf_0.05, colSNPs=2:2000, info=inf_maf_0.05, sep='')
myDat.HapMap_2<-setupSNP(output_maf_0.05, colSNPs=2000:3744, info=inf_maf_0.05, sep="")
#test<-setupSNP(output_maf_0.05, colSNPs=2:3745, info=inf_maf_0.05, sep="")

resHapMap_1<-WGassociation(casecontrol, data=myDat.HapMap_1)
resHapMap_2<-WGassociation(casecontrol, data=myDat.HapMap_2)
#resHapMap_2<-WGassociation(casecontrol, data=test)
resHapMap_0.05 = rbind(resHapMap_1, resHapMap_2)

# write.csv(resHapMap_0.05,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_0.05.csv', row.names=TRUE)


#### dla 0.01
myDat.HapMap_1<-setupSNP(output_maf_0.01, colSNPs=2:2000, info=inf_maf_0.01, sep="")
myDat.HapMap_2<-setupSNP(output_maf_0.01, colSNPs=2000:3999, info=inf_maf_0.01, sep="")
myDat.HapMap_3<-setupSNP(output_maf_0.01, colSNPs=3999:4050, info=inf_maf_0.01, sep="")
#test<-setupSNP(output_maf_0.01, colSNPs=2:4051, info=inf_maf_0.01, sep="")

resHapMap_1<-WGassociation(casecontrol, data=myDat.HapMap_1)
resHapMap_2<-WGassociation(casecontrol, data=myDat.HapMap_2)
resHapMap_3<-WGassociation(casecontrol, data=myDat.HapMap_3)

resHapMap_0.01 = rbind(resHapMap_1, resHapMap_2, resHapMap_3)

# write.csv(resHapMap_0.01,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_0.01.csv', row.names=TRUE)


# Manhattan_plots ---------------------------------------------------------
sig_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_out_0.05.csv", header = TRUE, sep =',')
sig_0.01 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_out_0.05.csv", header = TRUE, sep =',')

p_value = 0.00001
for (x in names(sig_0.05)[5:9]) {
  dane_wzor = sig_0.05[1:5000,2:4]
  dane_wzor[x] = sig_0.05[1:5000,][x]
  dane_wzor = na.omit(dane_wzor)
  manhattan(dane_wzor, chr="CHR", bp="BP", snp="SNP", p=x, xaxt = 'n', annotateTop = FALSE,
            xlab ='',main = paste("Manhattan Plot for 0.05 for", x, "model",sep = " "),
            annotatePval = p_value)  
}

for (x in names(sig_0.01)[5:9]) {
  dane_wzor = sig_0.01[1:5000,2:4]
  dane_wzor[x] = sig_0.01[1:5000,][x]
  dane_wzor = na.omit(dane_wzor)
  manhattan(dane_wzor, chr="CHR", bp="BP", snp="SNP", p=x, xaxt = 'n', annotateTop = FALSE,
            xlab ='',main = paste("Manhattan Plot for 0.01 for", x, "model",sep = " "), 
            annotatePval = p_value)
}
# Assosiaction ------------------------------------------------------------
sig_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_out_0.05.csv", header = TRUE, sep =',')
sig_0.01 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_out_0.05.csv", header = TRUE, sep =',')

p_value = 0.00001 
sig_dane_0.05 <- data.frame()

for (x in names(sig_0.05)[5:9]) {
  dane_wzor = sig_0.05[1:50000,2:4]
  dane_wzor[x] = sig_0.05[1:5000,][x]
  dane_wzor = na.omit(dane_wzor)
  dane_wzor  <- dane_wzor[dane_wzor[,x] < p_value,]
  dane_wzor$model = x
  sig_dane_0.05 = bind_rows(sig_dane_0.05, dane_wzor)
}


output_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.05.csv", header = TRUE, sep =',')
output_maf_0.01 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.01.csv", header = TRUE, sep =',')
inf_maf_0.01 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.01.map", header = FALSE)
inf_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.map", header = FALSE)
getSignificantSNPs(resHapMap_0.05,15, model = 'log')
association(casecontrol~snp(TIGRP2P2,sep=""), data=output_maf_0.05)

