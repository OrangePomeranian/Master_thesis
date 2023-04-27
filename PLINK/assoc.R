daneassoc001 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/PLINK/assoc_test_maf_0.01.assoc", header=TRUE)
#install.packages("qqman")
library(qqman)

########################################################
summary(daneassoc001)
assoc <- manhattan(daneassoc001,main = "Manhattan Plot for 0.01 maf", chr="CHR", bp="BP", snp="SNP", p="P", 
                   annotatePval = 0.01, xlab ='')

################################################################################################################

daneassoc005 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/PLINK/assoc_test_maf_0.05.assoc", header=TRUE)
#install.packages("qqman")
library(qqman)

########################################################
summary(daneassoc005)
assoc <- manhattan(daneassoc005,main = "Manhattan Plot for 0.05 maf", chr="CHR", bp="BP", snp="SNP", p="P", 
                   annotatePval = 0.01, xlab ='')
