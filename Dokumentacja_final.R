# W pliku znajdują się kody wykorzystane w pracy magisterskiej pt. 
# "Analiza markerów genetycznych nowotworu listwy mlecznej u psa domowego (Canis lupus familiaris) 
# w obrębie chromosomu 15" autorstwa Darii Plewy, kierunek Bioinformatyka, studia II stopnia, 
# pod kierunkiem dr inż. Anny Muchy z Katedry Genetyki na Wydziale Biologii i Hodowli 
# Zwierząt Uniwersytetu Przyrodniczego we Wrocławiu.


# Niniejszy plik zawiera zestawienie komend użytych do przetwarzania i analizy danych genetycznych
# dotyczących markerów SNP, mających związek z nowotworem listwy mlecznej u psów. Analiza opiera się na
# narzędziu PLINK i zawiera kroki takie jak filtrowanie SNP-ów, testy równowagi Hardy'ego-Weinberga (HWE),
# testy asocjacji genotypów z fenotypem, a także tworzenie bloków haplotypowych.

# Do przetwarzania plików wykorzystano również narzędzie AWK, które służyło do wyciągania i filtrowania 
# odpowiednich danych, takich jak identyfikatory SNP-ów. Było to szczególnie przydatne w etapach przygotowywania
# danych wejściowych do dalszych analiz.

# Wizualizacja wyników analizy, w tym wykresy Manhattan związane z asocjacją SNP-ów z fenotypem, 
# zostały wykonana za pomocą języka R, który pozwala na graficzne przedstawienie istotnych wyników.

# Do identyfikacji haplotypów wykorzystano R w połączeniu z algorytmem EM, 
# co umożliwiło precyzyjną identyfikację prawdopodobnych haplotypów oraz ocenę ich związku z fenotypem 
# nowotworowym u psów.



# Plik został podzielony na 7 części, z których każda zawiera szczegółowy opis zastosowanego oprogramowania 
# oraz komend, wraz z komentarzami, co pozwala na lepsze zrozumienie poszczególnych etapów analizy.



# Ważne: Z pliku najlepiej korzystać z RStudio, które obsługuje opcję zwijania kodu, 
# ułatwiając tym samym poruszanie się między różnymi etapami analizy.


# Terminal/Wiersz Poleceń - Plink i AWK - Filtracja danych ------------------------------------------------

# Filtruje SNP-y z brakującymi danymi powyżej 5%
./plink --file out --1 --recode --geno 0.05 --allow-no-sex --out geno0.05

# Usuwa SNP-y o częstości allelu mniejszej niz 5% 
./plink --file geno0.05 --recode --maf 0.05 --out maf0.05                

# Filtruje SNP-y na podstawie testu równowagi Hardy'ego-Weinberga (HWE) z p < 1e-6 
./plink --file maf0.05  --recode --hwe 1e-6 --out hwe_1_part             

# Filtruje SNP-y na HWE z bardziej restrykcyjnym progiem p < 1e-10
./plink --file hwe_1_part --recode --hwe 1e-10 --hwe-all --out hwe_2_part

# Przeprowadza test asocjacji genotypów z fenotypem na przefiltrowanych danych (MAF 0.05) 
./plink --file hwe_2_part  --recode --assoc --allow-no-sex --out assoc_test_maf_0.05

# Wydobywa identyfikatory SNP-ów z drugiej kolumny pliku .map i zapisuje je
awk '{print $2}' hwe_2_part.map > selected_variants.txt


# R - Manhattan Plot ------------------------------------------------------

# Ładowanie bibliotek
library(qqman)

# Wczytanie danych
daneassoc005 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/PLINK/assoc_test_maf_0.05.assoc", header=TRUE)

# Przedstawienie wyników
summary(daneassoc005)
assoc <- manhattan(daneassoc005,main = "Manhattan Plot for 0.05 maf", chr="CHR", bp="BP", snp="SNP", p="P", 
                   annotatePval = 0.01, xlab ='')


# R - SNPassoc - Przygotowanie Plików ----------------------------------------------------------------

# Załadowanie bibliotek
library(dplyr)
library(SNPassoc)
library(qqman)

# Wczytanie plików
chr = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.ped", header = F) 
map = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.map", header = F) 
casecontrol = read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/casecontrol.csv",sep = ',', header = F)
#casecontrol = data.frame(casecontrol$X1)

# Połączenie plików i zapisanie pliku końcowego zawierającego kolumny z numerem osobnika, wartością case/control 
# oraz identyfikatorami SNP 
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
write.csv(temp.chr,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.05.csv', row.names=TRUE)

# Sprawdzenie pliku
sum(temp.chr[,2:ncol(temp.chr)]=='00')
sum(temp=='TRUE')

# R - SNPassoc - Associacja ----------------------------------------------------------------

# Ładowanie bibliotek
library(dplyr)
library(SNPassoc)

# Wczytanie danych
output_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.05.csv", header = TRUE, sep =',')
inf_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.map", header = FALSE)

# Selekcja kolumn
output_maf_0.05 = output_maf_0.05[,2:3745]
inf_maf_0.05 = inf_maf_0.05[,2:4]

# Podział danych na 2 części
myDat.HapMap_1<-setupSNP(output_maf_0.05, colSNPs=2:2000, info=inf_maf_0.05, sep='')
myDat.HapMap_2<-setupSNP(output_maf_0.05, colSNPs=2000:3744, info=inf_maf_0.05, sep="")

# Analiza asocjacyjna
resHapMap_1<-WGassociation(casecontrol, data=myDat.HapMap_1)
resHapMap_2<-WGassociation(casecontrol, data=myDat.HapMap_2)

# Łaczenie wyników
resHapMap_0.05 = rbind(resHapMap_1, resHapMap_2)

# Zapisanie wyników
write.csv(resHapMap_0.05,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_0.05.csv', row.names=TRUE)

# Tworzenie tabeli z wynikami
sig_out = as.data.frame(table(rownames(resHapMap_0.05), 
                              resHapMap_0.05$codominant))

# Dodawanie kolumn do ramki danych
sig_out$SNP = rownames(resHapMap_0.05)
sig_out$CHR = 15
sig_out$BP = inf_maf_0.05$V4

# Dodawanie wyników dla różnych modeli asocjacyjnych
sig_out$codominant = resHapMap_0.05$codominant
sig_out$dominant = resHapMap_0.05$dominant
sig_out$recessive = resHapMap_0.05$recessive
sig_out$overdominant = resHapMap_0.05$overdominant
sig_out$`log-additive` = resHapMap_0.05$`log-additive`

# Przeksztalcenie struktury danych
sig_out = sig_out[ , c("SNP","CHR","BP","codominant","dominant","recessive","overdominant","log-additive")]

# Zapisanie wyników
write.csv(sig_out,file='/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_out_0.05.csv', row.names=TRUE)


# R - SNPassoc - Associacja - Manhattan ----------------------------------------------------------------
# 
# Wczytywanie pliku
sig_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/sig_out_0.05.csv", header = TRUE, sep =',')


# Wczytanie bibliotek
library(ggplot2)
library(ggrepel) # Dodanie pakietu do obsługi etykiet

# Ustawienie wartości p-value
p_value <- 0.00001
threshold <- -log10(p_value)

# Stworzenie folderu "Wykresy_manhattan", jeśli jeszcze nie istnieje
if (!dir.exists("Wykresy_manhattan")) {
  dir.create("Wykresy_manhattan")
}

# Iteracja przez wybrane kolumny
for (x in names(sig_0.05)[5:9]) {
  # Przygotowanie danych do wykresu
  dane_wzor <- sig_0.05[1:5000, c("CHR", "BP", "SNP")]
  dane_wzor$p_value <- sig_0.05[1:5000, x]
  dane_wzor <- na.omit(dane_wzor)
  
  # Transformacja danych: -log10(p-value)
  dane_wzor$logP <- -log10(dane_wzor$p_value)
  
  # Filtracja punktów powyżej poziomu p-value
  dane_wzor$highlight <- dane_wzor$logP > threshold
  dane_wzor$label <- ifelse(dane_wzor$highlight, dane_wzor$SNP, NA)
  
  # Tworzenie wykresu Manhattan
  plot <- ggplot(dane_wzor, aes(x = BP, y = logP)) +
    geom_point(aes(color = highlight), size = 0.7, show.legend = FALSE) + # Punkty z kolorowaniem i bez legendy
    scale_color_manual(values = c("black", "blue")) + # Czarny dla punktów ponizej, niebieski dla powyzej
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") + # Linia p-value
    geom_text_repel(aes(label = label), size = 3, na.rm = TRUE, max.overlaps = Inf) + # Etykiety z unikanie kolizji
    theme_minimal(base_size = 14) + # Minimalistyczny styl
    theme(
      panel.grid.major = element_line(color = "grey"), # Dodanie gridu
      panel.grid.minor = element_line(color = "lightgrey"), 
      panel.background = element_rect(fill = "white"),
      plot.title = element_blank(), # Usuniecie tytułu
      axis.title.x = element_blank(), # Usuniecie tytułu osi X
      axis.text.x = element_blank(), # Usuniecie tekstu osi X
      axis.ticks.x = element_blank(), # Usuniecie znaczników osi X
      legend.position = "none" # Usuniecie legendy
    ) +
    labs(
      x = NULL, # Brak etykiety osi X
      y = expression(-log[10](italic(p))) # Etykieta osi Y
    )
  
  # Zapisywanie wykresu jako plik PNG w folderze "Wykresy_manhattan"
  ggsave(
    filename = file.path("Wykresy_manhattan", paste0("Manhattan_Plot_", x, ".png")),
    plot = plot,
    width = 12,
    height = 8
  )
}
# R - SNPassoc - Associacja cd ----------------------------------------------------------------

# Ustalenie progu istotnosci
p_value = 0.00001 # 1^-5

# Utworzenie pustej ramki danych
sig_dane_0.05 <- data.frame()

# Pęlta do filtracji SNP
for (x in names(sig_0.05)[5:9]) {
  dane_wzor = sig_0.05[1:50000,2:4]
  dane_wzor[x] = sig_0.05[1:5000,][x]
  dane_wzor = na.omit(dane_wzor)
  dane_wzor  <- dane_wzor[dane_wzor[,x] < p_value,]
  dane_wzor$model = x
  sig_dane_0.05 = bind_rows(sig_dane_0.05, dane_wzor)
}

# Wczytywanie danych genotypowych i mapy SNP
output_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/WYNIKI/SNPassoc/output_maf_0.05.csv", header = TRUE, sep =',')
inf_maf_0.05 <- read.table("/Users/daria/Desktop/1_semestr/Praca_magisterska/Plink_2/maf0.05.map", header = FALSE)

# Wyodrębnianie istotnych SNP
getSignificantSNPs(resHapMap_0.05,15, model = 'log')

# Przedstwienie wyników dla wybranych istotnych SNP
association(casecontrol~snp(TIGRP2P201649_rs8752112,sep=""), data=output_maf_0.05)
association(casecontrol~snp(BICF2S23334099,sep=""), data=output_maf_0.05)
association(casecontrol~snp(BICF2P352914,sep=""), data=output_maf_0.05)


# Terminal/Wiersz Polecen - Plink i AWK - bloki -----------------------------------------------------------

# Oblicza współczynniki LD (r²) dla SNP-ów z pliku selected_variants.txt
# filtrując wyniki dla SNP-ów o wartości r² powyzej 0.9
./plink --file chr15 --allow-no-sex --extract selected_variants.txt --r2 --ld-window-r2 0.9 --out ld_results

# Oblicza współczynniki LD (r²) dla SNP-ów z pliku selected_variants.txt bez żadnych filtrów na r²
./plink --file chr15 --allow-no-sex --extract selected_variants.txt --r2 --out ld_results_wszystkie

# Wyodrębnia liczbę bloków haplotypowych
awk '{print  NF-1}' plink.blocks > liczby_blokow.txt

# Zlicza liczbę wystąpień różnych liczebności bloków haplotypowych
awk '{counts[$1]++} END {for (number in counts) print number, counts[number]}' liczby_blokow.txt | sort -n > liczby_blokow_final.txt

# Wyodrębnia SNP-y o wartości r² większej niż 0.9 z pliku 
awk '$7 > 0.9 {print $3; print $6}' ld_results.ld | sort | uniq > snps_for_haplotypes.txt

# Tworzy bloki haplotypowe na podstawie SNP-ów 
./plink --file chr15 --extract  snps_for_haplotypes.txt --allow-no-sex --blocks 

# Przekształca genotypy SNP-ów na format BIMBAM
./plink --file chr15 --extract  snps_for_haplotypes.txt --recode bimbam --out bimbam

# Przekształca genotypy SNP-ów na format Haploview (do dodatkowej wizualizacji)
./plink --file chr15 --extract snps_for_haplotypes.txt --recodeHV --out haplo_to_view

# Zapisuje plik blocks do odpowiedniego formatu wykorzystanego w algorytmie EM
awk '{sub(/^\* /, ""); gsub(/ +/, ",", $0); sub(/,$/, ""); print}' plink.blocks > cleaned_blocks.csv



# R - haplotypy - Przygotowywanie Plikow -----------------------------------------------------------

# Zaladowanie bibliotek
library(data.table)
library(dplyr)

# Wczytanie pliku
bloki <- read.csv('/Users/daria/Desktop/Final/4_Estymacja/bimbam.recode.geno.txt',sep=",")

# Modyfikacja pliku
bloki = transpose(bloki)
colnames(bloki) <- unlist(bloki[1, ])
bloki <- bloki[-1, ] 

# Czyszczenie z ??
bloki <- bloki %>%
  mutate(across(everything(), ~na_if(.x, "??")))

# Dodawanie trait
bloki_2 <- read.csv('/Users/daria/Desktop/Final/4_Estymacja/bimbam.recode.pheno.txt',header = F, sep=",")
colnames(bloki_2) = 'trait'
result <- cbind(bloki_2, bloki)

# Dodawanie IND
bloki <- read.csv('/Users/daria/Desktop/Final/4_Estymacja/bimbam.recode.geno.txt',header = T, sep=",")
column_names <- names(bloki)
new_df <- data.frame(IND = column_names)
new_df <- new_df[-1, ] 
new_df = data_frame(new_df)
names(new_df) <- c("IND")

# Laczenie plikow w jeden
final <- cbind(new_df, result)

# Zapisywanie plikow
write.csv(final, "/Users/daria/Desktop/Final/4_Estymacja/tohaplo.csv")

bloki <- read.csv('/Users/daria/Desktop/Final/Haplo/tohaplo.csv',sep=",")
write.csv(bloki, "/Users/daria/Desktop/Final/Haplo/tohaplo2.csv.csv")

# R - haplotypy - Algorytm EM -----------------------------------------------------------

# Załadowanie bibliotek
library(stringr)
library(stringr)
library(SNPassoc)
library(haplo.stats)
library(data.table)

# Wczytanie danych
bloki <- read.csv('/Users/daria/Desktop/Final/4_Estymacja/cleaned_blocks.csv', header=F, sep="+")
cancer <- read.csv('/Users/daria/Desktop/Final/4_Estymacja/tohaplo2.csv', sep=",")

# Modyfikacja danych
cancer <- cancer[, -1] 
cancer <- cancer[, -1] 
cancer.s <- setupSNP(data=cancer, colSNPs=2:ncol(cancer), sep="")


# Sprawdzenie i filtrowanie SNP
valid_snps <- colnames(cancer)[-1]  # Wyodrębnia nazwy wszystkich kolumn SNP (pomijając kolumnę trait)

# Dla każdej linii w pliku `bloki` rozdziela SNP-y, filtruje te zgodne z danymi `cancer`, 
# i zwraca je jako ciąg oddzielony przecinkami
bloki <- apply(bloki, 1, function(row) {
  snps <- str_split_1(row, ", ")
  snps <- snps[snps %in% valid_snps]
  paste(snps, collapse=", ")
})

# Pętla przetwarząjca każdy blok SNP-ów
for (i in 1:length(bloki)) {
  wektor_snp <- c(str_split_1(bloki[i], ", "))
  
  if (length(wektor_snp) < 2) {
    cat("Blok ", i, " ma mniej niż dwa SNP. Pominięcie tego bloku.\n")
    next
  }
  
  snpsH <- wektor_snp
  genoH <- make.geno(cancer.s, snpsH)
  trait <- cancer.s$trait
  

  # Analiza haplotypów za pomocą modelu regresji logistycznej
  mod <- haplo.glm(trait ~ genoH,           
                   family="binomial", 
                   locus.label=snpsH,
                   allele.lev=attributes(genoH)$unique.alleles,
                   control = haplo.glm.control(haplo.freq.min=0.25))   
  
  # Przetwarzanie wyników
  df1 <- intervals(mod)
  
  if (nrow(df1) < 2) {
    cat("Blok ", i, " nie ma wystarczającej liczby haplotypów. Pominięcie tego bloku.\n")
    next
  }
  
  df1 <- df1[-1,c('freq', 'or', '95%', 'C.I.', 'P-val')]
  df1 <- as.data.frame(round(df1, 3))
  df1 <- df1[order(df1$freq, decreasing = T),] 
  
  df1$`P-val`[1] <- 1
  df1 <- df1[!(is.na(df1$`P-val`)),]
  
  rownames(df1) <- 1:nrow(df1)
  
  
  df2 <- as.data.frame(mod$haplo.unique)
  df2['freq'] = round(mod$haplo.freq, 3)
  df2 <- df2[order(df2$freq, decreasing = T),] 
  df2 <- subset(df2, select = -c(freq))
  rownames(df2) <- 1:nrow(df2)
  
  
  # Łaczenie tabel 
  t_w <- merge(df2, df1, by=0)
  t_w$`95%` <- paste(t_w$`95%`, t_w$C.I., sep=' - ')
  t_w <- subset(t_w, select = -c(Row.names, C.I.))
  
  # Zamiana nazw kolumn
  colnames(t_w)[which(names(t_w) == "95%")] <- "95% CI"
  colnames(t_w)[which(names(t_w) == "or")] <- "OR"
  colnames(t_w)[which(names(t_w) == "P-val")] <- "P"
  colnames(t_w)[which(names(t_w) == "freq")] <- "Frekwencja"
  
  # Usunięcie rzadkich haplotypow
  t_w <- t_w[!(is.na(t_w$Frekwencja) | t_w$Frekwencja <= 0.25),]
  
  if (nrow(t_w) == 0) {
    cat("Blok ", i, " nie ma wystarczającej liczby częstości haplotypów. Pominięcie tego bloku.\n")
    next
  }
  
  # Przygotowanie p-value do sprawdzenia w kolejnych krokach
  p_values <- as.vector(t_w$P)
  
  # Usunięcie przedziału i p-value dla haplotypu referencyjnego
  t_w$`95% CI`[1] <- '-'
  t_w$`P`[1] <- '-'
  
  
  # Zapisanie plików do odpowiednich folderów po sprawdzeniu istotności wyniku
  if(any(p_values <= 0.05)) {
    sink(file = str_glue("/Users/daria/Desktop/Final/przyjete/{i}przyjexte.txt"))
    print(t_w)
    sink(file = NULL)
  } else {
    sink(file = str_glue("/Users/daria/Desktop/Final/odrzucone/{i}odrzucone.txt"))
    print(t_w)
    sink(file = NULL)
  }
}

