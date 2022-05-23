#Calculate rates (TPR, FPR, TNR, FNR) for different cut-offs in Nanopore Data 
#VCF file with all samples (Nanopore and Illumina)
#allele_depth.txt file from joint VCF 

df <- read.table("allele_depth.txt",sep="\t",header=F,stringsAsFactors=F)
#allele_depth.txt looks like
#V1 V2  V3  V4  V5  V6  V7
#Chromosome 108 C A ERR221573 0/0 60,0
#With genotype and allele depth

df$V8 <- unlist(apply(df$V7,function(x) strsplit(x,",")[[1]][1]))
df$V9 <- unlist(apply(df$V7,function(x) strsplit(x,",")[[1]][2]))

#remove indels and multiallelic positions
poly <- grep(",+.+,",df$V7)
indels <- which(nchar(df$V3)>1|nchar(df$V4)>1) 
remove <- unique(c(poly,indels))
df2 <- df[-remove,]

df2$V9 <- as.numeric(df2$V9)
df2$V8 <- as.numeric(df2$V8)
df2$V10 <- df2$V9/(df2$V8+df2$V9)

#V10 is frequency, now stablish the -F at which this variant is called 
df2$F1 <- 0
df2$F2 <- 0
df2$F3 <- 0
df2$F4 <- 0
df2$F5 <- 0
df2$F6 <- 0
df2$F7 <- 0
df2$F8 <- 0
df2$F9 <- 0
df2$F10 <- 0

cols <- colnames(df2[11:20])

for(i in 1:nrow(df2)){
  if(df2$V10[i]!=0&df2$V10[i]!=1&!is.na(df2$V10[i])){
  x <- as.numeric(strsplit(strsplit(as.character(df2$V10[i]),"[.]")[[1]][2],"")[[1]][1])
    if(x==9){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
      df2$F4[i] <- 1
      df2$F5[i] <- 1
      df2$F6[i] <- 1
      df2$F7[i] <- 1
      df2$F8[i] <- 1
      df2$F9[i] <- 1
    }
    if(x==8){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
      df2$F4[i] <- 1
      df2$F5[i] <- 1
      df2$F6[i] <- 1
      df2$F7[i] <- 1
      df2$F8[i] <- 1
    }
    if(x==7){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
      df2$F4[i] <- 1
      df2$F5[i] <- 1
      df2$F6[i] <- 1
      df2$F7[i] <- 1
    }
    if(x==6){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
      df2$F4[i] <- 1
      df2$F5[i] <- 1
      df2$F6[i] <- 1
    }
    if(x==5){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
      df2$F4[i] <- 1
      df2$F5[i] <- 1
    }
    if(x==4){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
      df2$F4[i] <- 1
    }
    if(x==3){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
      df2$F3[i] <- 1
    }
    if(x==2){
      df2$F1[i] <- 1
      df2$F2[i] <- 1
    }
    if(x<1){
      df2$F1[i] <- 1
    }
  }
  if(df2$V10[i]==1&!is.na(df2$V10[i])){
    df2$F1[i] <- 1
    df2$F2[i] <- 1
    df2$F3[i] <- 1
    df2$F4[i] <- 1
    df2$F5[i] <- 1
    df2$F6[i] <- 1
    df2$F7[i] <- 1
    df2$F8[i] <- 1
    df2$F9[i] <- 1
  }
  if(is.na(df2$V10[i])){
    df2$F1[i] <- NA
    df2$F2[i] <- NA
    df2$F3[i] <- NA 
    df2$F4[i] <- NA
    df2$F5[i] <- NA
    df2$F6[i] <- NA
    df2$F7[i] <- NA 
    df2$F8[i] <- NA
    df2$F9[i] <- NA
    df2$F10[i] <- NA
  }
}

illu <- df2[grep("ERR",df2$V5),]
nano <- df2[grep("barcode",df2$V5),]
illu2 <- illu
nano2 <- nano
meta <- read.table("./metadata.txt",sep="\t",header=T,stringsAsFactors=F) #metadata file with
#SAMPLE_ID, RUN_ID, BARCODE, ACCESSION, LINEAGE, SUBLINEAGE
#BARCODE = nanopore ID, ACCESSION = Illumina ID

for(i in 1:nrow(illu2)){
  id <- which(meta$ACCESSION==illu2$V5[i])
  illu2[i,"V5"] <- meta[id,"BARCODE"]
}
idx2 <- c()
for(i in 1:nrow(nano2)){
  if(length(which(illu2$V2==nano2$V2[i]&illu2$V5==nano2$V5[i]))>0){
  idx2[i] <-which(illu2$V2==nano2$V2[i]&illu2$V5==nano2$V5[i])
  }else{
  idx2[i] <- NA
  }
}

idx2 <- na.omit(idx2)
illu2 <- illu2[idx2,]

library(reshape2)
illu2_melt <- melt(illu2[11:20])

idx <- which(illu2_melt$value==1)
illu2_melt$group <- "no_variant"
illu2_melt$group[idx] <- "variant"


nano2_melt <- melt(nano2[11:20])
nano2_melt$group <- illu2_melt$group

colnames(nano2_melt) <- c("F","value","group")
nano2_melt <- nano2_melt[c(3,2,1)]
nano2_melt$F <- as.character(nano2_melt$F)
idx <- which(nano2_melt$F=="F1")
nano2_melt[idx,"F"] <- 0.1
idx <- which(nano2_melt$F=="F2")
nano2_melt[idx,"F"] <- 0.2
idx <- which(nano2_melt$F=="F3")
nano2_melt[idx,"F"] <- 0.3
idx <- which(nano2_melt$F=="F4")
nano2_melt[idx,"F"] <- 0.4
idx <- which(nano2_melt$F=="F5")
nano2_melt[idx,"F"] <- 0.5
idx <- which(nano2_melt$F=="F6")
nano2_melt[idx,"F"] <- 0.6
idx <- which(nano2_melt$F=="F7")
nano2_melt[idx,"F"] <- 0.7
idx <- which(nano2_melt$F=="F8")
nano2_melt[idx,"F"] <- 0.8
idx <- which(nano2_melt$F=="F9")
nano2_melt[idx,"F"] <- 0.9
idx <- which(nano2_melt$F=="F10")
nano2_melt[idx,"F"] <- 1

colnames(illu2_melt)[1] <- "F"
illu2_melt$F <- as.character(illu2_melt$F)
idx <- which(illu2_melt$F=="F1")
illu2_melt[idx,"F"] <- 0.1
idx <- which(illu2_melt$F=="F2")
illu2_melt[idx,"F"] <- 0.2
idx <- which(illu2_melt$F=="F3")
illu2_melt[idx,"F"] <- 0.3
idx <- which(illu2_melt$F=="F4")
illu2_melt[idx,"F"] <- 0.4
idx <- which(illu2_melt$F=="F5")
illu2_melt[idx,"F"] <- 0.5
idx <- which(illu2_melt$F=="F6")
illu2_melt[idx,"F"] <- 0.6
idx <- which(illu2_melt$F=="F7")
illu2_melt[idx,"F"] <- 0.7
idx <- which(illu2_melt$F=="F8")
illu2_melt[idx,"F"] <- 0.8
idx <- which(illu2_melt$F=="F9")
illu2_melt[idx,"F"] <- 0.9
idx <- which(illu2_melt$F=="F10")
illu2_melt[idx,"F"] <- 1

illu2_melt$pos <- rep(rep(unique(illu2$V2),each=10),10)
nano2_melt$pos <- rep(rep(unique(nano2$V2),each=10),10)
illu2_melt$sample <- rep(rep(unique(illu2$V5),10),length(unique(illu2$V2)))
nano2_melt$sample <- rep(rep(unique(nano2$V5),10),length(unique(nano2$V2)))


rates <- data.frame(matrix(ncol=4,nrow=10))
colnames(rates) <- c("TP","FP","TN","FN")
rownames(rates) <- c("F01","F02","F03","F04","F05","F06","F07","F08","F09","F1")

#Change variant or no variant based on fixed for illumina at least until 0.7

idx <- which(nano2_melt$F<0.7)
for(i in 1:length(idx)){
  id <- which(illu2_melt$pos==nano2_melt$pos[i]&illu2_melt$sample==nano2_melt$sample[i]&illu2_melt$F==0.7)
  nano2_melt$group[i] <- illu2_melt[id,"group"]
}

rates[1,1] <- length(which(nano2_melt[which(nano2_melt$F==0.1),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.1),"value"]==1))

rates[2,1] <- length(which(nano2_melt[which(nano2_melt$F==0.2),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.2),"value"]==1))

rates[3,1] <- length(which(nano2_melt[which(nano2_melt$F==0.3),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.3),"value"]==1))

rates[4,1] <- length(which(nano2_melt[which(nano2_melt$F==0.4),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.4),"value"]==1))

rates[5,1] <- length(which(nano2_melt[which(nano2_melt$F==0.5),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.5),"value"]==1))

rates[6,1] <- length(which(nano2_melt[which(nano2_melt$F==0.6),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.6),"value"]==1))

rates[7,1] <- length(which(nano2_melt[which(nano2_melt$F==0.7),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.7),"value"]==1))

rates[8,1] <- length(which(nano2_melt[which(nano2_melt$F==0.8),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.8),"value"]==1))

rates[9,1] <- length(which(nano2_melt[which(nano2_melt$F==0.9),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.9),"value"]==1))

rates[10,1] <- length(which(nano2_melt[which(nano2_melt$F==1),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==1),"value"]==1))
rates[1,2] <- length(which(nano2_melt[which(nano2_melt$F==0.1),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.1),"value"]==1))

rates[2,2] <- length(which(nano2_melt[which(nano2_melt$F==0.2),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.2),"value"]==1))

rates[3,2] <- length(which(nano2_melt[which(nano2_melt$F==0.3),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.3),"value"]==1))

rates[4,2] <- length(which(nano2_melt[which(nano2_melt$F==0.4),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.4),"value"]==1))

rates[5,2] <- length(which(nano2_melt[which(nano2_melt$F==0.5),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.5),"value"]==1))

rates[6,2] <- length(which(nano2_melt[which(nano2_melt$F==0.6),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.6),"value"]==1))

rates[7,2] <- length(which(nano2_melt[which(nano2_melt$F==0.7),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.7),"value"]==1))

rates[8,2] <- length(which(nano2_melt[which(nano2_melt$F==0.8),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.8),"value"]==1))

rates[9,2] <- length(which(nano2_melt[which(nano2_melt$F==0.9),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.9),"value"]==1))

rates[10,2] <- length(which(nano2_melt[which(nano2_melt$F==1),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==1),"value"]==1))

rates[1,3] <- length(which(nano2_melt[which(nano2_melt$F==0.1),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.1),"value"]==0))

rates[2,3] <- length(which(nano2_melt[which(nano2_melt$F==0.2),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.2),"value"]==0))

rates[3,3] <- length(which(nano2_melt[which(nano2_melt$F==0.3),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.3),"value"]==0))

rates[4,3] <- length(which(nano2_melt[which(nano2_melt$F==0.4),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.4),"value"]==0))

rates[5,3] <- length(which(nano2_melt[which(nano2_melt$F==0.5),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.5),"value"]==0))

rates[6,3] <- length(which(nano2_melt[which(nano2_melt$F==0.6),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.6),"value"]==0))

rates[7,3] <- length(which(nano2_melt[which(nano2_melt$F==0.7),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.7),"value"]==0))

rates[8,3] <- length(which(nano2_melt[which(nano2_melt$F==0.8),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.8),"value"]==0))

rates[9,3] <- length(which(nano2_melt[which(nano2_melt$F==0.9),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==0.9),"value"]==0))

rates[10,3] <- length(which(nano2_melt[which(nano2_melt$F==1),"group"]=="no_variant"&
nano2_melt[which(nano2_melt$F==1),"value"]==0))

rates[1,4] <- length(which(nano2_melt[which(nano2_melt$F==0.1),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.1),"value"]==0))

rates[2,4] <- length(which(nano2_melt[which(nano2_melt$F==0.2),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.2),"value"]==0))

rates[3,4] <- length(which(nano2_melt[which(nano2_melt$F==0.3),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.3),"value"]==0))

rates[4,4] <- length(which(nano2_melt[which(nano2_melt$F==0.4),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.4),"value"]==0))

rates[5,4] <- length(which(nano2_melt[which(nano2_melt$F==0.5),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.5),"value"]==0))

rates[6,4] <- length(which(nano2_melt[which(nano2_melt$F==0.6),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.6),"value"]==0))

rates[7,4] <- length(which(nano2_melt[which(nano2_melt$F==0.7),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.7),"value"]==0))

rates[8,4] <- length(which(nano2_melt[which(nano2_melt$F==0.8),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.8),"value"]==0))

rates[9,4] <- length(which(nano2_melt[which(nano2_melt$F==0.9),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==0.9),"value"]==0))

rates[10,4] <- length(which(nano2_melt[which(nano2_melt$F==1),"group"]=="variant"&
nano2_melt[which(nano2_melt$F==1),"value"]==0))


#TPR = TP / TP + FN
#TNR = TN / TN + FP

#FPR = FP / FP + TN
#FNR = FN / FN + TP
#PPV = TP / TP + FP

rates$TPR <- rates$TP / (rates$TP + rates$FN)
rates$TNR <- rates$TN / (rates$TN + rates$FP)
rates$FPR <- rates$FP / (rates$FP + rates$TN)
rates$FNR <- rates$FN / (rates$FN + rates$TP)
rates$PPV <- rates$TP / (rates$TP + rates$FP)



nano2_melt$F <- as.numeric(nano2_melt$F)

rates$F <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
