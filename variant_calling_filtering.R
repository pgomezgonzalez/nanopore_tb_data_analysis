#Pairs of samples extracted from combined VCF file

#####pairs extracted, called pair1.txt 
head1 <- c("CHROM","POS","REF","ALT","GT_nanopore","DP_nanopore","AD_nanopore","GT_illumina","DP_illumina","AD_illumina")
pair1 <- read.table("pair1.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair1) <- head1
pair1[11] <- NULL

pair1$AD1_ref <- unlist(lapply(pair1$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair1$AD1_alt <- unlist(lapply(pair1$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair1$AD1 <- as.numeric(pair1$AD1_alt)/(as.numeric(pair1$AD1_ref)+as.numeric(pair1$AD1_alt))
pair1$AD2_ref <- unlist(lapply(pair1$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair1$AD2_alt <- unlist(lapply(pair1$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair1$AD2 <- as.numeric(pair1$AD2_alt)/(as.numeric(pair1$AD2_ref)+as.numeric(pair1$AD2_alt))
pair1$DP_nanopore <- as.numeric(pair1$DP_nanopore)
pair1$DP_illumina <- as.numeric(pair1$DP_illumina)

##Filter by DP 
idx <- which(pair1$DP_nanopore<10|pair1$DP_illumina<10)
pair1_filtered <- pair1[-idx,]
#filter by AD 
#error in 1765 
pair1_filtered <- pair1_filtered[-1765,]
for(i in 1:nrow(pair1_filtered)){
  if(pair1_filtered$AD1[i]<0.7){
    if(pair1_filtered[i,"GT_nanopore"]=="1/1"|pair1_filtered[i,"GT_nanopore"]=="0/1"){
      pair1_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair1_filtered[i,"GT_nanopore"] <- pair1_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair1_filtered[i,"GT_nanopore"]=="0/1"|pair1_filtered[i,"GT_nanopore"]=="0/0"){
      pair1_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair1_filtered[i,"GT_nanopore"] <- pair1_filtered[i,"GT_nanopore"]
    }
  }
}
#error in 165 NA 
#no covered in illumina so delete 
pair1_filtered <- pair1_filtered[-165,]
for(i in 1:nrow(pair1_filtered)){
  if(pair1_filtered$AD2[i]<0.7){
    if(pair1_filtered[i,"GT_illumina"]=="1/1"|pair1_filtered[i,"GT_illumina"]=="0/1"){
      pair1_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair1_filtered[i,"GT_illumina"] <- pair1_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair1_filtered[i,"GT_illumina"]=="0/1"|pair1_filtered[i,"GT_illumina"]=="0/0"){
      pair1_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair1_filtered[i,"GT_illumina"] <- pair1_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair1_filtered$GT_nanopore)
if(length(idx)>0){
  pair1_filtered <- pair1_filtered[-idx,]
}
idx <- grep("[.]",pair1_filtered$GT_illumina)
if(length(idx)>0){
  pair1_filtered <- pair1_filtered[-idx,]
}
nrow(pair1_filtered) #3938

pair1_filtered$callNP <- unlist(lapply(pair1_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair1_filtered$callIL <- unlist(lapply(pair1_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair1_filtered$GT_nanopore)
unique(pair1_filtered$GT_illumina)
##Fix 2/2 and 1/2 2/1 
idx <- grep("2/2",pair1_filtered$GT_nanopore)
pair1_filtered[idx,]  #it is an indel, pos 2143327, so remove 
pair1_filtered <- pair1_filtered[-idx,]
idx <- grep("1/2",pair1_filtered$GT_nanopore)
pair1_filtered[idx,]
#AD doens't reach 0.7 in any allele, is POS 55553 discrepancy between nanopore and illumina
pair1_filtered[idx,"callNP"] <- 0

pair1_filtered$callNP <- as.numeric(pair1_filtered$callNP)
pair1_filtered$callIL <- as.numeric(pair1_filtered$callIL)

pair1upset <- upset(pair1_filtered[17:18])

table(c(pair1_filtered[c("callNP","callIL")]))
pair1_filtered[which(pair1_filtered$callNP==0&pair1_filtered$callIL==1),]


#pair2 
pair2 <- read.table("pair2.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair2) <- head1
pair2[11] <- NULL

pair2$AD1_ref <- unlist(lapply(pair2$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair2$AD1_alt <- unlist(lapply(pair2$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair2$AD1 <- as.numeric(pair2$AD1_alt)/(as.numeric(pair2$AD1_ref)+as.numeric(pair2$AD1_alt))
pair2$AD2_ref <- unlist(lapply(pair2$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair2$AD2_alt <- unlist(lapply(pair2$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair2$AD2 <- as.numeric(pair2$AD2_alt)/(as.numeric(pair2$AD2_ref)+as.numeric(pair2$AD2_alt))
pair2$DP_nanopore <- as.numeric(pair2$DP_nanopore)
pair2$DP_illumina <- as.numeric(pair2$DP_illumina)

##Filter by DP 
idx <- which(pair2$DP_nanopore<10|pair2$DP_illumina<10)
pair2_filtered <- pair2[-idx,]
#filter by AD 
#pos 3378828 (3031) remove because is deleted in both 
pair2_filtered <- pair2_filtered[-3031,]
for(i in 1:nrow(pair2_filtered)){
  if(pair2_filtered$AD1[i]<0.7){
    if(pair2_filtered[i,"GT_nanopore"]=="1/1"|pair2_filtered[i,"GT_nanopore"]=="0/1"){
      pair2_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair2_filtered[i,"GT_nanopore"] <- pair2_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair2_filtered[i,"GT_nanopore"]=="0/1"|pair2_filtered[i,"GT_nanopore"]=="0/0"){
      pair2_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair2_filtered[i,"GT_nanopore"] <- pair2_filtered[i,"GT_nanopore"]
    }
  }
}
#868 (pos 850206) deleted in Illumina 
pair2_filtered<-  pair2_filtered[-868,]
pair2_filtered<-  pair2_filtered[-868,]

for(i in 1:nrow(pair2_filtered)){
  if(pair2_filtered$AD2[i]<0.7){
    if(pair2_filtered[i,"GT_illumina"]=="1/1"|pair2_filtered[i,"GT_illumina"]=="0/1"){
      pair2_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair2_filtered[i,"GT_illumina"] <- pair2_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair2_filtered[i,"GT_illumina"]=="0/1"|pair2_filtered[i,"GT_illumina"]=="0/0"){
      pair2_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair2_filtered[i,"GT_illumina"] <- pair2_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair2_filtered$GT_nanopore)
if(length(idx)>0){
  pair2_filtered <- pair2_filtered[-idx,]
}
idx <- grep("[.]",pair2_filtered$GT_illumina)
if(length(idx)>0){
  pair2_filtered <- pair2_filtered[-idx,]
}
nrow(pair2_filtered) #3948

pair2_filtered$callNP <- unlist(lapply(pair2_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair2_filtered$callIL <- unlist(lapply(pair2_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair2_filtered$GT_nanopore)
unique(pair2_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("2/2",pair2_filtered$GT_nanopore)
pair2_filtered[idx,]  ##it is a deletion 
#if both nano and illu 2/2 <- considered it 1 and 1  pos 2143327, but is an indel, so delete 
pair2_filtered <- pair2_filtered[-idx,]
idx <- grep("1/2",pair2_filtered$GT_nanopore)
pair2_filtered[idx,]
#AD doens't reach 0.7 in any allele, pos 55553, same 
pair2_filtered[idx,"callNP"] <- 0

pair2_filtered$callNP <- as.numeric(pair2_filtered$callNP)
pair2_filtered$callIL <- as.numeric(pair2_filtered$callIL)

pair2upset <- upset(pair2_filtered[17:18])

table(c(pair2_filtered[c("callNP","callIL")]))
pair2_filtered[which(pair2_filtered$callNP==0&pair2_filtered$callIL==1),]

#pair3
pair3 <- read.table("pair3.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair3) <- head1
pair3[11] <- NULL

pair3$AD1_ref <- unlist(lapply(pair3$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair3$AD1_alt <- unlist(lapply(pair3$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair3$AD1 <- as.numeric(pair3$AD1_alt)/(as.numeric(pair3$AD1_ref)+as.numeric(pair3$AD1_alt))
pair3$AD2_ref <- unlist(lapply(pair3$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair3$AD2_alt <- unlist(lapply(pair3$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair3$AD2 <- as.numeric(pair3$AD2_alt)/(as.numeric(pair3$AD2_ref)+as.numeric(pair3$AD2_alt))
pair3$DP_nanopore <- as.numeric(pair3$DP_nanopore)
pair3$DP_illumina <- as.numeric(pair3$DP_illumina)

##Filter by DP 
idx <- which(pair3$DP_nanopore<10|pair3$DP_illumina<10)
pair3_filtered <- pair3[-idx,]
#filter by AD 
#3029, pos 3378828, not covered, delete 
pair3_filtered <- pair3_filtered[-3029,]
for(i in 1:nrow(pair3_filtered)){
  if(pair3_filtered$AD1[i]<0.7){
    if(pair3_filtered[i,"GT_nanopore"]=="1/1"|pair3_filtered[i,"GT_nanopore"]=="0/1"){
      pair3_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair3_filtered[i,"GT_nanopore"] <- pair3_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair3_filtered[i,"GT_nanopore"]=="0/1"|pair3_filtered[i,"GT_nanopore"]=="0/0"){
      pair3_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair3_filtered[i,"GT_nanopore"] <- pair3_filtered[i,"GT_nanopore"]
    }
  }
}
#867, remove 
pair3_filtered <- pair3_filtered[-867,]
pair3_filtered <- pair3_filtered[-867,]
pair3_filtered <- pair3_filtered[-3026,]
for(i in 1:nrow(pair3_filtered)){
  if(pair3_filtered$AD2[i]<0.7){
    if(pair3_filtered[i,"GT_illumina"]=="1/1"|pair3_filtered[i,"GT_illumina"]=="0/1"){
      pair3_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair3_filtered[i,"GT_illumina"] <- pair3_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair3_filtered[i,"GT_illumina"]=="0/1"|pair3_filtered[i,"GT_illumina"]=="0/0"){
      pair3_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair3_filtered[i,"GT_illumina"] <- pair3_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair3_filtered$GT_nanopore)
if(length(idx)>0){
  pair3_filtered <- pair3_filtered[-idx,]
}
idx <- grep("[.]",pair3_filtered$GT_illumina)
if(length(idx)>0){
  pair3_filtered <- pair3_filtered[-idx,]
}
nrow(pair3_filtered) #3940

pair3_filtered$callNP <- unlist(lapply(pair3_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair3_filtered$callIL <- unlist(lapply(pair3_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair3_filtered$GT_nanopore)
unique(pair3_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("2/2",pair3_filtered$GT_nanopore)
pair3_filtered[idx,]  ##it is a deletion, pos 2143327
pair3_filtered <- pair3_filtered[-idx,]
idx <- grep("1/2",pair3_filtered$GT_nanopore)
pair3_filtered[idx,]
#AD doens't reach 0.7 in any allele, pos 55553
pair3_filtered[idx,"callNP"] <- 0

pair3_filtered$callNP <- as.numeric(pair3_filtered$callNP)
pair3_filtered$callIL <- as.numeric(pair3_filtered$callIL)

pair3upset <- upset(pair3_filtered[17:18])

table(c(pair3_filtered[c("callNP","callIL")]))
pair3_filtered[which(pair3_filtered$callNP==0&pair3_filtered$callIL==1),]

#pair4
pair4 <- read.table("pair4.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair4) <- head1
pair4[11] <- NULL

pair4$AD1_ref <- unlist(lapply(pair4$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair4$AD1_alt <- unlist(lapply(pair4$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair4$AD1 <- as.numeric(pair4$AD1_alt)/(as.numeric(pair4$AD1_ref)+as.numeric(pair4$AD1_alt))
pair4$AD2_ref <- unlist(lapply(pair4$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair4$AD2_alt <- unlist(lapply(pair4$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair4$AD2 <- as.numeric(pair4$AD2_alt)/(as.numeric(pair4$AD2_ref)+as.numeric(pair4$AD2_alt))
pair4$DP_nanopore <- as.numeric(pair4$DP_nanopore)
pair4$DP_illumina <- as.numeric(pair4$DP_illumina)

##Filter by DP 
idx <- which(pair4$DP_nanopore<10|pair4$DP_illumina<10)
pair4_filtered <- pair4[-idx,]
#filter by AD
for(i in 1:nrow(pair4_filtered)){
  if(pair4_filtered$AD1[i]<0.7){
    if(pair4_filtered[i,"GT_nanopore"]=="1/1"|pair4_filtered[i,"GT_nanopore"]=="0/1"){
      pair4_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair4_filtered[i,"GT_nanopore"] <- pair4_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair4_filtered[i,"GT_nanopore"]=="0/1"|pair4_filtered[i,"GT_nanopore"]=="0/0"){
      pair4_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair4_filtered[i,"GT_nanopore"] <- pair4_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair4_filtered)){
  if(pair4_filtered$AD2[i]<0.7){
    if(pair4_filtered[i,"GT_illumina"]=="1/1"|pair4_filtered[i,"GT_illumina"]=="0/1"){
      pair4_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair4_filtered[i,"GT_illumina"] <- pair4_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair4_filtered[i,"GT_illumina"]=="0/1"|pair4_filtered[i,"GT_illumina"]=="0/0"){
      pair4_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair4_filtered[i,"GT_illumina"] <- pair4_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair4_filtered$GT_nanopore)
if(length(idx)>0){
  pair4_filtered <- pair4_filtered[-idx,]
}
idx <- grep("[.]",pair4_filtered$GT_illumina)
if(length(idx)>0){
  pair4_filtered <- pair4_filtered[-idx,]
}
nrow(pair4_filtered) #3945
pair4_filtered$callNP <- unlist(lapply(pair4_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair4_filtered$callIL <- unlist(lapply(pair4_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair4_filtered$GT_nanopore)
unique(pair4_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("2/2",pair4_filtered$GT_nanopore)
pair4_filtered[idx,]  ##it is a deletion 
pair4_filtered <- pair4_filtered[-idx,]
idx <- grep("1/2",pair4_filtered$GT_nanopore)
pair4_filtered[idx,]
#AD doens't reach 0.7 in any allele, 55553
pair4_filtered[idx,"callNP"] <- 0

pair4_filtered$callNP <- as.numeric(pair4_filtered$callNP)
pair4_filtered$callIL <- as.numeric(pair4_filtered$callIL)

pair4upset <- upset(pair4_filtered[17:18])

table(c(pair4_filtered[c("callNP","callIL")]))
pair4_filtered[which(pair4_filtered$callNP==0&pair4_filtered$callIL==1),]


#pair5
pair5 <- read.table("pair5.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair5) <- head1
pair5[11] <- NULL

pair5$AD1_ref <- unlist(lapply(pair5$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair5$AD1_alt <- unlist(lapply(pair5$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair5$AD1 <- as.numeric(pair5$AD1_alt)/(as.numeric(pair5$AD1_ref)+as.numeric(pair5$AD1_alt))
pair5$AD2_ref <- unlist(lapply(pair5$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair5$AD2_alt <- unlist(lapply(pair5$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair5$AD2 <- as.numeric(pair5$AD2_alt)/(as.numeric(pair5$AD2_ref)+as.numeric(pair5$AD2_alt))
pair5$DP_nanopore <- as.numeric(pair5$DP_nanopore)
pair5$DP_illumina <- as.numeric(pair5$DP_illumina)


##Filter by DP 
idx <- which(pair5$DP_nanopore<10|pair5$DP_illumina<10)
if(length(idx)>0){
pair5_filtered <- pair5[-idx,]
}else{
pair5_filtered <- pair5
}
#filter by AD 
for(i in 1:nrow(pair5_filtered)){
  if(pair5_filtered$AD1[i]<0.7){
    if(pair5_filtered[i,"GT_nanopore"]=="1/1"|pair5_filtered[i,"GT_nanopore"]=="0/1"){
      pair5_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair5_filtered[i,"GT_nanopore"] <- pair5_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair5_filtered[i,"GT_nanopore"]=="0/1"|pair5_filtered[i,"GT_nanopore"]=="0/0"){
      pair5_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair5_filtered[i,"GT_nanopore"] <- pair5_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair5_filtered)){
  if(pair5_filtered$AD2[i]<0.7){
    if(pair5_filtered[i,"GT_illumina"]=="1/1"|pair5_filtered[i,"GT_illumina"]=="0/1"){
      pair5_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair5_filtered[i,"GT_illumina"] <- pair5_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair5_filtered[i,"GT_illumina"]=="0/1"|pair5_filtered[i,"GT_illumina"]=="0/0"){
      pair5_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair5_filtered[i,"GT_illumina"] <- pair5_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair5_filtered$GT_nanopore)
if(length(idx)>0){
  pair5_filtered <- pair5_filtered[-idx,]
}
idx <- grep("[.]",pair5_filtered$GT_illumina)
if(length(idx)>0){
  pair5_filtered <- pair5_filtered[-idx,]
}
nrow(pair5_filtered) #3957

pair5_filtered$callNP <- unlist(lapply(pair5_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair5_filtered$callIL <- unlist(lapply(pair5_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair5_filtered$GT_nanopore)
unique(pair5_filtered$GT_illumina)

pair5_filtered$callNP <- as.numeric(pair5_filtered$callNP)
pair5_filtered$callIL <- as.numeric(pair5_filtered$callIL)

pair5upset <- upset(pair5_filtered[17:18])

table(c(pair5_filtered[c("callNP","callIL")]))



#pair6
pair6 <- read.table("pair6.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair6) <- head1
pair6[11] <- NULL

pair6$AD1_ref <- unlist(lapply(pair6$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair6$AD1_alt <- unlist(lapply(pair6$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair6$AD1 <- as.numeric(pair6$AD1_alt)/(as.numeric(pair6$AD1_ref)+as.numeric(pair6$AD1_alt))
pair6$AD2_ref <- unlist(lapply(pair6$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair6$AD2_alt <- unlist(lapply(pair6$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair6$AD2 <- as.numeric(pair6$AD2_alt)/(as.numeric(pair6$AD2_ref)+as.numeric(pair6$AD2_alt))
pair6$DP_nanopore <- as.numeric(pair6$DP_nanopore)
pair6$DP_illumina <- as.numeric(pair6$DP_illumina)


##Filter by DP 
idx <- which(pair6$DP_nanopore<10|pair6$DP_illumina<10)
if(length(idx)>0){
pair6_filtered <- pair6[-idx,]
}else{
pair6_filtered <- pair6
}
#filter by AD 
for(i in 1:nrow(pair6_filtered)){
  if(pair6_filtered$AD1[i]<0.7){
    if(pair6_filtered[i,"GT_nanopore"]=="1/1"|pair6_filtered[i,"GT_nanopore"]=="0/1"){
      pair6_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair6_filtered[i,"GT_nanopore"] <- pair6_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair6_filtered[i,"GT_nanopore"]=="0/1"|pair6_filtered[i,"GT_nanopore"]=="0/0"){
      pair6_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair6_filtered[i,"GT_nanopore"] <- pair6_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair6_filtered)){
  if(pair6_filtered$AD2[i]<0.7){
    if(pair6_filtered[i,"GT_illumina"]=="1/1"|pair6_filtered[i,"GT_illumina"]=="0/1"){
      pair6_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair6_filtered[i,"GT_illumina"] <- pair6_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair6_filtered[i,"GT_illumina"]=="0/1"|pair6_filtered[i,"GT_illumina"]=="0/0"){
      pair6_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair6_filtered[i,"GT_illumina"] <- pair6_filtered[i,"GT_illumina"]
    }
  }
}
idx<- grep("[.]",pair6_filtered$GT_nanopore)
if(length(idx)>0){
  pair6_filtered <- pair6_filtered[-idx,]
}
idx <- grep("[.]",pair6_filtered$GT_illumina)
if(length(idx)>0){
  pair6_filtered <- pair6_filtered[-idx,]
}
nrow(pair6_filtered) #3914

pair6_filtered$callNP <- unlist(lapply(pair6_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair6_filtered$callIL <- unlist(lapply(pair6_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair6_filtered$GT_nanopore)
unique(pair6_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair6_filtered$GT_nanopore)
pair6_filtered[idx,]  ##it is a deletion 

#Pos 55553 is a deletion in nanopore and SNP in illumina, so 0 1 
#Pos 2660319 is a SNP, second alternate in both nano and illumina, so 1 1 
pair6_filtered[idx,"callNP"] <- c(0,1)
pair6_filtered[idx,"callIL"] <- c(1,1)


pair6_filtered$callNP <- as.numeric(pair6_filtered$callNP)
pair6_filtered$callIL <- as.numeric(pair6_filtered$callIL)

pair6upset <- upset(pair6_filtered[17:18])

table(c(pair6_filtered[c("callNP","callIL")]))

pair6_filtered[which(pair6_filtered$callNP==0&pair6_filtered$callIL==1),]


#pair7
pair7 <- read.table("pair7.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair7) <- head1
pair7[11] <- NULL

pair7$AD1_ref <- unlist(lapply(pair7$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair7$AD1_alt <- unlist(lapply(pair7$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair7$AD1 <- as.numeric(pair7$AD1_alt)/(as.numeric(pair7$AD1_ref)+as.numeric(pair7$AD1_alt))
pair7$AD2_ref <- unlist(lapply(pair7$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair7$AD2_alt <- unlist(lapply(pair7$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair7$AD2 <- as.numeric(pair7$AD2_alt)/(as.numeric(pair7$AD2_ref)+as.numeric(pair7$AD2_alt))
pair7$DP_nanopore <- as.numeric(pair7$DP_nanopore)
pair7$DP_illumina <- as.numeric(pair7$DP_illumina)

##Filter by DP 
idx <- which(pair7$DP_nanopore<10|pair7$DP_illumina<10)
if(length(idx)>0){
pair7_filtered <- pair7[-idx,]
}else{
pair7_filtered <- pair7
}
#filter by AD 
for(i in 1:nrow(pair7_filtered)){
  if(pair7_filtered$AD1[i]<0.7){
    if(pair7_filtered[i,"GT_nanopore"]=="1/1"|pair7_filtered[i,"GT_nanopore"]=="0/1"){
      pair7_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair7_filtered[i,"GT_nanopore"] <- pair7_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair7_filtered[i,"GT_nanopore"]=="0/1"|pair7_filtered[i,"GT_nanopore"]=="0/0"){
      pair7_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair7_filtered[i,"GT_nanopore"] <- pair7_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair7_filtered)){
  if(pair7_filtered$AD2[i]<0.7){
    if(pair7_filtered[i,"GT_illumina"]=="1/1"|pair7_filtered[i,"GT_illumina"]=="0/1"){
      pair7_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair7_filtered[i,"GT_illumina"] <- pair7_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair7_filtered[i,"GT_illumina"]=="0/1"|pair7_filtered[i,"GT_illumina"]=="0/0"){
      pair7_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair7_filtered[i,"GT_illumina"] <- pair7_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair7_filtered$GT_nanopore)
if(length(idx)>0){
  pair7_filtered <- pair7_filtered[-idx,]
}
idx <- grep("[.]",pair7_filtered$GT_illumina)
if(length(idx)>0){
  pair7_filtered <- pair7_filtered[-idx,]
}
nrow(pair7_filtered) #3820

pair7_filtered$callNP <- unlist(lapply(pair7_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair7_filtered$callIL <- unlist(lapply(pair7_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair7_filtered$GT_nanopore)
unique(pair7_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair7_filtered$GT_nanopore)
pair7_filtered[idx,]  ##it is a deletion 
#Pos 26600319 is a SNP, second alternate in both nano and illumina, so 1 1 
pair7_filtered[idx,c("callNP","callIL")] <- c(1,1)

idx <- grep("1/2",pair7_filtered$GT_nanopore)
pair7_filtered[idx,]
#pos 55553 
#AD doens't reach 0.7 in any allele 
pair7_filtered[idx,"callNP"] <- 0

pair7_filtered$callNP <- as.numeric(pair7_filtered$callNP)
pair7_filtered$callIL <- as.numeric(pair7_filtered$callIL)

pair7upset <- upset(pair7_filtered[17:18])

table(c(pair7_filtered[c("callNP","callIL")]))
pair7_filtered[which(pair7_filtered$callNP==0&pair7_filtered$callIL==1),]

#pair8
pair8 <- read.table("pair8.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair8) <- head1
pair8[11] <- NULL

pair8$AD1_ref <- unlist(lapply(pair8$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair8$AD1_alt <- unlist(lapply(pair8$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair8$AD1 <- as.numeric(pair8$AD1_alt)/(as.numeric(pair8$AD1_ref)+as.numeric(pair8$AD1_alt))
pair8$AD2_ref <- unlist(lapply(pair8$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair8$AD2_alt <- unlist(lapply(pair8$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair8$AD2 <- as.numeric(pair8$AD2_alt)/(as.numeric(pair8$AD2_ref)+as.numeric(pair8$AD2_alt))
pair8$DP_nanopore <- as.numeric(pair8$DP_nanopore)
pair8$DP_illumina <- as.numeric(pair8$DP_illumina)

##Filter by DP 
idx <- which(pair8$DP_nanopore<10|pair8$DP_illumina<10)
if(length(idx)>0){
pair8_filtered <- pair8[-idx,]
}else{
pair8_filtered <- pair8
}
#filter by AD 
for(i in 1:nrow(pair8_filtered)){
  if(pair8_filtered$AD1[i]<0.7){
    if(pair8_filtered[i,"GT_nanopore"]=="1/1"|pair8_filtered[i,"GT_nanopore"]=="0/1"){
      pair8_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair8_filtered[i,"GT_nanopore"] <- pair8_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair8_filtered[i,"GT_nanopore"]=="0/1"|pair8_filtered[i,"GT_nanopore"]=="0/0"){
      pair8_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair8_filtered[i,"GT_nanopore"] <- pair8_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair8_filtered)){
  if(pair8_filtered$AD2[i]<0.7){
    if(pair8_filtered[i,"GT_illumina"]=="1/1"|pair8_filtered[i,"GT_illumina"]=="0/1"){
      pair8_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair8_filtered[i,"GT_illumina"] <- pair8_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair8_filtered[i,"GT_illumina"]=="0/1"|pair8_filtered[i,"GT_illumina"]=="0/0"){
      pair8_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair8_filtered[i,"GT_illumina"] <- pair8_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair8_filtered$GT_nanopore)
if(length(idx)>0){
  pair8_filtered <- pair8_filtered[-idx,]
}
idx <- grep("[.]",pair8_filtered$GT_illumina)
if(length(idx)>0){
  pair8_filtered <- pair8_filtered[-idx,]
}
nrow(pair8_filtered) #3913

pair8_filtered$callNP <- unlist(lapply(pair8_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair8_filtered$callIL <- unlist(lapply(pair8_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair8_filtered$GT_nanopore)
unique(pair8_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  


idx <- grep("1/2",pair8_filtered$GT_nanopore)
pair8_filtered[idx,]
#pos 55553 
#AD doens't reach 0.7 in any allele 
#pos 1831219 second allele is indel, so put to 0
pair8_filtered[idx,"callNP"] <- c(0,0)

pair8_filtered$callNP <- as.numeric(pair8_filtered$callNP)
pair8_filtered$callIL <- as.numeric(pair8_filtered$callIL)

pair8upset <- upset(pair8_filtered[17:18])

table(c(pair8_filtered[c("callNP","callIL")]))
pair8_filtered[which(pair8_filtered$callNP==0&pair8_filtered$callIL==1),]


#pair9
pair9 <- read.table("pair9.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair9) <- head1
pair9[11] <- NULL

pair9$AD1_ref <- unlist(lapply(pair9$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair9$AD1_alt <- unlist(lapply(pair9$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair9$AD1 <- as.numeric(pair9$AD1_alt)/(as.numeric(pair9$AD1_ref)+as.numeric(pair9$AD1_alt))
pair9$AD2_ref <- unlist(lapply(pair9$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair9$AD2_alt <- unlist(lapply(pair9$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair9$AD2 <- as.numeric(pair9$AD2_alt)/(as.numeric(pair9$AD2_ref)+as.numeric(pair9$AD2_alt))
pair9$DP_nanopore <- as.numeric(pair9$DP_nanopore)
pair9$DP_illumina <- as.numeric(pair9$DP_illumina)

##Filter by DP 
idx <- which(pair9$DP_nanopore<10|pair9$DP_illumina<10)
if(length(idx)>0){
pair9_filtered <- pair9[-idx,]
}else{
pair9_filtered <- pair9
}
#filter by AD 
for(i in 1:nrow(pair9_filtered)){
  if(pair9_filtered$AD1[i]<0.7){
    if(pair9_filtered[i,"GT_nanopore"]=="1/1"|pair9_filtered[i,"GT_nanopore"]=="0/1"){
      pair9_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair9_filtered[i,"GT_nanopore"] <- pair9_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair9_filtered[i,"GT_nanopore"]=="0/1"|pair9_filtered[i,"GT_nanopore"]=="0/0"){
      pair9_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair9_filtered[i,"GT_nanopore"] <- pair9_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair9_filtered)){
  if(pair9_filtered$AD2[i]<0.7){
    if(pair9_filtered[i,"GT_illumina"]=="1/1"|pair9_filtered[i,"GT_illumina"]=="0/1"){
      pair9_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair9_filtered[i,"GT_illumina"] <- pair9_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair9_filtered[i,"GT_illumina"]=="0/1"|pair9_filtered[i,"GT_illumina"]=="0/0"){
      pair9_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair9_filtered[i,"GT_illumina"] <- pair9_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair9_filtered$GT_nanopore)
if(length(idx)>0){
  pair9_filtered <- pair9_filtered[-idx,]
}
idx <- grep("[.]",pair9_filtered$GT_illumina)
if(length(idx)>0){
  pair9_filtered <- pair9_filtered[-idx,]
}
nrow(pair9_filtered) #3878

pair9_filtered$callNP <- unlist(lapply(pair9_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair9_filtered$callIL <- unlist(lapply(pair9_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair9_filtered$GT_nanopore)
unique(pair9_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair9_filtered$GT_nanopore)
pair9_filtered[idx,]  ##it is a deletion 
#pos 55553 same as before, and 1831219 same, indel in nanopore, SNP in illumina 

pair9_filtered[idx,c("callNP","callIL")] <- c(0,0,1,1)


pair9_filtered$callNP <- as.numeric(pair9_filtered$callNP)
pair9_filtered$callIL <- as.numeric(pair9_filtered$callIL)

pair9upset <- upset(pair9_filtered[17:18])

table(c(pair9_filtered[c("callNP","callIL")]))
pair9_filtered[which(pair9_filtered$callNP==0&pair9_filtered$callIL==1),]
pair9_filtered[which(pair9_filtered$callNP==1&pair9_filtered$callIL==0),]


#pair10
pair10 <- read.table("pair10.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair10) <- head1
pair10[11] <- NULL

pair10$AD1_ref <- unlist(lapply(pair10$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair10$AD1_alt <- unlist(lapply(pair10$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair10$AD1 <- as.numeric(pair10$AD1_alt)/(as.numeric(pair10$AD1_ref)+as.numeric(pair10$AD1_alt))
pair10$AD2_ref <- unlist(lapply(pair10$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair10$AD2_alt <- unlist(lapply(pair10$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair10$AD2 <- as.numeric(pair10$AD2_alt)/(as.numeric(pair10$AD2_ref)+as.numeric(pair10$AD2_alt))
pair10$DP_nanopore <- as.numeric(pair10$DP_nanopore)
pair10$DP_illumina <- as.numeric(pair10$DP_illumina)

##Filter by DP 
idx <- which(pair10$DP_nanopore<10|pair10$DP_illumina<10)
if(length(idx)>0){
pair10_filtered <- pair10[-idx,]
}else{
pair10_filtered <- pair10
}
#filter by AD 
##Error in pos 730077 but is an indel, delete 
pair10_filtered <- pair10_filtered[-737,]

for(i in 1:nrow(pair10_filtered)){
  if(pair10_filtered$AD1[i]<0.7){
    if(pair10_filtered[i,"GT_nanopore"]=="1/1"|pair10_filtered[i,"GT_nanopore"]=="0/1"){
      pair10_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair10_filtered[i,"GT_nanopore"] <- pair10_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair10_filtered[i,"GT_nanopore"]=="0/1"|pair10_filtered[i,"GT_nanopore"]=="0/0"){
      pair10_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair10_filtered[i,"GT_nanopore"] <- pair10_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair10_filtered)){
  if(pair10_filtered$AD2[i]<0.7){
    if(pair10_filtered[i,"GT_illumina"]=="1/1"|pair10_filtered[i,"GT_illumina"]=="0/1"){
      pair10_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair10_filtered[i,"GT_illumina"] <- pair10_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair10_filtered[i,"GT_illumina"]=="0/1"|pair10_filtered[i,"GT_illumina"]=="0/0"){
      pair10_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair10_filtered[i,"GT_illumina"] <- pair10_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair10_filtered$GT_nanopore)
if(length(idx)>0){
  pair10_filtered <- pair10_filtered[-idx,]
}
idx <- grep("[.]",pair10_filtered$GT_illumina)
if(length(idx)>0){
  pair10_filtered <- pair10_filtered[-idx,]
}
nrow(pair10_filtered) #3874

pair10_filtered$callNP <- unlist(lapply(pair10_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair10_filtered$callIL <- unlist(lapply(pair10_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair10_filtered$GT_nanopore)
unique(pair10_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("1/2",pair10_filtered$GT_nanopore)
pair10_filtered[idx,]  ##it is a deletion 
#pos 55553 same as before, and 1831219 same, indel in nanopore, SNP in illumina 

pair10_filtered[idx,c("callNP")] <- c(0,0)

idx <- grep("1/2",pair10_filtered$GT_illumina)
pair10_filtered[idx,]
#pos 549361 is a SNP in Nanopore, and a indel in Illumina 
pair10_filtered[idx,"callIL"] <- 0


pair10_filtered$callNP <- as.numeric(pair10_filtered$callNP)
pair10_filtered$callIL <- as.numeric(pair10_filtered$callIL)

pair10upset <- upset(pair10_filtered[17:18])

table(c(pair10_filtered[c("callNP","callIL")]))
pair10_filtered[which(pair10_filtered$callNP==0&pair10_filtered$callIL==1),]
pair10_filtered[which(pair10_filtered$callNP==1&pair10_filtered$callIL==0),]

pair1_filtered[which(pair1_filtered$GT_nanopore=="2/2"|pair1_filtered$GT_nanopore=="1/2"|pair1_filtered$GT_illumina=="2/2"|pair1_filtered$GT_illumina=="1/2"),]
pair2_filtered[which(pair2_filtered$GT_nanopore=="2/2"|pair2_filtered$GT_nanopore=="1/2"|pair2_filtered$GT_illumina=="2/2"|pair2_filtered$GT_illumina=="1/2"),]
pair3_filtered[which(pair3_filtered$GT_nanopore=="2/2"|pair3_filtered$GT_nanopore=="1/2"|pair3_filtered$GT_illumina=="2/2"|pair3_filtered$GT_illumina=="1/2"),]
pair4_filtered[which(pair4_filtered$GT_nanopore=="2/2"|pair4_filtered$GT_nanopore=="1/2"|pair4_filtered$GT_illumina=="2/2"|pair4_filtered$GT_illumina=="1/2"),]
pair5_filtered[which(pair5_filtered$GT_nanopore=="2/2"|pair5_filtered$GT_nanopore=="1/2"|pair5_filtered$GT_illumina=="2/2"|pair5_filtered$GT_illumina=="1/2"),]
pair6_filtered[which(pair6_filtered$GT_nanopore=="2/2"|pair6_filtered$GT_nanopore=="1/2"|pair6_filtered$GT_illumina=="2/2"|pair6_filtered$GT_illumina=="1/2"),]




#####pairs indels extracted, called pair1.indels.txt 
pair1indels <- read.table("pair1.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair1indels) <- head1
pair1indels[11] <- NULL

pair1indels$AD1_ref <- unlist(lapply(pair1indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair1indels$AD1_alt <- unlist(lapply(pair1indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair1indels$AD1 <- as.numeric(pair1indels$AD1_alt)/(as.numeric(pair1indels$AD1_ref)+as.numeric(pair1indels$AD1_alt))
pair1indels$AD2_ref <- unlist(lapply(pair1indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair1indels$AD2_alt <- unlist(lapply(pair1indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair1indels$AD2 <- as.numeric(pair1indels$AD2_alt)/(as.numeric(pair1indels$AD2_ref)+as.numeric(pair1indels$AD2_alt))
pair1indels$DP_nanopore <- as.numeric(pair1indels$DP_nanopore)
pair1indels$DP_illumina <- as.numeric(pair1indels$DP_illumina)

##Filter by DP 
idx <- which(pair1indels$DP_nanopore<10|pair1indels$DP_illumina<10)
pair1indels_filtered <- pair1indels[-idx,]
#filter by AD 
for(i in 1:nrow(pair1indels_filtered)){
  if(pair1indels_filtered$AD1[i]<0.7){
    if(pair1indels_filtered[i,"GT_nanopore"]=="1/1"|pair1indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair1indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair1indels_filtered[i,"GT_nanopore"] <- pair1indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair1indels_filtered[i,"GT_nanopore"]=="0/1"|pair1indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair1indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair1indels_filtered[i,"GT_nanopore"] <- pair1indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair1indels_filtered)){
  if(pair1indels_filtered$AD2[i]<0.7){
    if(pair1indels_filtered[i,"GT_illumina"]=="1/1"|pair1indels_filtered[i,"GT_illumina"]=="0/1"){
      pair1indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair1indels_filtered[i,"GT_illumina"] <- pair1indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair1indels_filtered[i,"GT_illumina"]=="0/1"|pair1indels_filtered[i,"GT_illumina"]=="0/0"){
      pair1indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair1indels_filtered[i,"GT_illumina"] <- pair1indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair1indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair1indels_filtered <- pair1indels_filtered[-idx,]
}
idx <- grep("[.]",pair1indels_filtered$GT_illumina)
if(length(idx)>0){
  pair1indels_filtered <- pair1indels_filtered[-idx,]
}
nrow(pair1indels_filtered) #295

pair1indels_filtered$callNP <- unlist(lapply(pair1indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair1indels_filtered$callIL <- unlist(lapply(pair1indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair1indels_filtered$GT_nanopore)
unique(pair1indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("1/2",pair1indels_filtered$GT_nanopore)
pair1indels_filtered[idx,]
#AD doens't reach 0.7 in 854252 for any of the alleles in nanopore, so 0
#AD doesn't reach 0.7 in 1365837 for any of the alleles in nanopore, so 0
pair1indels_filtered[idx,"callNP"] <- c(0,0)

pair1indels_filtered$callNP <- as.numeric(pair1indels_filtered$callNP)
pair1indels_filtered$callIL <- as.numeric(pair1indels_filtered$callIL)

pair1indelsupset <- upset(pair1indels_filtered[17:18])

table(c(pair1indels_filtered[c("callNP","callIL")]))
pair1indels_filtered[which(pair1indels_filtered$callNP==0&pair1indels_filtered$callIL==1),]
pair1indels_filtered[which(pair1indels_filtered$callNP==1&pair1indels_filtered$callIL==0),]


#pair2indels 
pair2indels <- read.table("pair2.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair2indels) <- head1
pair2indels[11] <- NULL

pair2indels$AD1_ref <- unlist(lapply(pair2indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair2indels$AD1_alt <- unlist(lapply(pair2indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair2indels$AD1 <- as.numeric(pair2indels$AD1_alt)/(as.numeric(pair2indels$AD1_ref)+as.numeric(pair2indels$AD1_alt))
pair2indels$AD2_ref <- unlist(lapply(pair2indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair2indels$AD2_alt <- unlist(lapply(pair2indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair2indels$AD2 <- as.numeric(pair2indels$AD2_alt)/(as.numeric(pair2indels$AD2_ref)+as.numeric(pair2indels$AD2_alt))
pair2indels$DP_nanopore <- as.numeric(pair2indels$DP_nanopore)
pair2indels$DP_illumina <- as.numeric(pair2indels$DP_illumina)

##Filter by DP 
idx <- which(pair2indels$DP_nanopore<10|pair2indels$DP_illumina<10)
pair2indels_filtered <- pair2indels[-idx,]
#filter by AD 
for(i in 1:nrow(pair2indels_filtered)){
  if(pair2indels_filtered$AD1[i]<0.7){
    if(pair2indels_filtered[i,"GT_nanopore"]=="1/1"|pair2indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair2indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair2indels_filtered[i,"GT_nanopore"] <- pair2indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair2indels_filtered[i,"GT_nanopore"]=="0/1"|pair2indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair2indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair2indels_filtered[i,"GT_nanopore"] <- pair2indels_filtered[i,"GT_nanopore"]
    }
  }
}

for(i in 1:nrow(pair2indels_filtered)){
  if(pair2indels_filtered$AD2[i]<0.7){
    if(pair2indels_filtered[i,"GT_illumina"]=="1/1"|pair2indels_filtered[i,"GT_illumina"]=="0/1"){
      pair2indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair2indels_filtered[i,"GT_illumina"] <- pair2indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair2indels_filtered[i,"GT_illumina"]=="0/1"|pair2indels_filtered[i,"GT_illumina"]=="0/0"){
      pair2indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair2indels_filtered[i,"GT_illumina"] <- pair2indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair2indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair2indels_filtered <- pair2indels_filtered[-idx,]
}
idx <- grep("[.]",pair2indels_filtered$GT_illumina)
if(length(idx)>0){
  pair2indels_filtered <- pair2indels_filtered[-idx,]
}
nrow(pair2indels_filtered) #295
pair2indels_filtered$callNP <- unlist(lapply(pair2indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair2indels_filtered$callIL <- unlist(lapply(pair2indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair2indels_filtered$GT_nanopore)
unique(pair2indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("1/2",pair2indels_filtered$GT_nanopore)
pair2indels_filtered[idx,]
#854252 AD doens't reach 0.7 in any allele, and same in 1365837 
pair2indels_filtered[idx,"callNP"] <- c(0,0)

pair2indels_filtered$callNP <- as.numeric(pair2indels_filtered$callNP)
pair2indels_filtered$callIL <- as.numeric(pair2indels_filtered$callIL)

pair2indelsupset <- upset(pair2indels_filtered[17:18])

table(c(pair2indels_filtered[c("callNP","callIL")]))
pair2indels_filtered[which(pair2indels_filtered$callNP==0&pair2indels_filtered$callIL==1),]
pair2indels_filtered[which(pair2indels_filtered$callNP==1&pair2indels_filtered$callIL==0),]


#pair3 indels 
pair3indels <- read.table("pair3.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair3indels) <- head1
pair3indels[11] <- NULL

pair3indels$AD1_ref <- unlist(lapply(pair3indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair3indels$AD1_alt <- unlist(lapply(pair3indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair3indels$AD1 <- as.numeric(pair3indels$AD1_alt)/(as.numeric(pair3indels$AD1_ref)+as.numeric(pair3indels$AD1_alt))
pair3indels$AD2_ref <- unlist(lapply(pair3indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair3indels$AD2_alt <- unlist(lapply(pair3indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair3indels$AD2 <- as.numeric(pair3indels$AD2_alt)/(as.numeric(pair3indels$AD2_ref)+as.numeric(pair3indels$AD2_alt))
pair3indels$DP_nanopore <- as.numeric(pair3indels$DP_nanopore)
pair3indels$DP_illumina <- as.numeric(pair3indels$DP_illumina)

##Filter by DP 
idx <- which(pair3indels$DP_nanopore<10|pair3indels$DP_illumina<10)
pair3indels_filtered <- pair3indels[-idx,]
#filter by AD 
for(i in 1:nrow(pair3indels_filtered)){
  if(pair3indels_filtered$AD1[i]<0.7){
    if(pair3indels_filtered[i,"GT_nanopore"]=="1/1"|pair3indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair3indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair3indels_filtered[i,"GT_nanopore"] <- pair3indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair3indels_filtered[i,"GT_nanopore"]=="0/1"|pair3indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair3indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair3indels_filtered[i,"GT_nanopore"] <- pair3indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair3indels_filtered)){
  if(pair3indels_filtered$AD2[i]<0.7){
    if(pair3indels_filtered[i,"GT_illumina"]=="1/1"|pair3indels_filtered[i,"GT_illumina"]=="0/1"){
      pair3indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair3indels_filtered[i,"GT_illumina"] <- pair3indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair3indels_filtered[i,"GT_illumina"]=="0/1"|pair3indels_filtered[i,"GT_illumina"]=="0/0"){
      pair3indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair3indels_filtered[i,"GT_illumina"] <- pair3indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair3indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair3indels_filtered <- pair3indels_filtered[-idx,]
}
idx <- grep("[.]",pair3indels_filtered$GT_illumina)
if(length(idx)>0){
  pair3indels_filtered <- pair3indels_filtered[-idx,]
}
nrow(pair3indels_filtered) #293

pair3indels_filtered$callNP <- unlist(lapply(pair3indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair3indels_filtered$callIL <- unlist(lapply(pair3indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair3indels_filtered$GT_nanopore)
unique(pair3indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("1/2",pair3indels_filtered$GT_nanopore)
pair3indels_filtered[idx,]
#AD in nanopore in 854252 and 1365837 doens't reach 0.7 in any allele 
pair3indels_filtered[idx,"callNP"] <- c(0,0)

pair3indels_filtered$callNP <- as.numeric(pair3indels_filtered$callNP)
pair3indels_filtered$callIL <- as.numeric(pair3indels_filtered$callIL)

pair3indelsupset <- upset(pair3indels_filtered[17:18])

table(c(pair3indels_filtered[c("callNP","callIL")]))
pair3indels_filtered[which(pair3indels_filtered$callNP==0&pair3indels_filtered$callIL==1),]
pair3indels_filtered[which(pair3indels_filtered$callNP==1&pair3indels_filtered$callIL==0),]

#pair4 indels 
pair4indels <- read.table("pair4.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair4indels) <- head1
pair4indels[11] <- NULL

pair4indels$AD1_ref <- unlist(lapply(pair4indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair4indels$AD1_alt <- unlist(lapply(pair4indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair4indels$AD1 <- as.numeric(pair4indels$AD1_alt)/(as.numeric(pair4indels$AD1_ref)+as.numeric(pair4indels$AD1_alt))
pair4indels$AD2_ref <- unlist(lapply(pair4indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair4indels$AD2_alt <- unlist(lapply(pair4indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair4indels$AD2 <- as.numeric(pair4indels$AD2_alt)/(as.numeric(pair4indels$AD2_ref)+as.numeric(pair4indels$AD2_alt))
pair4indels$DP_nanopore <- as.numeric(pair4indels$DP_nanopore)
pair4indels$DP_illumina <- as.numeric(pair4indels$DP_illumina)

##Filter by DP 
idx <- which(pair4indels$DP_nanopore<10|pair4indels$DP_illumina<10)
pair4indels_filtered <- pair4indels[-idx,]
#filter by AD 
for(i in 1:nrow(pair4indels_filtered)){
  if(pair4indels_filtered$AD1[i]<0.7){
    if(pair4indels_filtered[i,"GT_nanopore"]=="1/1"|pair4indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair4indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair4indels_filtered[i,"GT_nanopore"] <- pair4indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair4indels_filtered[i,"GT_nanopore"]=="0/1"|pair4indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair4indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair4indels_filtered[i,"GT_nanopore"] <- pair4indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair4indels_filtered)){
  if(pair4indels_filtered$AD2[i]<0.7){
    if(pair4indels_filtered[i,"GT_illumina"]=="1/1"|pair4indels_filtered[i,"GT_illumina"]=="0/1"){
      pair4indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair4indels_filtered[i,"GT_illumina"] <- pair4indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair4indels_filtered[i,"GT_illumina"]=="0/1"|pair4indels_filtered[i,"GT_illumina"]=="0/0"){
      pair4indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair4indels_filtered[i,"GT_illumina"] <- pair4indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair4indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair4indels_filtered <- pair4indels_filtered[-idx,]
}
idx <- grep("[.]",pair4indels_filtered$GT_illumina)
if(length(idx)>0){
  pair4indels_filtered <- pair4indels_filtered[-idx,]
}
nrow(pair4indels_filtered) #292

pair4indels_filtered$callNP <- unlist(lapply(pair4indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair4indels_filtered$callIL <- unlist(lapply(pair4indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair4indels_filtered$GT_nanopore)
unique(pair4indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1 
idx <- grep("1/2",pair4indels_filtered$GT_nanopore)
pair4indels_filtered[idx,]
#AD in 854252 and 1365837 doens't reach 0.7 in any allele in nanopore 
pair4indels_filtered[idx,"callNP"] <- c(0,0)

pair4indels_filtered$callNP <- as.numeric(pair4indels_filtered$callNP)
pair4indels_filtered$callIL <- as.numeric(pair4indels_filtered$callIL)

pair4indelsupset <- upset(pair4indels_filtered[17:18])

table(c(pair4indels_filtered[c("callNP","callIL")]))
pair4indels_filtered[which(pair4indels_filtered$callNP==0&pair4indels_filtered$callIL==1),]
pair4indels_filtered[which(pair4indels_filtered$callNP==1&pair4indels_filtered$callIL==0),]


#pair5 indels 
pair5indels <- read.table("pair5.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair5indels) <- head1
pair5indels[11] <- NULL

pair5indels$AD1_ref <- unlist(lapply(pair5indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair5indels$AD1_alt <- unlist(lapply(pair5indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair5indels$AD1 <- as.numeric(pair5indels$AD1_alt)/(as.numeric(pair5indels$AD1_ref)+as.numeric(pair5indels$AD1_alt))
pair5indels$AD2_ref <- unlist(lapply(pair5indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair5indels$AD2_alt <- unlist(lapply(pair5indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair5indels$AD2 <- as.numeric(pair5indels$AD2_alt)/(as.numeric(pair5indels$AD2_ref)+as.numeric(pair5indels$AD2_alt))
pair5indels$DP_nanopore <- as.numeric(pair5indels$DP_nanopore)
pair5indels$DP_illumina <- as.numeric(pair5indels$DP_illumina)

##Filter by DP 
idx <- which(pair5indels$DP_nanopore<10|pair5indels$DP_illumina<10)
if(length(idx)>0){
pair5indels_filtered <- pair5indels[-idx,]
}else{
pair5indels_filtered <- pair5indels
}
#filter by AD 
for(i in 1:nrow(pair5indels_filtered)){
  if(pair5indels_filtered$AD1[i]<0.7){
    if(pair5indels_filtered[i,"GT_nanopore"]=="1/1"|pair5indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair5indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair5indels_filtered[i,"GT_nanopore"] <- pair5indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair5indels_filtered[i,"GT_nanopore"]=="0/1"|pair5indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair5indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair5indels_filtered[i,"GT_nanopore"] <- pair5indels_filtered[i,"GT_nanopore"]
    }
  }
}
pair5indels_filtered[89,"AD2"] <- 1
for(i in 1:nrow(pair5indels_filtered)){
  if(pair5indels_filtered$AD2[i]<0.7){
    if(pair5indels_filtered[i,"GT_illumina"]=="1/1"|pair5indels_filtered[i,"GT_illumina"]=="0/1"){
      pair5indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair5indels_filtered[i,"GT_illumina"] <- pair5indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair5indels_filtered[i,"GT_illumina"]=="0/1"|pair5indels_filtered[i,"GT_illumina"]=="0/0"){
      pair5indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair5indels_filtered[i,"GT_illumina"] <- pair5indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair5indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair5indels_filtered <- pair5indels_filtered[-idx,]
}
idx <- grep("[.]",pair5indels_filtered$GT_illumina)
if(length(idx)>0){
  pair5indels_filtered <- pair5indels_filtered[-idx,]
}
nrow(pair5indels_filtered) #302

pair5indels_filtered$callNP <- unlist(lapply(pair5indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair5indels_filtered$callIL <- unlist(lapply(pair5indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))
unique(pair5indels_filtered$GT_nanopore)
unique(pair5indels_filtered$GT_illumina)


##Fix 2/2 and 1/2 2/1 
idx <- grep("2/2",pair5indels_filtered$GT_nanopore)
pair5indels_filtered[idx,]
#1365837 both are 2/2, same deletion, AD of 0.82 and 1 so put 1 1
pair5indels_filtered[idx,c("callNP","callIL")] <- c(1,1)
idx <- grep("1/2",pair5indels_filtered$GT_nanopore)
pair5indels_filtered[idx,]
#AD in 854252 doens't reach 0.7 in any allele in nanopore 
pair5indels_filtered[idx,"callNP"] <- 0


pair5indels_filtered$callNP <- as.numeric(pair5indels_filtered$callNP)
pair5indels_filtered$callIL <- as.numeric(pair5indels_filtered$callIL)

pair5indelsupset <- upset(pair5indels_filtered[17:18])

table(c(pair5indels_filtered[c("callNP","callIL")]))
pair5indels_filtered[which(pair5indels_filtered$callNP==0&pair5indels_filtered$callIL==1),]
pair5indels_filtered[which(pair5indels_filtered$callNP==1&pair5indels_filtered$callIL==0),]

#pair6 indels 
pair6indels <- read.table("pair6.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair6indels) <- head1
pair6indels[11] <- NULL

pair6indels$AD1_ref <- unlist(lapply(pair6indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair6indels$AD1_alt <- unlist(lapply(pair6indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair6indels$AD1 <- as.numeric(pair6indels$AD1_alt)/(as.numeric(pair6indels$AD1_ref)+as.numeric(pair6indels$AD1_alt))
pair6indels$AD2_ref <- unlist(lapply(pair6indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair6indels$AD2_alt <- unlist(lapply(pair6indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair6indels$AD2 <- as.numeric(pair6indels$AD2_alt)/(as.numeric(pair6indels$AD2_ref)+as.numeric(pair6indels$AD2_alt))
pair6indels$DP_nanopore <- as.numeric(pair6indels$DP_nanopore)
pair6indels$DP_illumina <- as.numeric(pair6indels$DP_illumina)


##Filter by DP 
idx <- which(pair6indels$DP_nanopore<10|pair6indels$DP_illumina<10)
if(length(idx)>0){
pair6indels_filtered <- pair6indels[-idx,]
}else{
pair6indels_filtered <- pair6indels
}
#filter by AD 
for(i in 1:nrow(pair6indels_filtered)){
  if(pair6indels_filtered$AD1[i]<0.7){
    if(pair6indels_filtered[i,"GT_nanopore"]=="1/1"|pair6indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair6indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair6indels_filtered[i,"GT_nanopore"] <- pair6indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair6indels_filtered[i,"GT_nanopore"]=="0/1"|pair6indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair6indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair6indels_filtered[i,"GT_nanopore"] <- pair6indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair6indels_filtered)){
  if(pair6indels_filtered$AD2[i]<0.7){
    if(pair6indels_filtered[i,"GT_illumina"]=="1/1"|pair6indels_filtered[i,"GT_illumina"]=="0/1"){
      pair6indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair6indels_filtered[i,"GT_illumina"] <- pair6indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair6indels_filtered[i,"GT_illumina"]=="0/1"|pair6indels_filtered[i,"GT_illumina"]=="0/0"){
      pair6indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair6indels_filtered[i,"GT_illumina"] <- pair6indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair6indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair6indels_filtered <- pair6indels_filtered[-idx,]
}
idx <- grep("[.]",pair6indels_filtered$GT_illumina)
if(length(idx)>0){
  pair6indels_filtered <- pair6indels_filtered[-idx,]
}
nrow(pair6indels_filtered) #292

pair6indels_filtered$callNP <- unlist(lapply(pair6indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair6indels_filtered$callIL <- unlist(lapply(pair6indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair6indels_filtered$GT_nanopore)
unique(pair6indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair6indels_filtered$GT_nanopore)
pair6indels_filtered[idx,]  ##it is a deletion 
#3015119 both 2/2, change for 1 1 
pair6indels_filtered[idx,c("callNP","callIL")] <- c(1,1)
idx <- grep("1/2",pair6indels_filtered$GT_nanopore)
pair6indels_filtered[idx,]  ##it is a deletion 
#854252 0 in nanopore, 1 in illumina 
#1365837 0 in nanopore, 1 in illumina 
pair6indels_filtered[idx,"callNP"] <- c(0,0)

pair6indels_filtered$callNP <- as.numeric(pair6indels_filtered$callNP)
pair6indels_filtered$callIL <- as.numeric(pair6indels_filtered$callIL)

pair6indelsupset <- upset(pair6indels_filtered[17:18])

table(c(pair6indels_filtered[c("callNP","callIL")]))

pair6indels_filtered[which(pair6indels_filtered$callNP==0&pair6indels_filtered$callIL==1),]
pair6indels_filtered[which(pair6indels_filtered$callNP==1&pair6indels_filtered$callIL==0),]


#pair7 indels 
pair7indels <- read.table("pair7.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair7indels) <- head1
pair7indels[11] <- NULL

pair7indels$AD1_ref <- unlist(lapply(pair7indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair7indels$AD1_alt <- unlist(lapply(pair7indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair7indels$AD1 <- as.numeric(pair7indels$AD1_alt)/(as.numeric(pair7indels$AD1_ref)+as.numeric(pair7indels$AD1_alt))
pair7indels$AD2_ref <- unlist(lapply(pair7indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair7indels$AD2_alt <- unlist(lapply(pair7indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair7indels$AD2 <- as.numeric(pair7indels$AD2_alt)/(as.numeric(pair7indels$AD2_ref)+as.numeric(pair7indels$AD2_alt))
pair7indels$DP_nanopore <- as.numeric(pair7indels$DP_nanopore)
pair7indels$DP_illumina <- as.numeric(pair7indels$DP_illumina)

##Filter by DP 
idx <- which(pair7indels$DP_nanopore<10|pair7indels$DP_illumina<10)
if(length(idx)>0){
pair7indels_filtered <- pair7indels[-idx,]
}else{
pair7indels_filtered <- pair7indels
}
#filter by AD 
for(i in 1:nrow(pair7indels_filtered)){
  if(pair7indels_filtered$AD1[i]<0.7){
    if(pair7indels_filtered[i,"GT_nanopore"]=="1/1"|pair7indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair7indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair7indels_filtered[i,"GT_nanopore"] <- pair7indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair7indels_filtered[i,"GT_nanopore"]=="0/1"|pair7indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair7indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair7indels_filtered[i,"GT_nanopore"] <- pair7indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair7indels_filtered)){
  if(pair7indels_filtered$AD2[i]<0.7){
    if(pair7indels_filtered[i,"GT_illumina"]=="1/1"|pair7indels_filtered[i,"GT_illumina"]=="0/1"){
      pair7indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair7indels_filtered[i,"GT_illumina"] <- pair7indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair7indels_filtered[i,"GT_illumina"]=="0/1"|pair7indels_filtered[i,"GT_illumina"]=="0/0"){
      pair7indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair7indels_filtered[i,"GT_illumina"] <- pair7indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair7indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair7indels_filtered <- pair7indels_filtered[-idx,]
}
idx <- grep("[.]",pair7indels_filtered$GT_illumina)
if(length(idx)>0){
  pair7indels_filtered <- pair7indels_filtered[-idx,]
}
nrow(pair7indels_filtered) #284

pair7indels_filtered$callNP <- unlist(lapply(pair7indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair7indels_filtered$callIL <- unlist(lapply(pair7indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair7indels_filtered$GT_nanopore)
unique(pair7indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("1/2",pair7indels_filtered$GT_nanopore)
pair7indels_filtered[idx,]
#pos 854252 should be 0 in nanopore, 1365837, same 
#AD doens't reach 0.7 in any allele 
pair7indels_filtered[idx,"callNP"] <- c(0,0)

pair7indels_filtered$callNP <- as.numeric(pair7indels_filtered$callNP)
pair7indels_filtered$callIL <- as.numeric(pair7indels_filtered$callIL)

pair7indelsupset <- upset(pair7indels_filtered[17:18])

table(c(pair7indels_filtered[c("callNP","callIL")]))
pair7indels_filtered[which(pair7indels_filtered$callNP==0&pair7indels_filtered$callIL==1),]
pair7indels_filtered[which(pair7indels_filtered$callNP==1&pair7indels_filtered$callIL==0),]

#pair8 indels 
pair8indels <- read.table("pair8.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair8indels) <- head1
pair8indels[11] <- NULL

pair8indels$AD1_ref <- unlist(lapply(pair8indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair8indels$AD1_alt <- unlist(lapply(pair8indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair8indels$AD1 <- as.numeric(pair8indels$AD1_alt)/(as.numeric(pair8indels$AD1_ref)+as.numeric(pair8indels$AD1_alt))
pair8indels$AD2_ref <- unlist(lapply(pair8indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair8indels$AD2_alt <- unlist(lapply(pair8indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair8indels$AD2 <- as.numeric(pair8indels$AD2_alt)/(as.numeric(pair8indels$AD2_ref)+as.numeric(pair8indels$AD2_alt))
pair8indels$DP_nanopore <- as.numeric(pair8indels$DP_nanopore)
pair8indels$DP_illumina <- as.numeric(pair8indels$DP_illumina)

##Filter by DP 
idx <- which(pair8indels$DP_nanopore<10|pair8indels$DP_illumina<10)
if(length(idx)>0){
pair8indels_filtered <- pair8indels[-idx,]
}else{
pair8indels_filtered <- pair8indels
}
#filter by AD 
for(i in 1:nrow(pair8indels_filtered)){
  if(pair8indels_filtered$AD1[i]<0.7){
    if(pair8indels_filtered[i,"GT_nanopore"]=="1/1"|pair8indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair8indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair8indels_filtered[i,"GT_nanopore"] <- pair8indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair8indels_filtered[i,"GT_nanopore"]=="0/1"|pair8indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair8indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair8indels_filtered[i,"GT_nanopore"] <- pair8indels_filtered[i,"GT_nanopore"]
    }
  }
}
pair8indels_filtered[159,"AD2"] <- 1
for(i in 1:nrow(pair8indels_filtered)){
  if(pair8indels_filtered$AD2[i]<0.7){
    if(pair8indels_filtered[i,"GT_illumina"]=="1/1"|pair8indels_filtered[i,"GT_illumina"]=="0/1"){
      pair8indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair8indels_filtered[i,"GT_illumina"] <- pair8indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair8indels_filtered[i,"GT_illumina"]=="0/1"|pair8indels_filtered[i,"GT_illumina"]=="0/0"){
      pair8indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair8indels_filtered[i,"GT_illumina"] <- pair8indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair8indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair8indels_filtered <- pair8indels_filtered[-idx,]
}
idx <- grep("[.]",pair8indels_filtered$GT_illumina)
if(length(idx)>0){
  pair8indels_filtered <- pair8indels_filtered[-idx,]
}
nrow(pair8indels_filtered) #280

pair8indels_filtered$callNP <- unlist(lapply(pair8indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair8indels_filtered$callIL <- unlist(lapply(pair8indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair8indels_filtered$GT_nanopore)
unique(pair8indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair8indels_filtered$GT_nanopore)
pair8indels_filtered[idx,]
#2536625 both are 2/2, put 1 1
pair8indels_filtered[idx,c("callNP","callIL")] <- c(1,1)

idx <- grep("1/2",pair8indels_filtered$GT_nanopore)
pair8indels_filtered[idx,]
#pos 854252 put 0 in nanopore, and 1365837 same  
#AD doens't reach 0.7 in any allele 
pair8indels_filtered[idx,"callNP"] <- c(0,0)

pair8indels_filtered$callNP <- as.numeric(pair8indels_filtered$callNP)
pair8indels_filtered$callIL <- as.numeric(pair8indels_filtered$callIL)

pair8indelsupset <- upset(pair8indels_filtered[17:18])

table(c(pair8indels_filtered[c("callNP","callIL")]))
pair8indels_filtered[which(pair8indels_filtered$callNP==0&pair8indels_filtered$callIL==1),]
pair8indels_filtered[which(pair8indels_filtered$callNP==1&pair8indels_filtered$callIL==0),]


#pair9 indels 
pair9indels <- read.table("pair9.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair9indels) <- head1
pair9indels[11] <- NULL

pair9indels$AD1_ref <- unlist(lapply(pair9indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair9indels$AD1_alt <- unlist(lapply(pair9indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair9indels$AD1 <- as.numeric(pair9indels$AD1_alt)/(as.numeric(pair9indels$AD1_ref)+as.numeric(pair9indels$AD1_alt))
pair9indels$AD2_ref <- unlist(lapply(pair9indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair9indels$AD2_alt <- unlist(lapply(pair9indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair9indels$AD2 <- as.numeric(pair9indels$AD2_alt)/(as.numeric(pair9indels$AD2_ref)+as.numeric(pair9indels$AD2_alt))
pair9indels$DP_nanopore <- as.numeric(pair9indels$DP_nanopore)
pair9indels$DP_illumina <- as.numeric(pair9indels$DP_illumina)

##Filter by DP 
idx <- which(pair9indels$DP_nanopore<10|pair9indels$DP_illumina<10)
if(length(idx)>0){
pair9indels_filtered <- pair9indels[-idx,]
}else{
pair9indels_filtered <- pair9indels
}
#filter by AD 
pair9indels_filtered[160,"AD1"] <- 1
pair9indels_filtered[160,"AD2"] <- 1
for(i in 1:nrow(pair9indels_filtered)){
  if(pair9indels_filtered$AD1[i]<0.7){
    if(pair9indels_filtered[i,"GT_nanopore"]=="1/1"|pair9indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair9indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair9indels_filtered[i,"GT_nanopore"] <- pair9indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair9indels_filtered[i,"GT_nanopore"]=="0/1"|pair9indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair9indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair9indels_filtered[i,"GT_nanopore"] <- pair9indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair9indels_filtered)){
  if(pair9indels_filtered$AD2[i]<0.7){
    if(pair9indels_filtered[i,"GT_illumina"]=="1/1"|pair9indels_filtered[i,"GT_illumina"]=="0/1"){
      pair9indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair9indels_filtered[i,"GT_illumina"] <- pair9indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair9indels_filtered[i,"GT_illumina"]=="0/1"|pair9indels_filtered[i,"GT_illumina"]=="0/0"){
      pair9indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair9indels_filtered[i,"GT_illumina"] <- pair9indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair9indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair9indels_filtered <- pair9indels_filtered[-idx,]
}
idx <- grep("[.]",pair9indels_filtered$GT_illumina)
if(length(idx)>0){
  pair9indels_filtered <- pair9indels_filtered[-idx,]
}
nrow(pair9indels_filtered) #281

pair9indels_filtered$callNP <- unlist(lapply(pair9indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair9indels_filtered$callIL <- unlist(lapply(pair9indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair9indels_filtered$GT_nanopore)
unique(pair9indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair9indels_filtered$GT_nanopore)
pair9indels_filtered[idx,]  ##it is a deletion 
#pos 2536625 call is 1 in NP and 1 in illumina 
pair9indels_filtered[idx,c("callNP","callIL")] <- c(1,1)
idx <- grep("1/2",pair9indels_filtered$GT_nanopore)
pair9indels_filtered[idx,]
##854252 0 in NP, 1 in IL 
#1365837 same 
pair9indels_filtered[idx,"callNP"] <- c(0,0)


pair9indels_filtered$callNP <- as.numeric(pair9indels_filtered$callNP)
pair9indels_filtered$callIL <- as.numeric(pair9indels_filtered$callIL)

pair9indelsupset <- upset(pair9indels_filtered[17:18])

table(c(pair9indels_filtered[c("callNP","callIL")]))
pair9indels_filtered[which(pair9indels_filtered$callNP==0&pair9indels_filtered$callIL==1),]
pair9indels_filtered[which(pair9indels_filtered$callNP==1&pair9indels_filtered$callIL==0),]


#pair10 indels
pair10indels <- read.table("pair10.indels.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(pair10indels) <- head1
pair10indels[11] <- NULL

pair10indels$AD1_ref <- unlist(lapply(pair10indels$AD_nanopore,function(x)strsplit(x,",")[[1]][1]))
pair10indels$AD1_alt <- unlist(lapply(pair10indels$AD_nanopore,function(x)strsplit(x,",")[[1]][2]))
pair10indels$AD1 <- as.numeric(pair10indels$AD1_alt)/(as.numeric(pair10indels$AD1_ref)+as.numeric(pair10indels$AD1_alt))
pair10indels$AD2_ref <- unlist(lapply(pair10indels$AD_illumina,function(x)strsplit(x,",")[[1]][1]))
pair10indels$AD2_alt <- unlist(lapply(pair10indels$AD_illumina,function(x)strsplit(x,",")[[1]][2]))
pair10indels$AD2 <- as.numeric(pair10indels$AD2_alt)/(as.numeric(pair10indels$AD2_ref)+as.numeric(pair10indels$AD2_alt))
pair10indels$DP_nanopore <- as.numeric(pair10indels$DP_nanopore)
pair10indels$DP_illumina <- as.numeric(pair10indels$DP_illumina)

##Filter by DP 
idx <- which(pair10indels$DP_nanopore<10|pair10indels$DP_illumina<10)
if(length(idx)>0){
pair10indels_filtered <- pair10indels[-idx,]
}else{
pair10indels_filtered <- pair10indels
}
#filter by AD 
pair10indels_filtered[53,"AD1"] <- 1
pair10indels_filtered[53,"AD2"] <- 1
pair10indels_filtered[160,"AD1"] <- 1

for(i in 1:nrow(pair10indels_filtered)){
  if(pair10indels_filtered$AD1[i]<0.7){
    if(pair10indels_filtered[i,"GT_nanopore"]=="1/1"|pair10indels_filtered[i,"GT_nanopore"]=="0/1"){
      pair10indels_filtered[i,"GT_nanopore"] <- "0/0"
    }else{
      pair10indels_filtered[i,"GT_nanopore"] <- pair10indels_filtered[i,"GT_nanopore"]
    }
  }else{
    if(pair10indels_filtered[i,"GT_nanopore"]=="0/1"|pair10indels_filtered[i,"GT_nanopore"]=="0/0"){
      pair10indels_filtered[i,"GT_nanopore"] <- "1/1"
    }else{
      pair10indels_filtered[i,"GT_nanopore"] <- pair10indels_filtered[i,"GT_nanopore"]
    }
  }
}
for(i in 1:nrow(pair10indels_filtered)){
  if(pair10indels_filtered$AD2[i]<0.7){
    if(pair10indels_filtered[i,"GT_illumina"]=="1/1"|pair10indels_filtered[i,"GT_illumina"]=="0/1"){
      pair10indels_filtered[i,"GT_illumina"] <- "0/0"
    }else{
      pair10indels_filtered[i,"GT_illumina"] <- pair10indels_filtered[i,"GT_illumina"]
    }
  }else{
    if(pair10indels_filtered[i,"GT_illumina"]=="0/1"|pair10indels_filtered[i,"GT_illumina"]=="0/0"){
      pair10indels_filtered[i,"GT_illumina"] <- "1/1"
    }else{
      pair10indels_filtered[i,"GT_illumina"] <- pair10indels_filtered[i,"GT_illumina"]
    }
  }
}
idx <- grep("[.]",pair10indels_filtered$GT_nanopore)
if(length(idx)>0){
  pair10indels_filtered <- pair10indels_filtered[-idx,]
}
idx <- grep("[.]",pair10indels_filtered$GT_illumina)
if(length(idx)>0){
  pair10indels_filtered <- pair10indels_filtered[-idx,]
}
nrow(pair10indels_filtered) #281

pair10indels_filtered$callNP <- unlist(lapply(pair10indels_filtered$GT_nanopore,function(x) strsplit(x,"/")[[1]][1]))
pair10indels_filtered$callIL <- unlist(lapply(pair10indels_filtered$GT_illumina,function(x) strsplit(x,"/")[[1]][1]))

unique(pair10indels_filtered$GT_nanopore)
unique(pair10indels_filtered$GT_illumina)

##Fix 2/2 and 1/2 2/1  
idx <- grep("2/2",pair10indels_filtered$GT_nanopore)
pair10indels_filtered[idx,]  ##
#pos 854252 and 2536625, are 2/2 in both, put 1, 1
pair10indels_filtered[idx,c("callNP","callIL")] <- c(1,1)
idx <- grep("2/2",pair10indels_filtered$GT_illumina)
pair10indels_filtered[idx,]  ##
#2338194, is 0/0 in NP and 2/2 in IL
#keep 1 in IL and 0 in NP 
pair10indels_filtered[146,c("callNP","callIL")] <- c(0,1)

idx <- grep("1/2",pair10indels_filtered$GT_nanopore)
pair10indels_filtered[idx,]
#pos 1365837 is 0 in NP and 1 in IL 
pair10indels_filtered[idx,"callNP"] <- 0


pair10indels_filtered$callNP <- as.numeric(pair10indels_filtered$callNP)
pair10indels_filtered$callIL <- as.numeric(pair10indels_filtered$callIL)

pair10indelsupset <- upset(pair10indels_filtered[17:18])

table(c(pair10indels_filtered[c("callNP","callIL")]))
pair10indels_filtered[which(pair10indels_filtered$callNP==0&pair10indels_filtered$callIL==1),]
pair10indels_filtered[which(pair10indels_filtered$callNP==1&pair10indels_filtered$callIL==0),]

pair1_filtered[which(pair1_filtered$GT_nanopore=="2/2"|pair1_filtered$GT_nanopore=="1/2"|pair1_filtered$GT_illumina=="2/2"|pair1_filtered$GT_illumina=="1/2"),]
pair2_filtered[which(pair2_filtered$GT_nanopore=="2/2"|pair2_filtered$GT_nanopore=="1/2"|pair2_filtered$GT_illumina=="2/2"|pair2_filtered$GT_illumina=="1/2"),]
pair3_filtered[which(pair3_filtered$GT_nanopore=="2/2"|pair3_filtered$GT_nanopore=="1/2"|pair3_filtered$GT_illumina=="2/2"|pair3_filtered$GT_illumina=="1/2"),]
pair4_filtered[which(pair4_filtered$GT_nanopore=="2/2"|pair4_filtered$GT_nanopore=="1/2"|pair4_filtered$GT_illumina=="2/2"|pair4_filtered$GT_illumina=="1/2"),]
pair5_filtered[which(pair5_filtered$GT_nanopore=="2/2"|pair5_filtered$GT_nanopore=="1/2"|pair5_filtered$GT_illumina=="2/2"|pair5_filtered$GT_illumina=="1/2"),]
pair6_filtered[which(pair6_filtered$GT_nanopore=="2/2"|pair6_filtered$GT_nanopore=="1/2"|pair6_filtered$GT_illumina=="2/2"|pair6_filtered$GT_illumina=="1/2"),]

     
