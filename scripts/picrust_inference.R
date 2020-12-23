library(reshape2)
library(ggpubr)
library(VennDiagram)

#china
pi=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/metagenome_predictions.txt",header=TRUE,sep="\t",row.names=1)
pi=pi[,sort(colnames(pi))]
pi_p1=t(t(pi)/colSums(pi))[,1:40]
pi_p=pi_p1[rowSums(pi_p1)!=0,]
adonis(t(pi_p)~factor(c(rep(1,20),rep(2,20))))
"Terms added sequentially (first to last)

Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
factor(c(rep(1, 20), rep(2, 20)))  1   0.04823 0.048234  3.8877 0.09281  0.008 **
Residuals                         38   0.47147 0.012407         0.90719          
Total                             39   0.51970                  1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
wg=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/kegg_complete_perc.txt",sep="\t",stringsAsFactors=F,header=T,row.names=1)
wg=wg[,sort(colnames(wg))]
adonis(t(wg)~factor(c(rep(1,20),rep(2,20))))
"
                                  Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
factor(c(rep(1, 20), rep(2, 20)))  1   0.03899 0.038988  4.2326 0.10022  0.001 ***
Residuals                         38   0.35003 0.009211         0.89978           
Total                             39   0.38902                  1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"


otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/otu_table_42000.txt",header=T,row.names=1)
otu=otu[,sort(colnames(otu))]
otu_p1=t(t(otu)/colSums(otu))[,1:40]
otu_p=otu_p1[rowSums(otu_p1)!=0,]
adonis(t(otu_p)~factor(c(rep(1,20),rep(2,20))))
"                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
factor(c(rep(1, 20), rep(2, 20)))  1    0.7666 0.76661  2.9314 0.07162  0.004 **
Residuals                         38    9.9378 0.26152         0.92838          
Total                             39   10.7044                 1.00000 "

#orthology
dim(wg)[1] #13880
dim(pi_p)[1] #5612
comm2=Reduce(union, list(rownames(wg),rownames(pi_p)))
comm3=Reduce(intersect, list(rownames(wg),rownames(pi_p)))
uni=setdiff(rownames(wg),rownames(pi_p))
length(comm2)#13918
length(comm3)#5574

Kpv=list()
Ppv=list()
abdK1=list()
abdP1=list()
abdK1_1=list()
abdP1_1=list()
i=1
for (m in comm2){
  if (m %in% rownames(wg)){
    if (t.test(wg[m,1:20],wg[m,21:40])$statistic>0){
      Kpv2=-log10(t.test(wg[m,1:20],wg[m,21:40])$p.value)
    }else{
      Kpv2=log10(t.test(wg[m,1:20],wg[m,21:40])$p.value)
    }
  }else{
    Kpv2=NA
  }
  if (m %in% rownames(pi_p)){
    if (t.test(pi_p[m,1:20],pi_p[m,21:40])$statistic>0){
      Ppv2=-log10(t.test(pi_p[m,1:20],pi_p[m,21:40])$p.value)
    }else{
      Ppv2=log10(t.test(pi_p[m,1:20],pi_p[m,21:40])$p.value)
    }
  }else{
    Ppv2=NA
  }
  Kpv[i]=Kpv2
  Ppv[i]=Ppv2
  if (m %in% rownames(wg) & m %in% rownames(pi_p)){
    abdK1[i]=as.numeric(t.test(wg[m,1:20],wg[m,21:40])$estimate[1]+t.test(wg[m,1:20],wg[m,21:40])$estimate[2])
    abdP1[i]=as.numeric(t.test(pi_p[m,1:20],pi_p[m,21:40])$estimate[1]+t.test(pi_p[m,1:20],pi_p[m,21:40])$estimate[2])
    abdK1_1[i]=as.numeric(t.test(wg[m,1:20],wg[m,21:40])$estimate[1]/t.test(wg[m,1:20],wg[m,21:40])$estimate[2])
    abdP1_1[i]=as.numeric(t.test(pi_p[m,1:20],pi_p[m,21:40])$estimate[1]/t.test(pi_p[m,1:20],pi_p[m,21:40])$estimate[2])
  }
  i=i+1
}

china_K=cbind(unlist(Kpv),unlist(Ppv))
rownames(china_K)=comm2
colnames(china_K)=c("metagenome","PICRUSt")
write.csv(china_K,file="/Users/Shansun/Google\ Drive/picrust_paper/china_K.csv")


dim(wg)[1] #13880
dim(pi_p)[1] #5612
comm3=Reduce(intersect, list(rownames(wg),rownames(pi_p)))
length(comm3)#5574

pdf("/Users/shansun/Google\ Drive/picrust_paper/china_K_venn.pdf", width = 8, height = 8)
grid.newpage()
venn.plot=draw.pairwise.venn(13880, 5612, 5574, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cex=3,cat.cex=3)
grid.draw(venn.plot)
dev.off()

china_K1=china_K[china_K[,1]!=0&china_K[,2]!=0,]
china_K2=china_K[china_K[,1]==0|china_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/china_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(china_K1[,1],china_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Human 1",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

cor.test(china_K1[,1],china_K1[,2],method="spearman")
"	Pearson's product-moment correlation

data:  china_K1[, 1] and china_K1[, 2]
t = 53.811, df = 5572, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
0.5672350 0.6017941
sample estimates:
cor 
0.5847798 "

"	Spearman's rank correlation rho

data:  china_K1[, 1] and china_K1[, 2]
S = 1.4707e+10, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.4904512 "


adjP3=cbind(unlist(Kpv),unlist(Ppv))
dim(adjP3)[1]#13918
dim(adjP3[adjP3[,1]>0&adjP3[,2]>0,])[1]#1814
dim(adjP3[adjP3[,1]<0&adjP3[,2]<0,])[1]#1886
dim(adjP3[adjP3[,1]>0&adjP3[,2]<0,])[1]#673
dim(adjP3[adjP3[,1]<0&adjP3[,2]>0,])[1]#1201
dim(adjP3[adjP3[,1]==0,])[1]#38
dim(adjP3[adjP3[,2]==0,])[1]#8306
1814/(1814+673)#taxa more abundant in rural, 0.7293928 right
1886/(1886+1201)#taxa more abundant in urban, 0.6109491 right

1201/(1201+1814)#rural 0.398
673/(1886+673)#urban 0.263
(1814+1886)/13918 #0.266
(1814+1886)/(13918-38-8306)#0.664
(8306+38)/13918 #0.597


cor.test(as.numeric(unlist(Kpv)),as.numeric(unlist(Ppv)))
#right ones, both are significant and the same direction #509
correct=comm2[c(c(1:13918)[adjP3[,1]>1.30103&adjP3[,2]>1.30103],c(1:13918)[adjP3[,1]<(-1.30103)&adjP3[,2]<(-1.30103)])]
#both are insignificant #10994
insig=comm2[c(1:13918)[adjP3[,1]<1.30103&adjP3[,1]>(-1.30103)&adjP3[,2]>(-1.30103)&adjP3[,2]<1.30103]]
#both are significant but wrong direction or one is significant one is not
13918-10994-509 #2415
incorr=comm2[setdiff(setdiff(c(1:13918),insig),correct)]

#the ranks of functions are quite consistent between samples

pi_p_perm=randomizeMatrix(pi_p,null.model = "richness",iterations = 100)
wg_perm=randomizeMatrix(wg,null.model = "richness",iterations = 100)
wg_rank=apply(wg,2,rank)
cor.test(wg[,1],wg_perm[,1],method="spearman")

comm3=intersect(rownames(pi_p),rownames(wg))
cor_china=vector()
i=1
for (m in 1:dim(wg)[2]){
  cor1=cor.test(pi_p[comm3,m],wg[comm3,m],method="spearman")
  cor_china[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(wg,null.model = "richness",iterations = 100)
  for (m in 1:dim(wg)[2]){
    cor1=cor.test(pi_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_china[i]=cor1$estimate
    i=i+1
  }
}
boxplot(cor_china[41:4040])
points(rep(1,40),cor_china[1:40],col="red")
t.test(cor_china[1:40],cor_china[41:4040])

#randomize the metagenome and test our method
comm2=union(rownames(pi_p),rownames(wg))
Kpv_perm=matrix(nrow=length(comm2),ncol=100)
for (n in 1:100){
  wg_perm=randomizeMatrix(wg,null.model = "richness",iterations = 100)
  i=1
  for (m in comm2){
    if (m %in% rownames(wg)){
      if (t.test(wg_perm[m,1:20],wg_perm[m,21:40])$statistic>0){
        Kpv2=-log10(t.test(wg_perm[m,1:20],wg_perm[m,21:40])$p.value)
      }else{
        Kpv2=log10(t.test(wg_perm[m,1:20],wg_perm[m,21:40])$p.value)
      }
    }else{
      Kpv2=NA
    }
    Kpv_perm[i,n]=Kpv2
    i=i+1
  }
  print(n)
}
write.csv(Kpv_perm,file="/Users/shansun/Google\ Drive/picrust_paper/china_perm.csv")

#geography paper
mal_P1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/picrust_perc.txt",sep="\t")
mal_K1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/wgs_perc.txt",sep="\t")
mal_otu_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/otu_perc.txt",sep="\t")

adonis(t(mal_P1_p)~factor(c(rep(1,36),rep(2,48))))
"
                                  Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
factor(c(rep(1, 36), rep(2, 48)))  1    0.2799 0.279921  7.9314 0.08819  0.001 ***
Residuals                         82    2.8940 0.035293         0.91181           
Total                             83    3.1739                  1.00000           "

adonis(t(mal_K1_p)~factor(c(rep(1,36),rep(2,48))))

"                                  Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
factor(c(rep(1, 36), rep(2, 48)))  1    0.1954 0.195443  3.6116 0.04219  0.006 **
Residuals                         82    4.4374 0.054115         0.95781          
Total                             83    4.6329                  1.00000  "

adonis(t(mal_otu_p)~factor(c(rep(1,36),rep(2,48))))

"                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
factor(c(rep(1, 36), rep(2, 48)))  1    1.6288 1.62878  4.8627 0.05598  0.001 ***
Residuals                         82   27.4662 0.33495         0.94402           
Total                             83   29.0950                 1.00000 "


knames1=union(rownames(mal_P1_p),rownames(mal_K1_p))
Kpv=list()
Ppv=list()
i=1
for (m in knames1){
  if (m %in% rownames(mal_K1_p)){
    if (all(t.test(mal_K1_p[m,1:36],mal_K1_p[m,37:84])$statistic>0)){
      Kpv2=-log10(t.test(mal_K1_p[m,1:36],mal_K1_p[m,37:84])$p.value)
    }else{
      Kpv2=log10(t.test(mal_K1_p[m,1:36],mal_K1_p[m,37:84])$p.value)
    }
  }else{
    Kpv2=NA
  }
  if (m %in% rownames(mal_P1_p)){
    if (all(t.test(mal_P1_p[m,1:36],mal_P1_p[m,37:84])$statistic>0)){
      Ppv2=-log10(t.test(mal_P1_p[m,1:36],mal_P1_p[m,37:84])$p.value)
    }else{
      Ppv2=log10(t.test(mal_P1_p[m,1:36],mal_P1_p[m,37:84])$p.value)
    }
  }else{
    Ppv2=NA
  }
  Kpv[i]=Kpv2
  Ppv[i]=Ppv2
  i=i+1
}

ggs_K=cbind(unlist(Kpv),unlist(Ppv))
rownames(ggs_K)=knames1
colnames(ggs_K)=c("metagenome","PICRUSt")
write.csv(ggs_K,file="/Users/Shansun/Google\ Drive/picrust_paper/ggs_K.csv")

ggs_K1=ggs_K[ggs_K[,1]!=0&ggs_K[,2]!=0,]
ggs_K2=ggs_K[ggs_K[,1]==0|ggs_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/ggs_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(ggs_K1[,1],ggs_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Human 2",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

cor.test(ggs_K1[,1],ggs_K1[,2])
"
Pearson's product-moment correlation

data:  ggs_K1[, 1] and ggs_K1[, 2]
t = 40.464, df = 5185, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
0.4689336 0.5103078
sample estimates:
cor 
0.4898965  "

cor.test(ggs_K1[,1],ggs_K1[,2],method = "spearman")
"	Spearman's rank correlation rho

data:  ggs_K1[, 1] and ggs_K1[, 2]
S = 8.454e+09, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.553729 "

dim(ggs_K[ggs_K[,1]!=0&ggs_K[,2]!=0,])[1] #5187
dim(ggs_K[ggs_K[,1]!=0,])[1] #6022
dim(ggs_K[ggs_K[,2]!=0,])[1] #6230

pdf("/Users/shansun/Google\ Drive/picrust_paper/ggs_K_venn.pdf", width = 8, height = 8)
grid.newpage()
draw.pairwise.venn(6022, 6230, 5187, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.03, 0.07),cex=3,cat.cex=3)
dev.off()

#test spearman with the randomized metagenome data
comm3=intersect(rownames(mal_P1_p),rownames(mal_K1_p))
cor_ggs=vector()
i=1
for (m in 1:dim(mal_K1_p)[2]){
  cor1=cor.test(mal_P1_p[comm3,m],mal_K1_p[comm3,m],method="spearman")
  cor_ggs[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(mal_K1_p,null.model = "richness",iterations = 100)
  for (m in 1:dim(mal_K1_p)[2]){
    cor1=cor.test(mal_P1_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_ggs[i]=cor1$estimate
    i=i+1
  }
}
length(cor_ggs)
boxplot(cor_ggs[dim(mal_K1_p)[2]+1:length(cor_ggs)],ylim=c(min(cor_ggs)-0.2,max(cor_ggs)+0.2))
points(rep(1,dim(mal_K1_p)[2]),cor_ggs[1:dim(mal_K1_p)[2]],col="red")
t.test(cor_ggs[1:dim(mal_K1_p)[2]],cor_ggs[dim(mal_K1_p)[2]+1:length(cor_ggs)])

#test our method with randomized metagenome
comm2=union(rownames(mal_P1_p),rownames(mal_K1_p))
ggs_perm=matrix(nrow=length(comm2),ncol=100)
i=1
for (n in 1:100){
  wg_perm=randomizeMatrix(mal_K1_p,null.model = "richness",iterations = 100)
  i=1
  for (m in comm2){
    if (m %in% rownames(mal_K1_p)){
      if (all(t.test(wg_perm[m,1:36],wg_perm[m,37:84])$statistic>0)){
        Kpv2=-log10(t.test(wg_perm[m,1:36],wg_perm[m,37:84])$p.value)
      }else{
        Kpv2=log10(t.test(wg_perm[m,1:36],wg_perm[m,37:84])$p.value)
      }
    }else{
      Kpv2=NA
    }
    ggs_perm[i,n]=Kpv2
    i=i+1
  }
  print(n)
}
write.csv(ggs_perm,file="/Users/Shansun/Google\ Drive/picrust_paper/ggs_perm.csv")

#gorilla
gorK1=read.table(file='/Users/Shansun/Google Drive/picrust/gorilla/genefamilies.tsv',sep="\t",header=TRUE,row.names=1)
gorP1=read.table(file="/Users/Shansun/Google\ Drive/picrust/gorilla/metagenome_predictions.txt",sep="\t",header=TRUE,row.names=1)
gor_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust/gorilla/map_gorilla.txt",sep="\t",header=TRUE,row.names=1)
colnames(gorK1)=as.character(gor_meta$Sample_Name[match(colnames(gorK1),gor_meta$Run)])
gorK1=gorK1[,colnames(gorP1)]
gor_meta1=gor_meta[match(colnames(gorP1),as.character(gor_meta$Sample_Name)),]
table(sapply(strsplit(as.character(gor_meta1$collection_date),"-"),"[",2))

gorK1=gorK1[rowSums(gorK1)!=0,]
gorP1=gorP1[rowSums(gorP1)!=0,]
gorK1_p=t(t(gorK1)/colSums(gorK1))
gorP1_p=t(t(gorP1)/colSums(gorP1))
gorOTU=read.table(file="/Users/shansun/Google\ Drive/picrust/gorilla/otu_table_2640.txt",row.names=1,head=T)
gorOTU_p=t(t(gorOTU)/colSums(gorOTU))

meta_n=rep(1,14)
meta_n[c(3,5:8,10,12:14)]=2
adonis(t(gorK1_p)~factor(meta_n))
"              Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)
factor(meta_n)  1 0.0001028 0.00010284  0.2521 0.02058  0.538
Residuals      12 0.0048953 0.00040794         0.97942       
Total          13 0.0049982                    1.00000   "
adonis(t(gorP1_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
factor(meta_n)  1   0.01784 0.017844 0.52534 0.04194  0.754
Residuals      12   0.40759 0.033966         0.95806       
Total          13   0.42543                  1.00000  "
adonis(t(gorOTU_p)~factor(meta_n))
"               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
factor(meta_n)  1    0.3359 0.33591 0.99531 0.07659  0.437
Residuals      12    4.0500 0.33750         0.92341       
Total          13    4.3859                 1.00000  "


knames=union(rownames(gorK1_p),rownames(gorP1_p))
knames2=intersect(rownames(gorK1_p),rownames(gorP1_p))
length(knames)#5081
m=1
p_co=vector()
s_co=vector()
p_co_log=vector()
s_co_log=vector()
for (j in knames){
  if (!j %in% rownames(gorK1_p)){
    s_co_log[m]=NA
  }else{
    a1=t.test(gorK1_p[j,c(1,2,4,9,11)],gorK1_p[j,c(3,5:8,10,12:14)])
    s_co[m]=a1$p.value
    if (as.numeric(a1$statistic)<0){
      s_co_log[m]=-log10(s_co[m])
    }else{
      s_co_log[m]=log10(s_co[m])
    }
  }
  if (!j %in% rownames(gorP1_p)){
    p_co_log[m]=NA
  }else{
    a2=t.test(gorP1_p[j,c(1,2,4,9,11)],gorP1_p[j,c(3,5:8,10,12:14)])
    p_co[m]=a2$p.value
    if (as.numeric(a2$statistic)<0){
      p_co_log[m]=-log10(p_co[m])
    }else{
      p_co_log[m]=log10(p_co[m])
    }
  }
  m=m+1
}

gor_K=cbind(unlist(s_co_log),unlist(p_co_log))
rownames(gor_K)=knames
colnames(gor_K)=c("metagenome","PICRUSt")
write.csv(gor_K,file="/Users/Shansun/Google\ Drive/picrust_paper/gor_K.csv")
gor_K=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/gor_K.csv",row.names=1)

gor_K1=gor_K[gor_K[,1]!=0&gor_K[,2]!=0,]
gor_K2=gor_K[gor_K[,1]==0|gor_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/gor_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(gor_K1[,1],gor_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Gorilla",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

dim(gor_K[gor_K[,1]!=0&gor_K[,2]!=0,])[1] #3539
dim(gor_K[gor_K[,1]!=0,])[1] #3830
dim(gor_K[gor_K[,2]!=0,])[1] #4790

pdf("/Users/shansun/Google\ Drive/picrust_paper/gor_K_venn.pdf", width = 8, height = 8)
grid.newpage()
draw.pairwise.venn(3830, 4790, 3539, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.03, 0.05),cex=3,cat.cex=3)
dev.off()

cor.test(gor_K1[,1],gor_K1[,2])
"	Pearson's product-moment correlation

data:  gor_K1[, 1] and gor_K1[, 2]
t = 15.082, df = 3537, p-value < 2.2e-16
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
0.2146061 0.2765252
sample estimates:
cor 
0.2458164 "
cor.test(gor_K1[,1],gor_K1[,2],method="spearman")
"	Spearman's rank correlation rho

data:  gor_K1[, 1] and gor_K1[, 2]
S = 5076200000, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.3128601 "

comm3=intersect(rownames(gorK1_p),rownames(gorP1_p))
cor_gor=vector()
i=1
for (m in 1:dim(gorK1_p)[2]){
  cor1=cor.test(gorP1_p[comm3,m],gorK1_p[comm3,m],method="spearman")
  cor_gor[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(gorK1_p,null.model = "richness",iterations = 100)
  for (m in 1:dim(gorK1_p)[2]){
    cor1=cor.test(gorP1_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_gor[i]=cor1$estimate
    i=i+1
  }
}
length(cor_gor)
boxplot(cor_gor[dim(gorK1_p)[2]+1:length(cor_gor)],ylim=c(min(cor_gor)-0.2,max(cor_gor)+0.2))
points(rep(1,dim(gorK1_p)[2]),cor_gor[1:dim(gorK1_p)[2]],col="red")
t.test(cor_gor[1:dim(gorK1_p)[2]],cor_gor[dim(gorK1_p)[2]+1:length(cor_gor)])

comm3=union(rownames(gorK1_p),rownames(gorP1_p))
gor_perm=matrix(nrow=length(comm3),ncol=100)
for (n in 1:100){
  wg_perm=randomizeMatrix(gorK1_p,null.model = "richness",iterations = 100)
  i=1
  for (j in comm3){
    if (!j %in% rownames(gorK1_p)){
      Kpv2=NA
    }else{
      a1=t.test(wg_perm[j,c(1,2,4,9,11)],wg_perm[j,c(3,5:8,10,12:14)])
      Kpv=a1$p.value
      if (as.numeric(a1$statistic)<0){
        Kpv2=-log10(Kpv)
      }else{
        Kpv2=log10(Kpv)
      }
    }
    gor_perm[i,n]=Kpv2
    i=i+1
  }
}
write.csv(gor_perm,file="/Users/Shansun/Google\ Drive/picrust_paper/gor_perm.csv")


#mice 
pmeta_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/mice/picrust_perc.txt",sep="\t")
miceK_p=read.table(file='/Users/Shansun/Google Drive/picrust/mice/wgs_perc.txt',sep="\t")

otu_p=read.table(file='/Users/Shansun/Google Drive/picrust/mice/otu_perc.txt',sep="\t")

meta_n=rep(1,11)
meta_n[c(3:6,9)]=2
adonis(t(otu_p)~factor(meta_n))
"               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
factor(meta_n)  1   0.68338 0.68338  2.5677 0.22197  0.011 *
Residuals       9   2.39533 0.26615         0.77803         
Total          10   3.07870                 1.00000    "
adonis(t(pmeta_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
factor(meta_n)  1   0.09110 0.091100  5.9877 0.39951  0.004 **
Residuals       9   0.13693 0.015215         0.60049          
Total          10   0.22803                  1.00000   "
adonis(t(miceK_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
factor(meta_n)  1   0.04148 0.041485  1.1253 0.11114  0.291
Residuals       9   0.33178 0.036865         0.88886       
Total          10   0.37327                  1.00000   "

knames=union(rownames(miceK_p),rownames(pmeta_p))
m=1
p_co=vector()
s_co=vector()
p_co_log=vector()
s_co_log=vector()
for (j in knames){
  if (!j %in% rownames(miceK_p)){
    s_co_log[m]=NA
  }else{
    a1=t.test(miceK_p[j,c(1,2,7,8,10,11)],miceK_p[j,c(3:6,9)])
    s_co[m]=a1$p.value
    if (as.numeric(a1$statistic)<0){
      s_co_log[m]=-log10(s_co[m])
    }else{
      s_co_log[m]=log10(s_co[m])
    }
  }
  if (!j %in% rownames(pmeta_p)){
    p_co_log[m]=NA
  }else{
    a2=t.test(pmeta_p[j,c(1,2,7,8,10,11)],pmeta_p[j,c(3:6,9)])
    p_co[m]=a2$p.value
    if (as.numeric(a2$statistic)<0){
      p_co_log[m]=-log10(p_co[m])
    }else{
      p_co_log[m]=log10(p_co[m])
    }
  }
  m=m+1
}

mouse_K=cbind(unlist(s_co_log),unlist(p_co_log))
rownames(mouse_K)=knames
colnames(mouse_K)=c("metagenome","PICRUSt")
write.csv(mouse_K,file="/Users/Shansun/Google\ Drive/picrust_paper/mouse_K.csv")
mouse_K=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/mouse_K.csv",row.names=1)


mouse_K1=mouse_K[mouse_K[,1]!=0&mouse_K[,2]!=0,]
mouse_K2=mouse_K[mouse_K[,1]==0|mouse_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/mouse_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(mouse_K1[,1],mouse_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Mouse",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

cor.test(mouse_K1[,1],mouse_K1[,2])
"	Pearson's product-moment correlation

data:  mouse_K1[, 1] and mouse_K1[, 2]
t = 6.1575, df = 1362, p-value = 9.707e-10
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
0.1124760 0.2157639
sample estimates:
cor 
0.1645711 "

cor.test(mouse_K1[,1],mouse_K1[,2],method="spearman")
"	Spearman's rank correlation rho

data:  mouse_K1[, 1] and mouse_K1[, 2]
S = 335890000, p-value = 1.614e-14
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.2058562 "

dim(mouse_K[mouse_K[,1]!=0&mouse_K[,2]!=0,])[1] #1364
dim(mouse_K[mouse_K[,1]!=0,])[1] #1937
dim(mouse_K[mouse_K[,2]!=0,])[1] #4140

pdf("/Users/shansun/Google\ Drive/picrust_paper/mouse_K_venn.pdf", width = 8, height = 8)
grid.newpage()
venn.plot=draw.pairwise.venn(1937, 4140, 1364, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),cex=3,cat.cex=3)
grid.draw(venn.plot)
dev.off()

comm3=intersect(rownames(miceK_p),rownames(pmeta_p))
cor_mouse=vector()
i=1
for (m in 1:dim(miceK_p)[2]){
  cor1=cor.test(pmeta_p[comm3,m],miceK_p[comm3,m],method="spearman")
  cor_mouse[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(miceK_p,null.model = "richness",iterations = 100)
  for (m in 1:dim(miceK_p)[2]){
    cor1=cor.test(pmeta_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_mouse[i]=cor1$estimate
    i=i+1
  }
}
length(cor_mouse)
boxplot(cor_mouse[dim(miceK_p)[2]+1:length(cor_mouse)],ylim=c(min(cor_mouse)-0.2,max(cor_mouse)+0.2))
points(rep(1,dim(miceK_p)[2]),cor_mouse[1:dim(miceK_p)[2]],col="red")
t.test(cor_mouse[1:dim(miceK_p)[2]],cor_mouse[dim(miceK_p)[2]+1:length(cor_mouse)])

comm3=union(rownames(miceK_p),rownames(pmeta_p))
mouse_perm=matrix(nrow=length(comm3),ncol=100)
for (n in 1:100){
  wg_perm=randomizeMatrix(miceK_p,null.model = "richness",iterations = 100)
  i=1
  for (j in comm3){
    if (!j %in% rownames(miceK_p)){
      Kpv2=NA
    }else{
      a1=t.test(wg_perm[j,c(1,2,7,8,10,11)],wg_perm[j,c(3:6,9)])
      Kpv=a1$p.value
      if (as.numeric(a1$statistic)<0){
        Kpv2=-log10(Kpv)
      }else{
        Kpv2=log10(Kpv)
      }
    }
    mouse_perm[i,n]=Kpv2
    i=i+1
  }
}
write.csv(mouse_perm,file="/Users/Shansun/Google\ Drive/picrust_paper/mouse_perm.csv")

#chicken
chiK1=read.table(file='/Users/Shansun/Google Drive/picrust/chicken/genefamilies.tsv',sep="\t",header=TRUE,row.names=1)
chiP1=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/metagenome_predictions.txt",sep="\t",header=TRUE,row.names=1)
dim(chiP1)#6909   25
dim(chiK1)#4097   29
chi_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/map_chicken.txt",sep="\t",header=TRUE,row.names=1)
colnames(chiK1)=as.character(chi_meta$Sample_Name[match(colnames(chiK1),chi_meta$Run)])
chi_correct=read.csv(file="/Users/Shansun/Google\ Drive/picrust/chicken/sample_correct.csv")
colnames(chiP1)=as.character(chi_correct$sam[match(colnames(chiP1),chi_correct$seq)])

chiK1=chiK1[,sort(intersect(colnames(chiP1),colnames(chiK1)))]
chiP1=chiP1[,sort(intersect(colnames(chiP1),colnames(chiK1)))]
chi_meta1=chi_meta[match(colnames(chiP1),as.character(chi_meta$Sample_Name)),]

chiK1=chiK1[rowSums(chiK1)!=0,]
chiP1=chiP1[rowSums(chiP1)!=0,]
chiK1_p=t(t(chiK1)/colSums(chiK1))
chiP1_p=t(t(chiP1)/colSums(chiP1))

chiOTU=read.table(file="/Users/shansun/Google\ Drive/picrust/chicken/otu_table_1000.txt",sep="\t",header=TRUE,row.names=1)
colnames(chiOTU)=as.character(chi_correct$sam[match(colnames(chiOTU),chi_correct$seq)])
chiOTU_p=t(t(chiOTU)/colSums(chiOTU))

chiOTU_p=chiOTU_p[,order(colnames(chiOTU_p))]

meta_n=rep(1,25)
meta_n[c(1,12,15:23)]=2
adonis(t(chiOTU_p)~factor(meta_n))

"               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
factor(meta_n)  1    0.6801 0.68007  2.6652 0.10385  0.002 **
Residuals      23    5.8688 0.25517         0.89615          
Total          24    6.5489                 1.00000  "
adonis(t(chiP1_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
factor(meta_n)  1   0.05972 0.059719  3.6957 0.13844  0.012 *
Residuals      23   0.37165 0.016159         0.86156         
Total          24   0.43137                  1.00000 "

adonis(t(chiK1_p)~factor(meta_n))
"               Df SumsOfSqs    MeanSqs  F.Model      R2 Pr(>F)
factor(meta_n)  1  0.000071 0.00007084 0.046225 0.00201  0.945
Residuals      23  0.035248 0.00153251          0.99799       
Total          24  0.035319                     1.00000    "



dim(chiP1_p)#4864   25
dim(chiK1_p)#4085   25

knames=union(rownames(chiK1_p),rownames(chiP1_p))
knames2=intersect(rownames(chiK1_p),rownames(chiP1_p))
length(knames)#5224
m=1
p_co=vector()
s_co=vector()
p_co_log=vector()
s_co_log=vector()
for (j in knames){
  if (!j %in% rownames(chiK1_p)){
    s_co_log[m]=NA
  }else{
    a1=t.test(chiK1_p[j,c(2:11,13,14,24,25)],chiK1_p[j,c(1,12,15:23)])
    s_co[m]=a1$p.value
    if (as.numeric(a1$statistic)<0){
      s_co_log[m]=-log10(s_co[m])
    }else{
      s_co_log[m]=log10(s_co[m])
    }
  }
  if (!j %in% rownames(chiP1_p)){
    p_co_log[m]=NA
  }else{
    a2=t.test(chiP1_p[j,c(2:11,13,14,24,25)],chiP1_p[j,c(1,12,15:23)])
    p_co[m]=a2$p.value
    if (as.numeric(a2$statistic)<0){
      p_co_log[m]=-log10(p_co[m])
    }else{
      p_co_log[m]=log10(p_co[m])
    }
  }
  m=m+1
}

chi_K=cbind(unlist(s_co_log),unlist(p_co_log))
rownames(chi_K)=knames
colnames(chi_K)=c("metagenome","PICRUSt")
write.csv(chi_K,file="/Users/Shansun/Google\ Drive/picrust_paper/chi_K.csv")

chi_K1=chi_K[chi_K[,1]!=0&chi_K[,2]!=0,]
chi_K2=chi_K[chi_K[,1]==0|chi_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/chi_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(chi_K1[,1],chi_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Chicken",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

dim(chi_K[chi_K[,1]!=0&chi_K[,2]!=0,])[1] #3725
dim(chi_K[chi_K[,1]!=0,])[1] #4085
dim(chi_K[chi_K[,2]!=0,])[1] #4864

pdf("/Users/shansun/Google\ Drive/picrust_paper/chi_K_venn.pdf", width = 8, height = 8)
grid.newpage()
draw.pairwise.venn(4085, 4864, 3725, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.03, 0.05),cex=3,cat.cex=3)
dev.off()

cor.test(chi_K1[,1],chi_K1[,2])

"	Pearson's product-moment correlation

data:  chi_K1[, 1] and chi_K1[, 2]
t = -1.1402, df = 3723, p-value = 0.2543
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
-0.05076748  0.01344051
sample estimates:
cor 
-0.01868275 "

cor.test(chi_K1[,1],chi_K1[,2],method="spearman")
"	Spearman's rank correlation rho

data:  chi_K1[, 1] and chi_K1[, 2]
S = 8893500000, p-value = 0.04805
alternative hypothesis: true rho is not equal to 0
sample estimates:
        rho 
-0.03239407 "

comm3=intersect(rownames(chiK1_p),rownames(chiP1_p))
cor_chi=vector()
i=1
for (m in 1:dim(chiK1_p)[2]){
  cor1=cor.test(chiP1_p[comm3,m],chiK1_p[comm3,m],method="spearman")
  cor_chi[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(chiK1_p,null.model = "richness",iterations = 100)
  for (m in 1:dim(chiK1_p)[2]){
    cor1=cor.test(chiP1_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_chi[i]=cor1$estimate
    i=i+1
  }
}
length(cor_chi)
boxplot(cor_chi[dim(chiK1_p)[2]+1:length(cor_chi)],ylim=c(min(cor_chi)-0.2,max(cor_chi)+0.2))
points(rep(1,dim(chiK1_p)[2]),cor_chi[1:dim(chiK1_p)[2]],col="red")
t.test(cor_chi[1:dim(chiK1_p)[2]],cor_chi[dim(chiK1_p)[2]+1:length(cor_chi)])

comm2=union(rownames(chiK1_p),rownames(chiP1_p))
chi_perm=matrix(nrow=length(comm2),ncol=100)
for (n in 1:100){
  wg_perm=randomizeMatrix(chiK1_p,null.model = "richness",iterations = 100)
  i=1
  for (j in comm2){
    if (!j %in% rownames(chiK1_p)){
      Kpv2=NA
    }else{
      a1=t.test(wg_perm[j,c(2:11,13,14,24,25)],wg_perm[j,c(1,12,15:23)])
      Kpv=a1$p.value
      if (as.numeric(a1$statistic)<0){
        Kpv2=-log10(Kpv)
      }else{
        Kpv2=log10(Kpv)
      }
    }
    chi_perm[i,n]=Kpv2
    i=i+1
  }
}
write.csv(chi_perm,file="/Users/Shansun/Google\ Drive/picrust_paper/chicken_perm.csv")


#rhizosphere
rhizoK1_p=read.table(file='/Users/Shansun/Google Drive/picrust/rhizo/wgs_perc.txt',sep="\t")
rhizo_p1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/rhizo/picrust_perc.txt",sep="\t")
rhizo_otu_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/rhizo/otu_perc.txt",sep="\t")

meta_n=rep(1,27)
meta_n[c(4:6,17:27)]=2
adonis(t(rhizoK1_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
factor(meta_n)  1  0.035026 0.035026  6.5789 0.20833  0.001 ***
Residuals      25  0.133102 0.005324         0.79167           
Total          26  0.168129                  1.00000   "
adonis(t(rhizo_otu_p)~factor(meta_n))
"               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
factor(meta_n)  1    2.6861 2.68614  37.683 0.60117  0.001 ***
Residuals      25    1.7821 0.07128         0.39883           
Total          26    4.4682                 1.00000 "
adonis(t(rhizo_p1_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
factor(meta_n)  1   0.08934 0.089339  8.4973 0.25367  0.001 ***
Residuals      25   0.26284 0.010514         0.74633           
Total          26   0.35218                  1.00000           
"


knames=union(rownames(rhizoK1_p),rownames(rhizo_p1_p))
length(knames)#7976
m=1
p_co=vector()
s_co=vector()
p_co_log=vector()
s_co_log=vector()
for (j in knames){
  if (!j %in% rownames(rhizoK1_p)){
    s_co_log[m]=NA
  }else{
    a1=t.test(rhizoK1_p[j,c(1:3,7:16)],rhizoK1_p[j,c(4:6,17:27)])
    s_co[m]=a1$p.value
    if (as.numeric(a1$statistic)<0){
      s_co_log[m]=-log10(s_co[m])
    }else{
      s_co_log[m]=log10(s_co[m])
    }
  }
  if (!j %in% rownames(rhizo_p1_p)){
    p_co_log[m]=NA
  }else{
    a2=t.test(rhizo_p1_p[j,c(1:3,7:16)],rhizo_p1_p[j,c(4:6,17:27)])
    p_co[m]=a2$p.value
    if (as.numeric(a2$statistic)<0){
      p_co_log[m]=-log10(p_co[m])
    }else{
      p_co_log[m]=log10(p_co[m])
    }
  }
  m=m+1
}

rhizo_K=cbind(unlist(s_co_log),unlist(p_co_log))
rownames(rhizo_K)=knames
colnames(rhizo_K)=c("metagenome","PICRUSt")
write.csv(rhizo_K,file="/Users/Shansun/Google\ Drive/picrust_paper/rhizo_K.csv")

rhizo_K1=rhizo_K[rhizo_K[,1]!=0&rhizo_K[,2]!=0,]
rhizo_K2=rhizo_K[rhizo_K[,1]==0|rhizo_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/rhizo_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(rhizo_K1[,1],rhizo_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Soil 1",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

cor.test(rhizo_K1[,1],rhizo_K1[,2])
"	Pearson's product-moment correlation

data:  rhizo_K1[, 1] and rhizo_K1[, 2]
t = 7.6506, df = 2303, p-value = 2.921e-14
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
0.1173602 0.1969948
sample estimates:
cor 
0.1574334 "

cor.test(rhizo_K1[,1],rhizo_K1[,2],method="spearman")
"	Spearman's rank correlation rho

data:  rhizo_K1[, 1] and rhizo_K1[, 2]
S = 1.786e+09, p-value = 1.744e-09
alternative hypothesis: true rho is not equal to 0
sample estimates:
     rho 
0.124966 
"


dim(rhizo_K[rhizo_K[,1]!=0&rhizo_K[,2]!=0,])[1] #2305
dim(rhizo_K[rhizo_K[,1]!=0,])[1] #3518
dim(rhizo_K[rhizo_K[,2]!=0,])[1] #6051

pdf("/Users/shansun/Google\ Drive/picrust_paper/rhizo_K_venn.pdf", width = 8, height = 8)
grid.newpage()
draw.pairwise.venn(3518, 6051, 2305, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.03, 0.05),cex=3,cat.cex=3)
dev.off()

comm3=intersect(rownames(rhizoK1_p),rownames(rhizo_p1_p))
cor_rhizo=vector()
i=1
for (m in 1:dim(rhizoK1_p)[2]){
  cor1=cor.test(rhizo_p1_p[comm3,m],rhizoK1_p[comm3,m],method="spearman")
  cor_rhizo[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(rhizoK1_p,null.model = "richness",iterations = 100)
  for (m in 1:dim(rhizoK1_p)[2]){
    cor1=cor.test(rhizo_p1_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_rhizo[i]=cor1$estimate
    i=i+1
  }
}
length(cor_rhizo)
boxplot(cor_rhizo[dim(rhizoK1_p)[2]+1:length(cor_rhizo)],ylim=c(min(cor_rhizo)-0.2,max(cor_rhizo)+0.2))
points(rep(1,dim(rhizoK1_p)[2]),cor_rhizo[1:dim(rhizoK1_p)[2]],col="red")
t.test(cor_rhizo[1:dim(rhizoK1_p)[2]],cor_rhizo[dim(rhizoK1_p)[2]+1:length(cor_rhizo)])

pdf(file="/Users/Shansun/Google\ Drive/picrust_paper/soil1_exam.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(rhizo_p1_p[comm3,3],rhizoK1_p[comm3,3],xlab="Relative abundance from PICRUSt",ylab="Relative abundance from metagenome",pch=16,col="red",main="Unpermuted\nSpearman's rho = 0.851\n P < 2.2e-16")
wg_perm=randomizeMatrix(rhizoK1_p,null.model = "richness",iterations = 10000000000)
plot(rhizo_p1_p[comm3,3],wg_perm[comm3,3],xlab="Relative abundance from PICRUSt",ylab="Permuted relative abundance from metagenome",pch=16,col="blue",main="Permuted\nSpearman's rho = 0.837\n P < 2.2e-16")
dev.off()
cor.test(rhizo_p1_p[comm3,3],rhizoK1_p[comm3,3],method="spearman")
cor.test(rhizo_p1_p[comm3,3],wg_perm[comm3,3],method="spearman")

comm2=union(rownames(rhizoK1_p),rownames(rhizo_p1_p))
rhizo_perm=matrix(nrow=length(comm2),ncol=100)
for (n in 1:100){
  wg_perm=randomizeMatrix(rhizoK1_p,null.model = "richness",iterations = 100)
  i=1
  for (j in comm2){
    if (!j %in% rownames(rhizoK1_p)){
      Kpv2=NA
    }else{
      a1=t.test(wg_perm[j,c(1:3,7:16)],wg_perm[j,c(4:6,17:27)])
      Kpv=a1$p.value
      if (as.numeric(a1$statistic)<0){
        Kpv2=-log10(Kpv)
      }else{
        Kpv2=log10(Kpv)
      }
    }
    rhizo_perm[i,n]=Kpv2
    i=i+1
  }
}
write.csv(rhizo_perm,file="/Users/Shansun/Google\ Drive/picrust_paper/rhizo_perm.csv")




#deforestation
deforK1_p=read.table(file='/Users/Shansun/Google\ Drive/picrust/defor/metagenome_perc.txt',sep="\t")
defor_p1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/defor/picrust_perc.txt",sep="\t")

defor_otu_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/defor/otu_perc.txt",sep="\t")
meta_n=rep(1,14)
meta_n[c(7:14)]=2
adonis(t(defor_otu_p)~factor(meta_n))
"               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
factor(meta_n)  1   0.56258 0.56258  8.8463 0.42436  0.001 ***
Residuals      12   0.76314 0.06360         0.57564           
Total          13   1.32572                 1.00000 "

adonis(t(defor_p1_p)~factor(meta_n))
"
               Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
factor(meta_n)  1 0.0079842 0.0079842  7.5969 0.38766  0.001 ***
Residuals      12 0.0126116 0.0010510         0.61234           
Total          13 0.0205958                   1.00000 "

adonis(t(deforK1_p)~factor(meta_n))
"               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
factor(meta_n)  1   0.04497 0.044967  1.0376 0.07959  0.447
Residuals      12   0.52005 0.043337         0.92041       
Total          13   0.56502                  1.00000 "


knames=union(rownames(deforK1_p),rownames(defor_p1_p))
knames2=intersect(rownames(deforK1_p),rownames(defor_p1_p))
length(knames)#9386
m=1
p_co=vector()
s_co=vector()
p_co_log=vector()
s_co_log=vector()
for (j in knames){
  if (!j %in% rownames(deforK1_p)){
    s_co_log[m]=NA
  }else{
    a1=t.test(deforK1_p[j,c(1:6)],deforK1_p[j,c(7:14)])
    s_co[m]=a1$p.value
    if (as.numeric(a1$statistic)<0){
      s_co_log[m]=-log10(s_co[m])
    }else{
      s_co_log[m]=log10(s_co[m])
    }
  }
  if (!j %in% rownames(defor_p1_p)){
    p_co_log[m]=NA
  }else{
    a2=t.test(defor_p1_p[j,c(1:6)],defor_p1_p[j,c(7:14)])
    p_co[m]=a2$p.value
    if (as.numeric(a2$statistic)<0){
      p_co_log[m]=-log10(p_co[m])
    }else{
      p_co_log[m]=log10(p_co[m])
    }
  }
  m=m+1
}

defor_K=cbind(unlist(s_co_log),unlist(p_co_log))
rownames(defor_K)=knames
colnames(defor_K)=c("metagenome","PICRUSt")
write.csv(defor_K,file="/Users/Shansun/Google\ Drive/picrust_paper/defor_K.csv")

defor_K1=defor_K[defor_K[,1]!=0&defor_K[,2]!=0,]
defor_K2=defor_K[defor_K[,1]==0|defor_K[,2]==0,]
pdf("/Users/shansun/Google\ Drive/picrust_paper/defor_K.pdf", width = 8, height = 8)
par(mfrow=c(1,1),mar=c(5,5,3,3))
plot(defor_K1[,1],defor_K1[,2],xlab="metagenome",ylab="PICRUSt",pch=16,col="blue",main="Soil 2",cex.lab=2,cex.main=2)
abline(h=0,v=0,cex=1)
dev.off()

cor.test(defor_K1[,1],defor_K1[,2])
"	Pearson's product-moment correlation

data:  defor_K1[, 1] and defor_K1[, 2]
t = 6.469, df = 2172, p-value = 1.216e-10
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
0.09600151 0.17849483
sample estimates:
cor 
0.1374866 "

cor.test(defor_K1[,1],defor_K1[,2],method="spearman")
"	Spearman's rank correlation rho

data:  defor_K1[, 1] and defor_K1[, 2]
S = 1495200000, p-value = 2.904e-09
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.1268888 
"

dim(defor_K[defor_K[,1]!=0&defor_K[,2]!=0,])[1] #2174
dim(defor_K[defor_K[,1]!=0,])[1] #5814
dim(defor_K[defor_K[,2]!=0,])[1] #5746

pdf("/Users/shansun/Google\ Drive/picrust_paper/defor_K_venn.pdf", width = 8, height = 8)
grid.newpage()
draw.pairwise.venn(5814, 5746, 2174, category = c("metagenome", "PICRUSt"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.03, 0.05),cex=3,cat.cex=3)
dev.off()

comm3=intersect(rownames(deforK1_p),rownames(defor_p1_p))
cor_defor=vector()
i=1
for (m in 1:dim(deforK1_p)[2]){
  cor1=cor.test(defor_p1_p[comm3,m],deforK1_p[comm3,m],method="spearman")
  cor_defor[i]=cor1$estimate
  i=i+1
}
print(i-1)
for (n in 1:100){
  wg_perm=randomizeMatrix(deforK1_p,null.model = "richness",iterations = 100)
  for (m in 1:dim(deforK1_p)[2]){
    cor1=cor.test(defor_p1_p[comm3,m],wg_perm[comm3,m],method="spearman")
    cor_defor[i]=cor1$estimate
    i=i+1
  }
}
length(cor_defor)
boxplot(cor_defor[dim(deforK1_p)[2]+1:length(cor_defor)],ylim=c(min(cor_defor)-0.2,max(cor_defor)+0.2))
points(rep(1,dim(deforK1_p)[2]),cor_defor[1:dim(deforK1_p)[2]],col="red")
t.test(cor_defor[1:dim(deforK1_p)[2]],cor_defor[dim(deforK1_p)[2]+1:length(cor_defor)])

comm2=union(rownames(deforK1_p),rownames(defor_p1_p))
defor_perm=matrix(nrow=length(comm2),ncol=100)
for (n in 1:100){
  wg_perm=randomizeMatrix(deforK1_p,null.model = "richness",iterations = 100)
  i=1
  for (j in comm2){
    if (!j %in% rownames(deforK1_p)){
      Kpv2=NA
    }else{
      a1=t.test(wg_perm[j,c(1:6)],wg_perm[j,c(7:14)])
      Kpv=a1$p.value
      if (as.numeric(a1$statistic)<0){
        Kpv2=-log10(Kpv)
      }else{
        Kpv2=log10(Kpv)
      }
    }
    defor_perm[i,n]=Kpv2
    i=i+1
  }
}
write.csv(defor_perm,file="/Users/Shansun/Google\ Drive/picrust_paper/defor_perm.csv")



dataset=c("china_K","ggs_K","gor_K","mouse_K","chi_K","rhizo_K","defor_K")
var_lm=vector()
j=1
for (n in dataset){
  filename=paste("/Users/Shansun/Google\ Drive/picrust_paper/",n,".csv",sep="")
  P_tab=read.csv(file=filename,row.names=1)
  P_tab=P_tab[P_tab[,1]!=0&P_tab[,2]!=0,]
  var_lm[j]=summary(lm(P_tab[,1]~P_tab[,2]))$adj.r.squared
  j=j+1
}

#calculate the FP, FN and R2 of correlations in each category
b=read.table(file="/Users/Shansun/Google\ Drive/picrust/kegg_cat.txt",sep="\t")
cat2=names(sort(table(b[,3]),decreasing = T))
cat3=as.character(unique(b[,4]))

dataset=c("china_K","ggs_K","gor_K","mouse_K","chi_K","rhizo_K","defor_K")
rsquared=matrix(nrow=length(cat2),ncol=7)
pvalue=matrix(nrow=length(cat2),ncol=7)
perc_FN=matrix(nrow=length(cat2),ncol=7)
perc_FP=matrix(nrow=length(cat2),ncol=7)
cat5len=matrix(nrow=length(cat2),ncol=7)
spearman_r=matrix(nrow=length(cat2),ncol=7)
spearman_p=matrix(nrow=length(cat2),ncol=7)
i=1
for (n in dataset){
  filename=paste("/Users/Shansun/Google\ Drive/picrust_paper/",n,".csv",sep="")
  P_tab=read.csv(file=filename,row.names=1)
  j=0
  for (m in cat2){
    j=j+1
    cat5=b[!is.na(match(as.character(b[,3]),m)),5]
    pos=na.omit(match(cat5,rownames(P_tab)))
    cat5len[j,i]=length(na.omit(match(cat5,rownames(P_tab[!is.na(P_tab[,1])&!is.na(P_tab[,2]),]))))
    if (all(is.na(pos))) next
    if (length(pos)==1) next
    perc_FN[j,i]=length(which(is.na(P_tab[pos,2])))/length(which(!is.na(P_tab[pos,1])))
    perc_FP[j,i]=length(which(is.na(P_tab[pos,1])))/length(which(!is.na(P_tab[pos,2])))
    adjP6=P_tab[pos,]
    adjP7=as.matrix(adjP6[!is.na(adjP6[,1]) & !is.na(adjP6[,2]),])
    if (is.null(dim(adjP7))) next
    if (dim(adjP7)[1]<3) next
    if (sd(unlist(adjP7[,1]))==0|sd(unlist(adjP7[,2]))==0) next
    rsquared[j,i]=summary(lm(adjP7[,1]~adjP7[,2]))$adj.r.squared
    pvalue[j,i]=summary(lm(adjP7[,1]~adjP7[,2]))$coefficients[2,4]
    spearman_r[j,i]=cor.test(adjP7[,1],adjP7[,2],method="spearman")$estimate
    spearman_p[j,i]=cor.test(adjP7[,1],adjP7[,2],method="spearman")$p.value
  }
  i=i+1
}

rownames(spearman_p)=cat2
rownames(spearman_r)=cat2
rownames(perc_FP)=cat2
rownames(perc_FN)=cat2

write.csv(cat5len,file="/Users/Shansun/Google\ Drive/picrust_paper/category_cat5len.csv")
write.csv(perc_FN,file="/Users/Shansun/Google\ Drive/picrust_paper/category_perc_FN.csv")
write.csv(perc_FP,file="/Users/Shansun/Google\ Drive/picrust_paper/category_perc_FP.csv")
write.csv(rsquared,file="/Users/Shansun/Google\ Drive/picrust_paper/category_rsquared.csv")
write.csv(pvalue,file="/Users/Shansun/Google\ Drive/picrust_paper/category_pvalue.csv")
write.csv(spearman_r,file="/Users/Shansun/Google\ Drive/picrust_paper/category_rho.csv")
write.csv(spearman_p,file="/Users/Shansun/Google\ Drive/picrust_paper/category_spearman_p.csv")


rownames(cat5len)=cat2
rownames(perc_FN)=cat2
rownames(perc_FP)=cat2
rownames(rsquared)=cat2
rownames(pvalue)=cat2

cat5len=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/category_cat5len.csv",row.names=1)
perc_FN=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/category_perc_FN.csv",row.names=1)
perc_FP=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/category_perc_FP.csv",row.names=1)
rsquared=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/category_rsquared.csv",row.names=1)
pvalue=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/category_pvalue.csv",row.names=1)

rsquared_sig=rsquared
rsquared_sig[pvalue>=0.05]=0

level1=b[match(cat2,b[,3]),2]
col_bar=c("pink","coral","plum2","yellow","lightgreen","lightblue","lightgray")[factor(level1)]

legend_names=substring(levels(factor(level1)),7)
pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_legend.pdf",height=3,width=20)
plot.new()
legend("center",legend = legend_names, fill=c("pink","coral","plum2","yellow","lightgreen","lightblue","lightgray"), cex=1, horiz = TRUE,bty = "n",text.width=c(0.1,0.12,0.11,-0.0,0.16,0.115,0.1))
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_rsquared_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(rsquared_sig[,1]),names.arg=rev(cat2),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1)
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_rsquared_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(rsquared_sig[,i]),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5)
}
dev.off()

spearman_r1=spearman_r
spearman_r1[spearman_p>0.05]=0
spearman_r1[spearman_r1<0]=0
pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_rho_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(spearman_r1[,1]),names.arg=rev(cat2),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1)
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_rho_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(spearman_r1[,i]),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5)
}
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_fp_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FP[,1]),names.arg=rev(cat2),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1)
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_fp_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(perc_FP[,i]),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5)
}
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust_paper/category_fn_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(perc_FN[,i]),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5)
}
dev.off()



spear_cor=matrix(nrow=8484,ncol=7)
spear_cor[,2]=cor_ggs
spear_cor[1:40,1]=cor_china[1:40]
spear_cor[85:(84+4000),1]=cor_china[41:4040]

spear_cor[1:14,3]=cor_gor[1:14]
spear_cor[85:(84+1400),3]=cor_gor[15:1414]

spear_cor[1:11,4]=cor_mouse[1:11]
spear_cor[85:(84+1100),4]=cor_mouse[12:1111]

spear_cor[1:25,5]=cor_chi[1:25]
spear_cor[85:(84+2500),5]=cor_chi[26:2525]

spear_cor[1:27,6]=cor_rhizo[1:27]
spear_cor[85:(84+2700),6]=cor_rhizo[28:2727]

spear_cor[1:14,7]=cor_defor[1:14]
spear_cor[85:(84+1400),7]=cor_defor[15:1414]

colnames(spear_cor)=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")
write.csv(spear_cor,file="/Users/Shansun/Google\ Drive/picrust_paper/spearman_cor_random.csv")
spear_cor=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/spearman_cor_random.csv",row.names=1)
spear_cor1=cbind(spear_cor,c(rep("unpermuted",84),rep("permuted",8400)))
colnames(spear_cor1)[8]="group"
spear_cor_m=melt(spear_cor1)
colnames(spear_cor_m)=c("Group","Dataset","rho")
library(reshape2)
library(ggpubr)
pdf("/Users/Shansun/Google\ Drive/picrust_paper/spearman_cor_random.pdf",height=5,width=12)
ggboxplot(spear_cor_m, x = "Dataset", y = "rho", color = "Group", palette = c("blue","red"), add = "jitter",  ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),
                                     axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
dev.off()




dataset=c("china","ggs","gor","mouse","chi","rhizo","defor")
perm_cor=matrix(nrow=100,ncol=7)
perm_pvalue=matrix(nrow=100,ncol=7)
j=1
for (n in dataset){
  filename1=paste("/Users/Shansun/Google\ Drive/picrust_paper/perms/",n,"_perm.csv",sep="")
  filename2=paste("/Users/Shansun/Google\ Drive/picrust_paper/",n,"_K.csv",sep="")
  P_tab1=read.csv(file=filename1,row.names=1)
  P_tab2=read.csv(file=filename2,row.names=1)
  for (i in 1:100){
    perm_cor[i,j]=cor.test(as.numeric(as.character(P_tab1[,i])),as.numeric(as.character(P_tab2[,2])),method="spearman")$estimate
    perm_pvalue[i,j]=cor.test(as.numeric(as.character(P_tab1[,i])),as.numeric(as.character(P_tab2[,2])),method="spearman")$p.value
  }
  j=j+1
}

colnames(perm_cor)=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")
write.csv(perm_cor,file="/Users/Shansun/Google\ Drive/picrust_paper/pvalue_cor_spearman_random.csv")
colMeans(perm_cor)
"-0.0004855118 -0.0008332769 -0.0010257308 -0.0006016839  0.0178251801  0.0001459619 -0.0036350653 "

spear_cor1=cbind(spear_cor,c(rep("unpermuted",84),rep("permuted",8400)))
colnames(spear_cor1)[8]="group"
perm_cor_m=melt(perm_cor)
perm_cor_m=perm_cor_m[,-1]
colnames(perm_cor_m)=c("Dataset","rho")
pdf("/Users/Shansun/Google\ Drive/picrust_paper/pvalue_cor_random_spearman.pdf",height=6,width=10)
plot1=ggboxplot(perm_cor_m, x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
plot2=plot1+geom_point(aes(x="Human_KW",y=0.490, color="red"), show.legend=TRUE)+geom_point(aes(x="Human_TY",y=0.554, color="red"))+geom_point(aes(x="Gorilla",y=0.313, color="red"))+geom_point(aes(x="Mouse",y=0.206, color="red"))+geom_point(aes(x="Chicken",y=-0.0324, color="red"))+geom_point(aes(x="Soil_LWM",y=0.126, color="red"))+geom_point(aes(x="Soil_AAN",y=0.127, color="red"))
ggpar(plot2, ylim = c(-0.05,0.6))
dev.off()


pdf("/Users/Shansun/Google\ Drive/picrust_paper/pvalue_cor_random_spearman.pdf",height=6,width=8)
boxplot(perm_cor,ylim=c(0,0.6),col="blue",ylab="Spearman's rho of transformed P-values",xlab="Datasets")
points(c(1:7),c(0.490,0.554,0.313,0.206,-0.0324,0.125,0.127),col="red",pch=16)
legend("topright",c("permuted","unpermuted"),col=c("lightblue","red"),pch=16,bty="n",cex=1.5)
dev.off()


#barplot showing the consistency of class and function categories across all 7 datasets
#merge all taxonomy tables at class level
dataset=c("china","ggs","gorilla","mice","chicken","rhizo","defor")
P_tab1=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/otu_table_L3.txt",header=T)
P_tab1=P_tab1[,c(1,order(colnames(P_tab1)[2:dim(P_tab1)[2]])+1)]
P_tab1=cbind(as.character(P_tab1[,1]),t(t(P_tab1[,-1])/colSums(P_tab1[,-1])))
colnames(P_tab1)[1]="OTUID"
sample_n=vector()
i=1
sample_n[i]=dim(P_tab1)[2]
for (n in dataset[-1]){
  filename=paste("/Users/Shansun/Google\ Drive/picrust/",n,"/otu_table_L3.txt",sep="")
  P_tab_a=read.table(file=filename,header=T)
  P_tab_a=P_tab_a[,c(1,order(colnames(P_tab_a)[2:dim(P_tab_a)[2]])+1)]
  P_tab=cbind(as.character(P_tab_a[,1]),t(t(P_tab_a[,-1])/colSums(P_tab_a[,-1])))
  colnames(P_tab)[1]="OTUID"
  P_tab1=merge(P_tab1,P_tab,by="OTUID",all=T)
  print(dim(P_tab))
  print(dim(P_tab1))
  i=i+1
  sample_n[i]=dim(P_tab1)[2]
}
rownames(P_tab1)=P_tab1[,1]
P_tab1=P_tab1[,-1]
P_tab1[is.na(P_tab1)]=0

#subset the table of all 16S samples to only keep those also with metagenomes
china_names=c(1:40)

mal_PK=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/metagenome_predictions.txt",header=TRUE,sep="\t",row.names=1)
mal_K=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/GGS-kom110_kegg58_matrix.txt",header=TRUE,sep="\t",row.names=1)
cnames=sort(intersect(colnames(mal_PK),colnames(mal_K)))
colnames(P_tab1)[81:176]=sapply(strsplit(colnames(P_tab1)[81:176],"\\.4"), "[[", 1)
ggs_names=na.omit(match(cnames,colnames(P_tab1)))

gor_names=c(177:190)

mice_meta=read.table(file="/Users/shansun/Google\ Drive/picrust/mice/seqs_loc.txt",header=T)
colnames(P_tab1)[191:202]=sapply(strsplit(colnames(P_tab1)[191:202],"\\.fastq"), "[[", 1)
colnames(P_tab1)[191:202]=sapply(strsplit(as.character(mice_meta$library_name[match(colnames(P_tab1)[191:202],mice_meta$run_accession)]),"_W0"), "[[", 1)
  
chiK1=read.table(file='/Users/Shansun/Google Drive/picrust/chicken/genefamilies.tsv',sep="\t",header=TRUE,row.names=1)
chiP1=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/metagenome_predictions.txt",sep="\t",header=TRUE,row.names=1)
chi_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/map_chicken.txt",sep="\t",header=TRUE,row.names=1)
colnames(chiK1)=as.character(chi_meta$Sample_Name[match(colnames(chiK1),chi_meta$Run)])
chi_correct=read.csv(file="/Users/Shansun/Google\ Drive/picrust/chicken/sample_correct.csv")
colnames(chiP1)=as.character(chi_correct$sam[match(colnames(chiP1),chi_correct$seq)])
colnames(P_tab1)[203:227]=colnames(chiP1)
chicken_names=na.omit(match(intersect(colnames(chiP1),colnames(chiK1)),colnames(P_tab1)))

colnames(P_tab1)[228:231]=c("BulkAG1", "BulkAG2", "BulkAG3", "BulkAG4" )
rhizo_names=c(228:267)

deforK1=read.table(file='/Users/Shansun/Google Drive/picrust/defor/defor_mgK',sep="\t")
defor_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust/defor/map.txt",sep="\t",header=TRUE,row.names=1)
colnames(P_tab1)[268:285]=as.character(defor_meta$Alias[match(colnames(P_tab1)[268:285],rownames(defor_meta))])
com_row=intersect(colnames(deforK1),colnames(P_tab1)[268:285])
defor_names=match(com_row,colnames(P_tab1))

P_tab2=P_tab1[,c(1:40,ggs_names,gor_names,mice_names,chicken_names,rhizo_names,defor_names)]

#merge all functional profiles of the 7 datasets
wg=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/kegg_complete_perc.txt",sep="\t",stringsAsFactors=F,header=T,row.names=NULL)
mal_K1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/metagenome_perc.txt",sep="\t",row.names=NULL)
gorK1=read.table(file='/Users/Shansun/Google Drive/picrust/gorilla/genefamilies.tsv',sep="\t",header=TRUE,row.names=NULL)
miceK_p=read.table(file='/Users/Shansun/Google Drive/picrust/mice/metagenome_perc.txt',sep="\t",row.names=NULL)
chiK1=read.table(file='/Users/Shansun/Google Drive/picrust/chicken/genefamilies.tsv',sep="\t",header=TRUE,row.names=NULL)
rhizoK1_p=read.table(file='/Users/Shansun/Google Drive/picrust/rhizo/metagenome_perc.txt',sep="\t",row.names=NULL)
deforK1_p=read.table(file='/Users/Shansun/Google Drive/picrust/defor/metagenome_perc.txt',sep="\t",row.names=NULL)
dim(wg) #13880    41
dim(mal_K1_p) #6022   85
dim(gorK1) #3925   20
dim(miceK_p) #1937   12
dim(chiK1) #4097   30
dim(rhizoK1_p) #3518   28
dim(deforK1_p) #5814   15

all_K=Reduce(function(x, y) merge(x, y, all=TRUE,by=1), list(wg, mal_K1_p, gorK1,miceK_p,chiK1,rhizoK1_p,deforK1_p))
rownames(all_K)=all_K[,1]
all_K=all_K[,-1]
all_K=all_K[-match("UNMAPPED",rownames(all_K)),]#remove unmapped
all_K[is.na(all_K)]=0
all_K1=t(t(all_K)/colSums(all_K))
dim(all_K1) #13937   224

colnames(all_K1)[1:40]=colnames(P_tab1)[1:40]

gor_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust/gorilla/map_gorilla.txt",sep="\t",header=TRUE,row.names=1)
colnames(all_K1)[125:143]=as.character(gor_meta$Sample_Name[match(colnames(all_K1)[125:143],gor_meta$Run)])

chi_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/map_chicken.txt",sep="\t",header=TRUE,row.names=1)
colnames(all_K1)[155:183]=as.character(chi_meta$Sample_Name[match(colnames(all_K1)[155:183],chi_meta$Run)])

pw_names=intersect(colnames(all_K1),colnames(P_tab2))
all_K1=all_K1[,pw_names]
P_tab2=P_tab2[,pw_names]
c(1:40,41:124,125:138,139:149,150:174,175:201,202:215)

P_tab3=apply(P_tab2,2,as.numeric)
colnames(P_tab3)=c(1:dim(P_tab3)[2])
rownames(P_tab3)=rownames(P_tab2)
P_tab3=P_tab3[order(rowMeans(P_tab3),decreasing = TRUE),]
P_tab4=rbind(P_tab3[1:20,],1-colSums(P_tab3[1:20,]))
rownames(P_tab4)[21]="Other"
P_tab4[,1:40]=P_tab4[,1:40][,order(P_tab4[1,1:40])]
P_tab4[,41:124]=P_tab4[,41:124][,order(P_tab4[1,41:124])]
P_tab4[,125:138]=P_tab4[,125:138][,order(P_tab4[1,125:138])]
P_tab4[,139:149]=P_tab4[,139:149][,order(P_tab4[1,139:149])]
P_tab4[,150:174]=P_tab4[,150:174][,order(P_tab4[1,150:174])]
P_tab4[,175:201]=P_tab4[,175:201][,order(P_tab4[1,175:201])]
P_tab4[,202:215]=P_tab4[,202:215][,order(P_tab4[1,202:215])]

#plot class barplot
pdf(file="/Users/shansun/Google\ Drive/picrust/class_bar.pdf",height=15,width=120)
par(mar=c(5,8,5,2),mfrow=c(1,2))
barplot(P_tab4[,-211],col=c(rainbow(1),rainbow(20,start=0.08)),space=0, border = NA,axes=F)
axis(side = 2, pos = -1,cex.axis=2)
plot.new()
legend("left",c("Other",sapply(strsplit(rev(rownames(P_tab4))[-1],"c__"),"[[",2)),fill=rev(c(rainbow(1),rainbow(20,start=0.08))),cex=2,bty="n")
dev.off()

all_K2=cbind(all_K1,as.character(b[match(rownames(all_K1),b[,5]),3]))
colnames(all_K2)[dim(all_K2)[2]]="cat2"
all_K3=aggregate(all_K1, list(all_K2[,216]), FUN=sum)
rownames(all_K3)=all_K3[,1]
all_K3=all_K3[,-1]
colnames(all_K3)=c(1:215)

all_K4=apply(all_K3,2,as.numeric)
rownames(all_K4)=rownames(all_K3)
all_K4=all_K4[order(rowMeans(all_K4),decreasing = TRUE),]
all_K5=rbind(all_K4[1:15,],1-colSums(all_K4[1:15,]))
rownames(all_K5)[16]="Other"

all_K5[,1:40]=all_K5[,1:40][,order(all_K5[1,1:40])]
all_K5[,41:124]=all_K5[,41:124][,order(all_K5[1,41:124])]
all_K5[,125:138]=all_K5[,125:138][,order(all_K5[1,125:138])]
all_K5[,139:149]=all_K5[,139:149][,order(all_K5[1,139:149])]
all_K5[,150:174]=all_K5[,150:174][,order(all_K5[1,150:174])]
all_K5[,175:201]=all_K5[,175:201][,order(all_K5[1,175:201])]
all_K5[,202:215]=all_K5[,202:215][,order(all_K5[1,202:215])]


#barplot of the functional category at level2
pdf(file="/Users/shansun/Google\ Drive/picrust/cat2_bar.pdf",height=15,width=120)
par(mar=c(5,8,5,2),mfrow=c(1,2))
barplot(all_K5[,-202],col=rainbow(16),space=0, border = NA,offset=0,axes=F)
axis(side = 2, pos = -1,cex.axis=2)
plot.new()
legend("left",rev(rownames(all_K5)),fill=rev(rainbow(16)),cex=2,bty="n")
dev.off()

#sequencing depth
#China
otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/otu_table_L3.txt",header=T,row.names=1)
mean(colSums(otu)) #67776.9
metagenome=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/kegg_complete_file.txt",header=T,row.names=1)
mean(colSums(metagenome)) #5055905

otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/otu_table_L3.txt",header=T,row.names=1)
mean(colSums(otu)) #1789817
#metagenome sequence depth from literature 155,890

otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/gorilla/otu_table.txt",header=T,row.names=1)
mean(colSums(otu)) #11905.19
metagenome=read.table(file="/Users/Shansun/Google\ Drive/picrust/gorilla/genefamilies.tsv",header=T,row.names=1)
mean(colSums(metagenome)) #27825038

otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/mice/otu_table.txt",header=T,row.names=1)
mean(colSums(otu)) #2328.667
metagenome=read.table(file="/Users/Shansun/Google\ Drive/picrust/mice/function.txt",sep="\t")
metagenome1=metagenome[-c(1:4),4:16]
metagenome2=apply(metagenome1,2,as.numeric)
mean(colSums(metagenome2))#184088.5

otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/otu_table.txt",header=T,row.names=1)
mean(colSums(otu)) #4869.345
metagenome=read.table(file="/Users/Shansun/Google\ Drive/picrust/chicken/genefamilies.tsv",header=T,row.names=1)
mean(colSums(metagenome)) #31138616

otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/rhizo/otu_table.txt",header=T,row.names=1)
mean(colSums(otu)) #25981.5
metagenome=read.table(file="/Users/Shansun/Google\ Drive/picrust/rhizo/fun.txt",sep="\t",header=T)
mean(colSums(metagenome[,-c(1:4)])) #207174.5

otu=read.table(file="/Users/Shansun/Google\ Drive/picrust/defor/otu_table.txt",header=T,row.names=1)
mean(colSums(otu)) #8824.778
metagenome=read.table(file="/Users/Shansun/Google\ Drive/picrust/defor/fun.txt",sep="\t",header=T)
mean(colSums(metagenome[,-c(1:4)])) #1498014







