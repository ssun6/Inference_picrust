#china
pi=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/metagenome_predictions.txt",header=TRUE,sep="\t",row.names=1)
pi=pi[,sort(colnames(pi))]
pi_p1=t(t(pi)/colSums(pi))[,1:40]
pi_p=pi_p1[rowSums(pi_p1)!=0,]
wg=read.table(file="/Users/Shansun/Google\ Drive/picrust/china/kegg_complete_perc.txt",sep="\t",stringsAsFactors=F,header=T,row.names=1)
wg=wg[,sort(colnames(wg))]

#sample size 40:20+20

#orthology
dim(wg)[1] #13880
dim(pi_p)[1] #5612
comm2=Reduce(union, list(rownames(wg),rownames(pi_p)))
comm3=Reduce(intersect, list(rownames(wg),rownames(pi_p)))
uni=setdiff(rownames(wg),rownames(pi_p))
length(comm2)#13918
length(comm3)#5574

Kpv=matrix(nrow=length(comm2),ncol=100)
Ppv=matrix(nrow=length(comm2),ncol=100)
for (j in 1:100){
  i=1
  x=c(1:20)[sample.int(20, 5, replace = FALSE)]
  y=c(21:40)[sample.int(20, 5, replace = FALSE)]
  for (m in comm2){
    if ((m %in% rownames(wg))){
      if (any(wg[m,x]!=0)&any(wg[m,y]!=0)){
        if (t.test(wg[m,x],wg[m,y])$statistic>0){
          Kpv2=-log10(t.test(wg[m,x],wg[m,y])$p.value)
        }else{
          Kpv2=log10(t.test(wg[m,x],wg[m,y])$p.value)
        }
      }else{
        Kpv2=NA
      }
    }else{
      Kpv2=NA
    }
    if ((m %in% rownames(pi_p))){
      if (any(pi_p[m,x]!=0)&any(pi_p[m,y]!=0)){
        if (t.test(pi_p[m,x],pi_p[m,y])$statistic>0){
          Ppv2=-log10(t.test(pi_p[m,x],pi_p[m,y])$p.value)
        }else{
          Ppv2=log10(t.test(pi_p[m,x],pi_p[m,y])$p.value)
        }
      }else{
        Ppv2=NA
      }
    }else{
      Ppv2=NA
    }
    Kpv[i,j]=Kpv2
    Ppv[i,j]=Ppv2
    i=i+1
  }
}
write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/china_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/china_P_subsample.csv")

Kpv=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/china_K_subsample.csv",row.names=1)
Ppv=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/china_P_subsample.csv",row.names=1)

p_t_china=vector()
r_t_china=vector()
for (n in 1:100){
  p_t_china[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_china[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}


#geography paper
#sample size 84:36+48
mal_P1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/picrust_perc.txt",sep="\t")
mal_K1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/ggs/wgs_perc.txt",sep="\t")

knames1=union(rownames(mal_P1_p),rownames(mal_K1_p))

Kpv=matrix(nrow=length(knames1),ncol=100)
Ppv=matrix(nrow=length(knames1),ncol=100)
for (j in 1:100){
  i=1
  x=c(1:36)[sample.int(36, 5, replace = FALSE)]
  y=c(37:84)[sample.int(48, 5, replace = FALSE)]
  for (m in knames1){
    if (m %in% rownames(mal_K1_p)){
      if (any(mal_K1_p[m,x]!=0)&any(mal_K1_p[m,y]!=0)){
        if (all(t.test(mal_K1_p[m,x],mal_K1_p[m,y])$statistic>0)){
          Kpv2=-log10(t.test(mal_K1_p[m,x],mal_K1_p[m,y])$p.value)
        }else{
          Kpv2=log10(t.test(mal_K1_p[m,x],mal_K1_p[m,y])$p.value)
        }
      }else{
        Kpv2=NA
      }
    }else{
      Kpv2=NA
    }
    if (m %in% rownames(mal_P1_p)){
      if (any(mal_P1_p[m,x]!=0)&any(mal_P1_p[m,y]!=0)){
        if (all(t.test(mal_P1_p[m,x],mal_P1_p[m,y])$statistic>0)){
          Ppv2=-log10(t.test(mal_P1_p[m,x],mal_P1_p[m,y])$p.value)
        }else{
          Ppv2=log10(t.test(mal_P1_p[m,x],mal_P1_p[m,y])$p.value)
        }
      }else{
        Ppv2=NA
      }
    }else{
      Ppv2=NA
    }
    Kpv[i,j]=Kpv2
    Ppv[i,j]=Ppv2
    i=i+1
  }
  print(j)
}

write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/ggs_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/ggs_P_subsample.csv")

p_t_ggs=vector()
r_t_ggs=vector()
for (n in 1:100){
  p_t_ggs[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_ggs[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}

#gorilla
#sample size 14:5+9
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

knames=union(rownames(gorK1_p),rownames(gorP1_p))
knames2=intersect(rownames(gorK1_p),rownames(gorP1_p))
length(knames)#5081

Kpv=matrix(nrow=length(knames),ncol=100)
Ppv=matrix(nrow=length(knames),ncol=100)
for (j in 1:100){
  i=1
  x=c(1,2,4,9,11)[sample.int(5, 5, replace = FALSE)]
  y=c(3,5:8,10,12:14)[sample.int(9, 5, replace = FALSE)]
  for (m in knames){
    if (!m %in% rownames(gorK1_p)){
      Kpv[i,j]=NA
    }else{
      if (any(gorK1_p[m,x]!=0)&any(gorK1_p[m,y]!=0)){
        a1=t.test(gorK1_p[m,c(1,2,4,9,11)],gorK1_p[m,c(3,5:8,10,12:14)])
        s_co=a1$p.value
        if (as.numeric(a1$statistic)<0){
          Kpv[i,j]=-log10(s_co)
        }else{
          Kpv[i,j]=log10(s_co)
        }
      }else{
        Kpv[i,j]=NA
      }
    }
    if (!m %in% rownames(gorP1_p)){
      Ppv[i,j]=NA
    }else{
      if (any(gorP1_p[m,x]!=0)&any(gorP1_p[m,y]!=0)){
        a2=t.test(gorP1_p[m,c(1,2,4,9,11)],gorP1_p[m,c(3,5:8,10,12:14)])
        p_co=a2$p.value
        if (as.numeric(a2$statistic)<0){
          Ppv[i,j]=-log10(p_co)
        }else{
          Ppv[i,j]=log10(p_co)
        }
      }else{
        Ppv[i,j]=NA
      }
    }
    i=i+1
  }
  print(j)
}

write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/gor_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/gor_P_subsample.csv")

Kpv=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/gor_K_subsample.csv",row.names=1)
Ppv=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/gor_P_subsample.csv",row.names=1)

p_t_gor=vector()
r_t_gor=vector()
for (n in 1:100){
  p_t_gor[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_gor[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}

  
#mouse sample size 11:6+5
pmeta_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/mice/picrust_perc.txt",sep="\t")
miceK_p=read.table(file='/Users/Shansun/Google Drive/picrust/mice/wgs_perc.txt',sep="\t")

knames=union(rownames(miceK_p),rownames(pmeta_p))

Kpv=matrix(nrow=length(knames),ncol=6)
Ppv=matrix(nrow=length(knames),ncol=6)
for (j in 1:6){
  i=1
  x=c(1,2,7,8,10,11)[sample.int(6, 5, replace = FALSE)]
  y=c(3:6,9)[sample.int(5, 5, replace = FALSE)]
  for (m in knames){
    if (!m%in% rownames(miceK_p)){
      Kpv[i,j]=NA
    }else{
      if (any(miceK_p[m,x]!=0)&any(miceK_p[m,y]!=0)){
        a1=t.test(miceK_p[m,x],miceK_p[m,y])
        s_co=a1$p.value
        if (as.numeric(a1$statistic)<0){
          Kpv[i,j]=-log10(s_co)
        }else{
          Kpv[i,j]=log10(s_co)
        }
      }else{
        Kpv[i,j]=NA
      }
    }
    if (!m %in% rownames(pmeta_p)){
      Ppv[i,j]=NA
    }else{
      if (any(pmeta_p[m,x]!=0)&any(pmeta_p[m,y]!=0)){
        a2=t.test(pmeta_p[m,x],pmeta_p[m,y])
        p_co=a2$p.value
        if (as.numeric(a2$statistic)<0){
          Ppv[i,j]=-log10(p_co)
        }else{
          Ppv[i,j]=log10(p_co)
        }
      }else{
        Ppv[i,j]=NA
      }
    }
    i=i+1
  }
}

write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/mouse_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/mouse_P_subsample.csv")

p_t_mouse=vector()
r_t_mouse=vector()
for (n in 1:6){
  p_t_mouse[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_mouse[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}


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

dim(chiP1_p)#4864   25
dim(chiK1_p)#4085   25

knames=union(rownames(chiK1_p),rownames(chiP1_p))
knames2=intersect(rownames(chiK1_p),rownames(chiP1_p))
length(knames)#5224

Kpv=matrix(nrow=length(knames),ncol=100)
Ppv=matrix(nrow=length(knames),ncol=100)
for (j in 61:100){
  i=1
  x=c(2:11,13,14,24,25)[sample.int(14, 5, replace = FALSE)]
  y=c(1,12,15:23)[sample.int(11, 5, replace = FALSE)]
  for (m in knames){
    if (!m %in% rownames(chiK1_p)){
      Kpv[i,j]=NA
    }else{
      if (any(chiK1_p[m,x]!=0)&any(chiK1_p[m,y]!=0)){
        a1=t.test(chiK1_p[m,x],chiK1_p[m,y])
        s_co=a1$p.value
        if (as.numeric(a1$statistic)<0){
          Kpv[i,j]=-log10(s_co)
        }else{
          Kpv[i,j]=log10(s_co)
        }
      }else{
        Kpv[i,j]=NA
      }
    }
    if (!m %in% rownames(chiP1_p)){
      Ppv[i,j]=NA
    }else{
      if (any(chiP1_p[m,x]!=0)&any(chiP1_p[m,y]!=0)){
        a2=t.test(chiP1_p[m,x],chiP1_p[m,y])
        p_co=a2$p.value
        if (as.numeric(a2$statistic)<0){
          Ppv[i,j]=-log10(p_co)
        }else{
          Ppv[i,j]=log10(p_co)
        }
      }else{
        Ppv[i,j]=NA
      }
    }
    i=i+1
  }
  print(j)
}
  
write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/chicken_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/chicken_P_subsample.csv")

p_t_chicken=vector()
r_t_chicken=vector()
for (n in 1:100){
  p_t_chicken[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_chicken[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}

#rhizosphere sample size 27 (13+14)
rhizoK1_p=read.table(file='/Users/Shansun/Google Drive/picrust/rhizo/wgs_perc.txt',sep="\t")
rhizo_p1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/rhizo/picrust_perc.txt",sep="\t")

knames=union(rownames(rhizoK1_p),rownames(rhizo_p1_p))
length(knames)#7976

Kpv=matrix(nrow=length(knames),ncol=100)
Ppv=matrix(nrow=length(knames),ncol=100)
for (j in 1:100){
  i=1
  x=c(1:3,7:16)[sample.int(13, 5, replace = FALSE)]
  y=c(4:6,17:27)[sample.int(14, 5, replace = FALSE)]
  for (m in knames){
    if (!m %in% rownames(rhizoK1_p)){
      Kpv[i,j]=NA
    }else{
      if (any(rhizoK1_p[m,x]!=0)&any(rhizoK1_p[m,y]!=0)){
        a1=t.test(rhizoK1_p[m,x],rhizoK1_p[m,y])
        s_co=a1$p.value
        if (as.numeric(a1$statistic)<0){
          Kpv[i,j]=-log10(s_co)
        }else{
          Kpv[i,j]=log10(s_co)
        }
      }else{
        Kpv[i,j]=NA
      }
    }
    if (!m %in% rownames(rhizo_p1_p)){
      Ppv[i,j]=NA
    }else{
      if (any(rhizo_p1_p[m,x]!=0)&any(rhizo_p1_p[m,y]!=0)){
        a2=t.test(rhizo_p1_p[m,x],rhizo_p1_p[m,y])
        p_co=a2$p.value
        if (as.numeric(a2$statistic)<0){
          Ppv[i,j]=-log10(p_co)
        }else{
          Ppv[i,j]=log10(p_co)
        }
      }else{
        Ppv[i,j]=NA
      }
    }
    i=i+1
  }
  print(j)
}

write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/rhizo_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/rhizo_P_subsample.csv")

p_t_rhizo=vector()
r_t_rhizo=vector()
for (n in 1:100){
  p_t_rhizo[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_rhizo[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}

#deforestation
deforK1_p=read.table(file='/Users/Shansun/Google Drive/picrust/defor/wgs_perc.txt',sep="\t")
defor_p1_p=read.table(file="/Users/Shansun/Google\ Drive/picrust/defor/picrust_perc.txt",sep="\t")
knames=union(rownames(deforK1_p),rownames(defor_p1_p))
length(knames)#9386


Kpv=matrix(nrow=length(knames),ncol=56)
Ppv=matrix(nrow=length(knames),ncol=56)
for (j in 1:56){
  i=1
  x=c(1:6)[sample.int(6, 5, replace = FALSE)]
  y=c(7:14)[sample.int(8, 5, replace = FALSE)]
  for (m in knames){
    if (!m %in% rownames(deforK1_p)){
      Kpv[i,j]=NA
    }else{
      if (any(deforK1_p[m,x]!=0)&any(deforK1_p[m,y]!=0)){
        a1=t.test(deforK1_p[m,x],deforK1_p[m,y])
        s_co=a1$p.value
        if (as.numeric(a1$statistic)<0){
          Kpv[i,j]=-log10(s_co)
        }else{
          Kpv[i,j]=log10(s_co)
        }
      }else{
        Kpv[i,j]=NA
      }
    }
    if (!m %in% rownames(defor_p1_p)){
      Ppv[i,j]=NA
    }else{
      if (any(defor_p1_p[m,x]!=0)&any(defor_p1_p[m,y]!=0)){
        a2=t.test(defor_p1_p[m,x],defor_p1_p[m,y])
        p_co=a2$p.value
        if (as.numeric(a2$statistic)<0){
          Ppv[i,j]=-log10(p_co)
        }else{
          Ppv[i,j]=log10(p_co)
        }
      }else{
        Ppv[i,j]=NA
      }
    }
    i=i+1
  }
  print(j)
}

write.csv(Kpv,file="/Users/Shansun/Google\ Drive/picrust_paper/defor_K_subsample.csv")
write.csv(Ppv,file="/Users/Shansun/Google\ Drive/picrust_paper/defor_P_subsample.csv")

p_t_defor=vector()
r_t_defor=vector()
for (n in 1:56){
  p_t_defor[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$p.value
  r_t_defor[n]=cor.test(Kpv[,n], Ppv[,n],method="spearman")$estimate
}



subsample_cor=matrix(nrow=100,ncol=7)
subsample_cor[1:100,1]=r_t_china
subsample_cor[1:100,2]=r_t_ggs
subsample_cor[1:100,3]=r_t_gor
subsample_cor[1:6,4]=r_t_mouse
subsample_cor[1:100,5]=r_t_chicken
subsample_cor[1:100,6]=r_t_rhizo
subsample_cor[1:56,7]=r_t_defor


colnames(subsample_cor)=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")
write.csv(subsample_cor,file="/Users/Shansun/Google\ Drive/picrust_paper/subsample_spearman_cor_random.csv")

pdf("/Users/Shansun/Google\ Drive/picrust_paper/subsample_pvalue_cor_random.pdf",height=6,width=8)
boxplot(subsample_cor,ylim=c(-0.2,1),col="lightblue",ylab="Spearman's rho of transformed P-values",xlab="Datasets")
points(c(1:7),c(0.490,0.554,0.313,0.206,-0.0324,0.125,0.127),col="red",pch=16)
legend("topright",c("All","subsampled"),col=c("red","lightblue"),pch=15,bty="n",cex=1.5)
dev.off()

perm_subsample_cor_m=melt(subsample_cor)
perm_subsample_cor_m=perm_subsample_cor_m[,-1]
colnames(perm_subsample_cor_m)=c("Dataset","rho")
pdf("/Users/Shansun/Google\ Drive/picrust_paper/subsample_pvalue_cor_random.pdf",height=6,width=10)
plot1=ggboxplot(perm_subsample_cor_m, x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
plot2=plot1+geom_point(aes(x="Human_KW",y=0.490, color="red"), show.legend=TRUE)+geom_point(aes(x="Human_TY",y=0.554, color="red"))+geom_point(aes(x="Gorilla",y=0.313, color="red"))+geom_point(aes(x="Mouse",y=0.206, color="red"))+geom_point(aes(x="Chicken",y=-0.0324, color="red"))+geom_point(aes(x="Soil_LWM",y=0.126, color="red"))+geom_point(aes(x="Soil_AAN",y=0.127, color="red"))
ggpar(plot2, ylim = c(-0.5,1))
dev.off()

spear_cor=read.csv(file="/Users/Shansun/Google\ Drive/picrust_paper/spearman_cor_random.csv",row.names=1)
spear_cor1=cbind(spear_cor,c(rep("real data",84),rep("randomized data",8400)))
colnames(spear_cor1)[8]="group"
spear_cor_m=melt(spear_cor1)
colnames(spear_cor_m)=c("Group","Dataset","rho")
library(reshape2)
library(ggpubr)
pdf("/Users/Shansun/Google\ Drive/picrust_paper/spearman_cor_random.pdf",height=5,width=12)
ggboxplot(spear_cor_m, x = "Dataset", y = "rho", color = "Group", palette = c("blue","red"), add = "jitter",  ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
       panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),
       axis.title=element_text(size=15),legend.text=element_text(size=15),,legend.title=element_text(size=16)))
dev.off()




