library(Tax4Fun)
library(VennDiagram)
library(ggplot2)
list1=c("china","ggs","gorilla","chicken","mouse","rhizo","defor")
for (n in list1[1:5]){
  a=importQIIMEData(paste0("/Users/shansun/Google Drive/picrust2/",n,"/otu_table.txt"))
  b=Tax4Fun(a,"/Users/shansun/Google Drive/picrust2/defor/4fun/SILVA123")
  write.csv(b$Tax4FunProfile,file=paste0("/Users/shansun/Google Drive/picrust2/",n,"/",n,"4fun.csv"))
}

#parse files for Tax4Fun
china4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/china/china4fun.csv",header=T,row.names=1)
china4fun=t(china4fun1)
china4fun=china4fun[,order(colnames(china4fun))]
china4fun=t(t(china4fun)/colSums(china4fun))[,1:40]
china4fun=china4fun[rowSums(china4fun)!=0,]
colnames(china4fun)=sub("first", "", colnames(china4fun))
colnames(china4fun)=sub("A", "", colnames(china4fun))
colnames(china4fun)=paste0("X", formatC(as.numeric(colnames(china4fun)), width=3, flag="0"))
china4fun=china4fun[,order(colnames(china4fun))]
write.csv(china4fun,file="/Users/shansun/Google Drive/picrust2/china/china4fun_n.csv")

ggs4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/ggs/ggs4fun.csv",header=T,row.names=1)
ggs4fun=t(ggs4fun1)
ggs4fun=ggs4fun[,order(colnames(ggs4fun))]
colnames(ggs4fun)=sapply(strsplit(colnames(ggs4fun),"\\.418"),"[[",1)
write.csv(ggs4fun,file="/Users/shansun/Google Drive/picrust2/ggs/ggs4fun_n.csv")

gor4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/gorilla/gorilla4fun.csv",header=T,row.names=1)
gor4fun=t(gor4fun1)
gor4fun=gor4fun[,order(colnames(gor4fun))]
write.csv(gor4fun,file="/Users/shansun/Google Drive/picrust2/gorilla/gorilla4fun_n.csv")

chi4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/chicken/chicken4fun.csv",header=T,row.names=1)
chi4fun=t(chi4fun1)
chi_correct=read.csv(file="/Users/Shansun/Google\ Drive/picrust/chicken/sample_correct.csv")
colnames(chi4fun)=as.character(chi_correct$sam[match(colnames(chi4fun),chi_correct$seq)])
chi4fun=chi4fun[,order(colnames(chi4fun))]
write.csv(chi4fun,file="/Users/shansun/Google Drive/picrust2/chicken/chicken4fun_n.csv")

mouse4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/mouse/mouse4fun.csv",header=T,row.names=1)
mouse4fun=t(mouse4fun1)
colnames(mouse4fun)=sapply(strsplit(colnames(mouse4fun),"\\.fastq"),"[[",1)
mouse_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust2/mouse/seqs_loc.txt",sep="\t",header=TRUE)
mouse_meta$sra=sapply(strsplit(as.character(mouse_meta$sra_galaxy),"/"),"[[",5)
match(colnames(mouse4fun),mouse_meta$sra)
colnames(mouse4fun)=sapply(strsplit(as.character(mouse_meta$library_name[match(colnames(mouse4fun),mouse_meta$sra)]),"_W0"),"[[",1)
write.csv(mouse4fun,file="/Users/shansun/Google Drive/picrust2/mouse/mouse4fun_n.csv")


rhizo4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/rhizo/rhizo4fun.csv",header=T,row.names=1)
rhizo4fun=t(rhizo4fun1)
rhizo4fun=rhizo4fun[,colnames(rhizoK1_p)]
write.csv(rhizo4fun,file="/Users/shansun/Google Drive/picrust2/rhizo/rhizo4fun_n.csv")

defor4fun1=read.csv(file="/Users/shansun/Google Drive/picrust2/defor/defor4fun.csv",header=T,row.names=1)
defor4fun=t(defor4fun1)
defor_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust2/defor/map.txt",sep="\t",header=TRUE,row.names=1)
colnames(defor4fun)=defor_meta$Alias[match(colnames(defor4fun),rownames(defor_meta))]
defor4fun=defor4fun[,sort(colnames(defor4fun))]
write.csv(defor4fun,file="/Users/shansun/Google Drive/picrust2/defor/defor4fun_n.csv")


#Parse files for Picrust2
china_p2=read.table(file="/Users/shansun/Google Drive/picrust2/china/china_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
china_p2=china_p2[,order(colnames(china_p2))]
china_p2=t(t(china_p2)/colSums(china_p2))[,1:40]
china_p2=china_p2[rowSums(china_p2)!=0,]
colnames(china_p2)=sub("first", "", colnames(china_p2))
colnames(china_p2)=sub("A", "", colnames(china_p2))
colnames(china_p2)=paste0("X", formatC(as.numeric(colnames(china_p2)), width=3, flag="0"))
china_p2=china_p2[,order(colnames(china_p2))]
write.csv(china_p2,file="/Users/shansun/Google Drive/picrust2/china/china_p2_n.csv")

ggs_p2=read.table(file="/Users/shansun/Google Drive/picrust2/ggs/ggs_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
ggs_p2=ggs_p2[,order(colnames(ggs_p2))]
colnames(ggs_p2)=sapply(strsplit(colnames(ggs_p2),"\\.418"),"[[",1)
write.csv(ggs_p2,file="/Users/shansun/Google Drive/picrust2/ggs/ggs_p2_n.csv")

gor_p2=read.table(file="/Users/shansun/Google Drive/picrust2/gorilla/gorilla_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
gor_p2=gor_p2[,order(colnames(gor_p2))]
write.csv(gor_p2,file="/Users/shansun/Google Drive/picrust2/gorilla/gorilla_p2_n.csv")

chi_p2=read.table(file="/Users/shansun/Google Drive/picrust2/chicken/chicken_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
chi_correct=read.csv(file="/Users/Shansun/Google\ Drive/picrust/chicken/sample_correct.csv")
colnames(chi_p2)=as.character(chi_correct$sam[match(colnames(chi_p2),chi_correct$seq)])
chi_p2=chi_p2[,order(colnames(chi_p2))]
write.csv(chi_p2,file="/Users/shansun/Google Drive/picrust2/chicken/chicken_p2_n.csv")

mouse_p2=read.table(file="/Users/shansun/Google Drive/picrust2/mouse/mouse_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
colnames(mouse_p2)=sapply(strsplit(colnames(mouse_p2),"\\.fastq"),"[[",1)
mouse_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust2/mouse/seqs_loc.txt",sep="\t",header=TRUE)
mouse_meta$sra=sapply(strsplit(as.character(mouse_meta$sra_galaxy),"/"),"[[",5)
match(colnames(mouse_p2),mouse_meta$sra)
colnames(mouse_p2)=sapply(strsplit(as.character(mouse_meta$library_name[match(colnames(mouse_p2),mouse_meta$sra)]),"_W0"),"[[",1)
write.csv(mouse_p2,file="/Users/shansun/Google Drive/picrust2/mouse/mouse_p2_n.csv")

rhizo_p2=read.csv(file="/Users/shansun/Google\ Drive/picrust2/rhizo/rhizo_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
rhizo_p2=rhizo_p2[,colnames(rhizo_p2)]
write.csv(rhizo_p2,file="/Users/shansun/Google Drive/picrust2/rhizo/rhizo_p2_n.csv")

defor_p2=read.table(file="/Users/shansun/Google Drive/picrust2/defor/defor_pred_metagenome_unstrat.tsv",header=T,row.names=1,sep="\t")
defor_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust2/defor/map.txt",sep="\t",header=TRUE,row.names=1)
colnames(defor_p2)=defor_meta$Alias[match(colnames(defor_p2),rownames(defor_meta))]
defor_p2=defor_p2[,sort(colnames(defor_p2))]
write.csv(defor_p2,file="/Users/shansun/Google Drive/picrust2/defor/defor_p2_n.csv")

#Parse files for picrust
china_p1=read.table(file="/Users/shansun/Google Drive/picrust2/china/metagenome_predictions.txt",header=T,row.names=1,sep="\t")
china_p1=china_p1[,order(colnames(china_p1))]
colnames(china_p1)=gsub("A","X",colnames(china_p1))
china_p1=t(t(china_p1)/colSums(china_p1))[,1:40]
china_p1=china_p1[rowSums(china_p1)!=0,]
write.csv(china_p1,file="/Users/shansun/Google Drive/picrust2/china/china_p1_n.csv")

ggs_p1=read.table(file="/Users/shansun/Google Drive/picrust2/ggs/metagenome_predictions.txt",header=T,row.names=1,sep="\t")
ggs_p1=ggs_p1[,order(colnames(ggs_p1))]
ggs_p1=t(t(ggs_p1)/colSums(ggs_p1))
write.csv(ggs_p1,file="/Users/shansun/Google Drive/picrust2/ggs/ggs_p1_n.csv")

gorP1=read.table(file="/Users/Shansun/Google\ Drive/picrust2/gorilla/metagenome_predictions.txt",sep="\t",header=TRUE,row.names=1)
gorP1=gorP1[,order(colnames(gorP1))]
gorP1=t(t(gorP1)/colSums(gorP1))
write.csv(gorP1,file="/Users/shansun/Google Drive/picrust2/gorilla/gorilla_p1_n.csv")

chi_p1=read.table(file="/Users/shansun/Google Drive/picrust2/chicken/metagenome_predictions.txt",header=T,row.names=1,sep="\t")
chi_correct=read.csv(file="/Users/Shansun/Google\ Drive/picrust/chicken/sample_correct.csv")
colnames(chi_p1)=as.character(chi_correct$sam[match(colnames(chi_p1),chi_correct$seq)])
chi_p1=chi_p1[,order(colnames(chi_p1))]
chi_p1=t(t(chi_p1)/colSums(chi_p1))
write.csv(chi_p1,file="/Users/shansun/Google Drive/picrust2/chicken/chicken_p1_n.csv")

mouse_p1=read.table(file="/Users/shansun/Google Drive/picrust2/mouse/metagenome_predictions.txt",header=T,row.names=1,sep="\t")
mouse_p1=mouse_p1[,order(colnames(mouse_p1))]
mouse_p1=t(t(mouse_p1)/colSums(mouse_p1))
write.csv(mouse_p1,file="/Users/shansun/Google Drive/picrust2/mouse/mouse_p1_n.csv")

rhizo_p1=read.csv(file="/Users/shansun/Google\ Drive/picrust2/rhizo/metagenome_predictions.txt",header=T,row.names=1,sep="\t")
rhizo_p1=rhizo_p1[,colnames(rhizo_p1)]
rhizo_p1=t(t(rhizo_p1)/colSums(rhizo_p1))
write.csv(rhizo_p1,file="/Users/shansun/Google Drive/picrust2/rhizo/rhizo_p1_n.csv")

defor_p1=read.table(file="/Users/shansun/Google Drive/picrust2/defor/metagenome_predictions.txt",header=T,row.names=1,sep="\t")
defor_meta=read.table(file="/Users/Shansun/Google\ Drive/picrust2/defor/map.txt",sep="\t",header=TRUE,row.names=1)
colnames(defor_p1)=defor_meta$Alias[match(colnames(defor_p1),rownames(defor_meta))]
defor_p1=defor_p1[,sort(colnames(defor_p1))]
defor_p1=t(t(defor_p1)/colSums(defor_p1))
write.csv(defor_p1,file="/Users/shansun/Google Drive/picrust2/defor/defor_p1_n.csv")



#inference calculation
setwd("/Users/shansun/Google Drive/picrust2/")

group1=list()
group2=list()
group1[["china"]]=c(1:20)
group2[["china"]]=c(21:40)
group1[["ggs"]]=c(1:36)
group2[["ggs"]]=c(37:84)
group1[["gorilla"]]=c(1:3,5:6)
group2[["gorilla"]]=c(4,7:14)
group1[["chicken"]]=c(2:11,13,14,24,25)
group2[["chicken"]]=c(1,12,15:23)
group1[["mouse"]]=c(1,2,7,8,10,11)
group2[["mouse"]]=c(3:6,9)
group1[["rhizo"]]=c(1:3,7:16)
group2[["rhizo"]]=c(4:6,17:27)
group1[["defor"]]=c(1:6)
group2[["defor"]]=c(7:14)

list2=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")

K_cor=matrix(ncol=9,nrow=21)
T_cor=matrix(ncol=9,nrow=21)

test_m=c("_p1","_p2","4fun")

j=1
for (n1 in test_m){
  for ( n in list1){
    n_P2_p=read.csv(file=paste0(n,"/",n,n1,"_n.csv"),header=T,row.names=1)
    rownames(n_P2_p)=sapply(strsplit(rownames(n_P2_p),"\\.\\."),"[[",1)
    n_K1_p=read.table(file=paste0(n,"/",n,"_wgs_p.txt"),sep="\t")
    
    cnames=sort(intersect(colnames(n_P2_p),colnames(n_K1_p)))
    knames=union(rownames(n_P2_p),rownames(n_K1_p))
    n_K1_p=n_K1_p[,cnames]
    n_P1_p=n_P2_p[,cnames]
    n_P1_p=n_P1_p[rowSums(n_P1_p)!=0,]
    n_K1_p=n_K1_p[rowSums(n_K1_p)!=0,]
    
    knames1=union(rownames(n_P1_p),rownames(n_K1_p))
    Kpv=list()
    Ppv=list()
    Kstat=list()
    Pstat=list()
    i=1
    for (m in knames1){
      if (m %in% rownames(n_K1_p)){
        if (all(t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$statistic>0)){
          Kpv2=-log10(t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$p.value)
          Kstat2=t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$statistic
        }else{
          Kpv2=log10(t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$p.value)
          Kstat2=t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$statistic
        }
      }else{
        Kpv2=NA
        Kstat2=NA
      }
      if (m %in% rownames(n_P1_p)){
        if (all(t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$statistic>0)){
          Ppv2=-log10(t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$p.value)
          Pstat2=t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$statistic
        }else{
          Ppv2=log10(t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$p.value)
          Pstat2=t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$statistic
        }
      }else{
        Ppv2=NA
        Pstat2=NA
      }
      Kpv[i]=Kpv2
      Ppv[i]=Ppv2
      Kstat[i]=Kstat2
      Pstat[i]=Pstat2
      i=i+1
    }
    
    n_K=cbind(unlist(Kpv),unlist(Ppv))
    n_T=cbind(unlist(Kstat),unlist(Pstat))
    
    rownames(n_K)=knames1
    colnames(n_K)=c("metagenome","prediction")
    write.csv(n_K,file=paste0("/Users/shansun/Google Drive/picrust2/",n,"/",n,n1,"_K.csv"))
    
    rownames(n_T)=knames1
    colnames(n_T)=c("metagenome","prediction")
    write.csv(n_T,file=paste0("/Users/shansun/Google Drive/picrust2/",n,"/",n,n1,"_T.csv"))
    
    n_K1=n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),]#2237
    n_K2=n_K[is.na(n_K[,1])|is.na(n_K[,2]),]#7939
    pdf(paste0("/Users/shansun/Google Drive/picrust2/",n,"/",n,n1,"_K.pdf"), width = 8, height = 8)
    par(mfrow=c(1,1),mar=c(5,5,3,3))
    plot(n_K1[,1],n_K1[,2],xlab="metagenome",ylab="prediction",pch=16,col="blue",main=list2[j],cex.lab=2,cex.main=2)
    abline(h=0,v=0,cex=1)
    dev.off()
    
    n_T1=n_T[!is.na(n_T[,1])&!is.na(n_T[,2]),]#2237
    n_T2=n_T[is.na(n_T[,1])|is.na(n_T[,2]),]#7939
    pdf(paste0("/Users/shansun/Google Drive/picrust2/",n,"/",n,n1,"_T.pdf"), width = 8, height = 8)
    par(mfrow=c(1,1),mar=c(5,5,3,3))
    plot(n_T1[,1],n_T1[,2],xlab="metagenome",ylab="prediction",pch=16,col="blue",main=list2[j],cex.lab=2,cex.main=2)
    abline(h=0,v=0,cex=1)
    dev.off()
    
    a=dim(n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),])[1]
    b=dim(n_K[!is.na(n_K[,1]),])[1]
    c=dim(n_K[!is.na(n_K[,2]),])[1] 
    pdf(paste0("/Users/shansun/Google Drive/picrust2/",n,"/",n,n1,"_venn.pdf"), width = 8, height = 8)
    grid.newpage()
    draw.pairwise.venn(b, c, a, category = c("metagenome", "prediction"), lty = rep("blank", 2), fill = c("lightblue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.03, 0.07),cex=3,cat.cex=3)
    dev.off()
    
    K_cor[j,1]=cor.test(n_K1[,1],n_K1[,2])$estimate
    K_cor[j,2]=cor.test(n_K1[,1],n_K1[,2],method='spearman')$estimate
    K_cor[j,3]=cor.test(n_K1[,1],n_K1[,2],method='kendall')$estimate
    K_cor[j,4]=cor.test(n_K1[,1],n_K1[,2])$p.value
    K_cor[j,5]=cor.test(n_K1[,1],n_K1[,2],method='spearman')$p.value
    K_cor[j,6]=cor.test(n_K1[,1],n_K1[,2],method='kendall')$p.value
    K_cor[j,7]=dim(n_K[!is.na(n_K[,1]),])[1]
    K_cor[j,8]=dim(n_K[!is.na(n_K[,2]),])[1]
    K_cor[j,9]=dim(n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),])[1]
    
    T_cor[j,1]=cor.test(n_T1[,1],n_T1[,2])$estimate
    T_cor[j,2]=cor.test(n_T1[,1],n_T1[,2],method='spearman')$estimate
    T_cor[j,3]=cor.test(n_T1[,1],n_T1[,2],method='kendall')$estimate
    T_cor[j,4]=cor.test(n_T1[,1],n_T1[,2])$p.value
    T_cor[j,5]=cor.test(n_T1[,1],n_T1[,2],method='spearman')$p.value
    T_cor[j,6]=cor.test(n_T1[,1],n_T1[,2],method='kendall')$p.value
    T_cor[j,7]=dim(n_T[!is.na(n_T[,1]),])[1]
    T_cor[j,8]=dim(n_T[!is.na(n_T[,2]),])[1]
    T_cor[j,9]=dim(n_T[!is.na(n_T[,1])&!is.na(n_T[,2]),])[1]
    
    j=j+1
  }
}

write.csv(T_cor,file="/Users/shansun/Google\ Drive/picrust2/ttest_rho.csv")

inter_num=data.frame(cbind(c(rep("PICRUSt",7),rep("PICRUSt2",7),rep("Tax4Fun",7)),c(rep(list2,3)),1-T_cor[,9]/T_cor[,7]))
colnames(inter_num)=c("method","dataset","Percentage")
inter_num[,3]=as.numeric(as.character(inter_num[,3]))
inter_num[,2]=factor(inter_num[,2],levels=list2)
pdf("/Users/shansun/Google\ Drive/picrust2/bar_false_neg.pdf", width = 10, height = 5)
ggplot(inter_num, aes(dataset, Percentage),ggtheme = theme_bw()) +  geom_bar(aes(fill = method), position = "dodge", stat="identity")+
  theme(strip.text.x = element_text(size = 20, angle = 0),
        legend.position="right",
        axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()

inter_num=data.frame(cbind(c(rep("PICRUSt",7),rep("PICRUSt2",7),rep("Tax4Fun",7)),c(rep(list2,3)),1-T_cor[,9]/T_cor[,8]))
colnames(inter_num)=c("method","dataset","Percentage")
inter_num[,3]=as.numeric(as.character(inter_num[,3]))
inter_num[,2]=factor(inter_num[,2],levels=list2)
pdf("/Users/shansun/Google\ Drive/picrust2/bar_false_pos.pdf", width = 10, height = 5)
ggplot(inter_num, aes(dataset, Percentage),ggtheme = theme_bw()) +  geom_bar(aes(fill = method), position = "dodge", stat="identity")+
  theme(strip.text.x = element_text(size = 20, angle = 0),
        legend.position="right",
        axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()


spear_group=matrix(nrow=21,ncol=4)
spear_group[1:7,1]=T_cor_p1[,2]
spear_group[8:14,1]=T_cor_p2[,2]
spear_group[15:21,1]=T_cor[,2]
spear_group[1:7,2]=rep("PICRUSt",7)
spear_group[8:14,2]=rep("PICRUSt2",7)
spear_group[15:21,2]=rep("Tax4Fun",7)
spear_group[,3]=rep(list2,3)
spear_group[1:7,4]=T_cor_p1[,5]
spear_group[8:14,4]=T_cor_p2[,5]
spear_group[15:21,4]=T_cor[,5]

spear_group=data.frame(spear_group)
colnames(spear_group)=c("rho","method","dataset","P")
spear_group[,3]=factor(spear_group[,3],levels=list2)
spear_group[,1]=as.numeric(as.character(spear_group[,1]))
spear_group[,4]=as.numeric(as.character(spear_group[,4]))
spear_group[,1][spear_group[,4]>0.05]=-0.003
spear_group[,1][spear_group[,1]<0]=-0.003

write.csv(spear_group,file="/Users/shansun/Google\ Drive/picrust2/inference_rho.csv")

spear_group=read.csv(file="/Users/shansun/Google\ Drive/picrust2/inference_rho.csv",row.names=1)
spear_group[,3]=factor(spear_group[,3],levels=c("Human_KW", "Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN"))

pdf("/Users/shansun/Google\ Drive/picrust2/bar_rho.pdf", width = 10, height = 5)
ggplot(spear_group, aes(dataset, rho),ggtheme = theme_bw()) +  geom_bar(aes(fill = method), position = "dodge", stat="identity")+
  theme(strip.text.x = element_text(size = 20, angle = 0),
        legend.position="right",
        axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
dev.off()





