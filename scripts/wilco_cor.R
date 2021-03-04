library(VennDiagram)
library(ggplot2)

#inference calculation
setwd("/Users/shansun/Google Drive/picrust2/")
setwd("/users/ssun5/rho_ram/")

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

list1=c("china","ggs","gorilla","mouse","chicken","rhizo","defor")
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
          Kpv2=-log10(wilcox.test(as.numeric(as.character(n_K1_p[m,group1[[n]]])),as.numeric(as.character(n_K1_p[m,group2[[n]]])))$p.value)
          Kstat2=t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$statistic
        }else{
          Kpv2=log10(wilcox.test(as.numeric(as.character(n_K1_p[m,group1[[n]]])),as.numeric(as.character(n_K1_p[m,group2[[n]]])))$p.value)
          Kstat2=t.test(n_K1_p[m,group1[[n]]],n_K1_p[m,group2[[n]]])$statistic
        }
      }else{
        Kpv2=NA
        Kstat2=NA
      }
      if (m %in% rownames(n_P1_p)){
        if (all(t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$statistic>0)){
          Ppv2=-log10(wilcox.test(as.numeric(as.character(n_P1_p[m,group1[[n]]])),as.numeric(as.character(n_P1_p[m,group2[[n]]])))$p.value)
          Pstat2=t.test(n_P1_p[m,group1[[n]]],n_P1_p[m,group2[[n]]])$statistic
        }else{
          Ppv2=log10(wilcox.test(as.numeric(as.character(n_P1_p[m,group1[[n]]])),as.numeric(as.character(n_P1_p[m,group2[[n]]])))$p.value)
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
    write.csv(n_K,file=paste0("wilcx/",n,"/",n,n1,"_K.csv"))
    
    rownames(n_T)=knames1
    colnames(n_T)=c("metagenome","prediction")
    write.csv(n_T,file=paste0("wilcx/",n,"/",n,n1,"_T.csv"))
    
    n_K1=n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),]#2237
    n_K2=n_K[is.na(n_K[,1])|is.na(n_K[,2]),]#7939
    pdf(paste0("wilcx/",n,"/",n,n1,"_K.pdf"), width = 8, height = 8)
    par(mfrow=c(1,1),mar=c(5,5,3,3))
    plot(n_K1[,1],n_K1[,2],xlab="metagenome",ylab="prediction",pch=16,col="blue",main=list2[j],cex.lab=2,cex.main=2)
    abline(h=0,v=0,cex=1)
    dev.off()
    
    n_T1=n_T[!is.na(n_T[,1])&!is.na(n_T[,2]),]#2237
    n_T2=n_T[is.na(n_T[,1])|is.na(n_T[,2]),]#7939
    pdf(paste0("wilcx/",n,"/",n,n1,"_T.pdf"), width = 8, height = 8)
    par(mfrow=c(1,1),mar=c(5,5,3,3))
    plot(n_T1[,1],n_T1[,2],xlab="metagenome",ylab="prediction",pch=16,col="blue",main=list2[j],cex.lab=2,cex.main=2)
    abline(h=0,v=0,cex=1)
    dev.off()
    
    a=dim(n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),])[1]
    b=dim(n_K[!is.na(n_K[,1]),])[1]
    c=dim(n_K[!is.na(n_K[,2]),])[1] 
    pdf(paste0("wilcx/",n,"/",n,n1,"_venn.pdf"), width = 8, height = 8)
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

write.csv(K_cor,file="wilcx/wilcx_rho.csv")
write.csv(T_cor,file="wilcx/ttest_rho.csv")

K_cor=read.csv(file="/Users/shansun/Google Drive/picrust2/wilcx/wilcx_rho.csv",row.names=1)


inter_num=data.frame(cbind(c(rep("PICRUSt",7),rep("PICRUSt2",7),rep("Tax4Fun",7)),c(rep(list2,3)),1-K_cor[,9]/K_cor[,7]))
colnames(inter_num)=c("method","dataset","Percentage")
inter_num[,3]=as.numeric(as.character(inter_num[,3]))
inter_num[,2]=factor(inter_num[,2],levels=list2)
pdf("wilcx/bar_false_neg.pdf", width = 10, height = 5)
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

inter_num=data.frame(cbind(c(rep("PICRUSt",7),rep("PICRUSt2",7),rep("Tax4Fun",7)),c(rep(list2,3)),1-K_cor[,9]/K_cor[,8]))
colnames(inter_num)=c("method","dataset","Percentage")
inter_num[,3]=as.numeric(as.character(inter_num[,3]))
inter_num[,2]=factor(inter_num[,2],levels=list2)
pdf("wilcx/bar_false_pos.pdf", width = 10, height = 5)
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
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  labs(x = "Dataset", y = "rho of inference P-values")
dev.off()


spear_group=matrix(nrow=21,ncol=4)
spear_group[,1]=K_cor[,2]
spear_group[1:7,2]=rep("PICRUSt",7)
spear_group[8:14,2]=rep("PICRUSt2",7)
spear_group[15:21,2]=rep("Tax4Fun",7)
spear_group[,3]=rep(list2,3)
spear_group[,4]=K_cor[,5]

spear_group=data.frame(spear_group)
colnames(spear_group)=c("rho","method","dataset","P")
write.csv(spear_group,file="wilcx/wilcx_inference_rho.csv")

spear_group=read.csv(file="wilcx/wilcx_inference_rho.csv",row.names=1)

spear_group[,3]=factor(spear_group[,3],levels=list2)
spear_group[,1]=as.numeric(as.character(spear_group[,1]))
spear_group[,4]=as.numeric(as.character(spear_group[,4]))
#spear_group[,1][which(spear_group[,4]>0.05)]=-0.003
#spear_group[,1][which(spear_group[,1]<0)]=-0.003
spear_group[,3]=factor(spear_group[,3],levels=c("Human_KW", "Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN"))

pdf("wilcx/wilcx_bar_rho.pdf", width = 10, height = 5)
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
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  labs(x = "Dataset", y = "rho of inference P-values")
dev.off()





