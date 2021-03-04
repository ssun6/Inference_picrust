library(picante)
library(ggplot2)
library(ggpubr)
library(reshape2)

setwd("/users/ssun5/rho_ram/")
list1=c("china","ggs","gorilla","chicken","mouse","rhizo","defor")
list2=c("Human_KW","Human_TY","Gorilla","Chicken","Mouse","Soil_LWM","Soil_AAN")

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



T_cor=matrix(ncol=21,nrow=100)
K_cor=matrix(ncol=21,nrow=100)
test_m=c("_p1","_p2","4fun")
m2=1
for (m1 in test_m){
  for ( n in list1){
    n_P2_p=read.csv(file=paste0(n,"/",n,m1,"_n.csv"),header=T,row.names=1)
    rownames(n_P2_p)=sapply(strsplit(rownames(n_P2_p),"\\.\\."),"[[",1)
    n_K1_p=read.table(file=paste0(n,"/",n,"_wgs_p.txt"),sep="\t")
    
    cnames=sort(intersect(colnames(n_P2_p),colnames(n_K1_p)))
    n_K1_p=n_K1_p[,cnames]
    n_P1_p=n_P2_p[,cnames]
    n_P1_p=n_P1_p[rowSums(n_P1_p)!=0,]
    n_K1_p=n_K1_p[rowSums(n_K1_p)!=0,]
    
    
    a1=group1[[n]]
    a2=group2[[n]]
    for (j in 1:100){
      b1=a1[sample.int(length(a1), 5, replace = FALSE)]
      b2=a2[sample.int(length(a2), 5, replace = FALSE)]
      
      knames1=union(rownames(n_P1_p),rownames(n_K1_p))
      Kpv=list()
      Ppv=list()
      Kstat=list()
      Pstat=list()
      i=1
      for (m in knames1){
        if (m %in% rownames(n_K1_p)){
          if (sd(n_K1_p[m,b1])!=0|sd(n_K1_p[m,b2])!=0){
            if (!is.na(t.test(n_K1_p[m,b1],n_K1_p[m,b2])$statistic)){
              if (all(t.test(n_K1_p[m,b1],n_K1_p[m,b2])$statistic>0)){
                Kpv2=-log10(wilcox.test(as.numeric(as.character(n_K1_p[m,b1])),as.numeric(as.character(n_K1_p[m,b2])))$p.value)
                Kstat2=t.test(n_K1_p[m,b1],n_K1_p[m,b2])$statistic
              }else{
                Kpv2=log10(wilcox.test(as.numeric(as.character(n_K1_p[m,b1])),as.numeric(as.character(n_K1_p[m,b2])))$p.value)
                Kstat2=t.test(n_K1_p[m,b1],n_K1_p[m,b2])$statistic
              }
            }else{
              Kpv2=NA
              Kstat2=NA
            }
          }else{
            Kpv2=NA
            Kstat2=NA
          }
        }else{
          Kpv2=NA
          Kstat2=NA
        }
        if (m %in% rownames(n_P1_p)){
          if (sd(n_P1_p[m,b1])!=0|sd(n_P1_p[m,b2])!=0){
            if (!is.na(t.test(n_P1_p[m,b1],n_P1_p[m,b2])$statistic)){
              if (all(t.test(n_P1_p[m,b1],n_P1_p[m,b2])$statistic>0)){
                Ppv2=-log10(wilcox.test(as.numeric(as.character(n_P1_p[m,b1])),as.numeric(as.character(n_P1_p[m,b2])))$p.value)
                Pstat2=t.test(n_P1_p[m,b1],n_P1_p[m,b2])$statistic
              }else{
                Ppv2=log10(wilcox.test(as.numeric(as.character(n_P1_p[m,b1])),as.numeric(as.character(n_P1_p[m,b2])))$p.value)
                Pstat2=t.test(n_P1_p[m,b1],n_P1_p[m,b2])$statistic
              }
            }else{
              Ppv2=NA
              Pstat2=NA
            }
          }else{
            Ppv2=NA
            Pstat2=NA
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
      
      n_T=cbind(unlist(Kstat),unlist(Pstat))
      n_T1=n_T[!is.na(n_T[,1])&!is.na(n_T[,2]),]
      T_cor[j,m2]=cor.test(n_T1[,1],n_T1[,2],method='spearman')$estimate
      
      n_K=cbind(unlist(Kpv),unlist(Ppv))
      n_K1=n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),]
      K_cor[j,m2]=cor.test(n_K1[,1],n_K1[,2],method='spearman')$estimate
    }
    print(m2)
    m2=m2+1 
  }
}

write.csv(T_cor,file="ttest_cor_subsample.csv")
write.csv(T_cor,file="wilcx_cor_subsample.csv")


K_cor=read.csv(file="wilcx/wilcx_cor_subsample.csv",row.names=1)
spear_group=read.csv(file="wilcx/wilcx_inference_rho.csv",row.names=1)


K_cor1=K_cor[,1:7]
colnames(K_cor1)=list2
K_cor1_m=melt(K_cor1)
colnames(K_cor1_m)=c("Dataset","rho")
spear_group1=spear_group[1:7,c(1,3)]
colnames(spear_group1)=c("rho","Dataset")
K_cor1_m[,1]=factor(K_cor1_m[,1],levels=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN"))
pdf("wilcx/wilcx_subsample_ttest_p1.pdf",height=6,width=10)
plot1=ggboxplot(K_cor1_m, x = "Dataset", y = "rho", color = "blue", add = "jitter",
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
plot2=plot1+geom_point(data=spear_group1,col="red", show.legend=TRUE)+labs(x = "Dataset", y = "Inference rho")
ggpar(plot2, ylim = c(-0.5,1))
dev.off()

K_cor2=K_cor[,8:14]
colnames(K_cor2)=list2
K_cor2_m=melt(K_cor2)
colnames(K_cor2_m)=c("Dataset","rho")
spear_group2=spear_group[8:14,c(1,3)]
colnames(spear_group2)=c("rho","Dataset")
K_cor2_m[,1]=factor(K_cor2_m[,1],levels=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN"))
pdf("wilcx/wilcx_subsample_ttest_p2.pdf",height=6,width=10)
plot1=ggboxplot(K_cor2_m, x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
plot2=plot1+geom_point(data=spear_group2,col="red", show.legend=TRUE)+labs(x = "Dataset", y = "Inference rho")
ggpar(plot2, ylim = c(-0.5,1))
dev.off()

K_cor3=K_cor[,15:21]
colnames(K_cor3)=list2
K_cor3_m=melt(K_cor3)
colnames(K_cor3_m)=c("Dataset","rho")
spear_group3=spear_group[15:21,c(1,3)]
colnames(spear_group3)=c("rho","Dataset")
K_cor3_m[,1]=factor(K_cor3_m[,1],levels=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN"))
pdf("wilcx/wilcx_subsample_ttest4fun.pdf",height=6,width=10)
plot1=ggboxplot(K_cor3_m, x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
plot2=plot1+geom_point(data=spear_group3,col="red", show.legend=TRUE)+labs(x = "Dataset", y = "Inference rho")
ggpar(plot2, ylim = c(-0.5,1))
dev.off()



