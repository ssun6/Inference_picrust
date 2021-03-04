library(picante)
library(ggplot2)
library(ggpubr)
library(reshape2)

#permute across samples in each gene and look at their correlations
#setwd("/Users/shansun/Google Drive/picrust2/")
setwd("/users/ssun5/rho_ram/")

list1=c("china","ggs","gorilla","mouse","chicken","rhizo","defor")
list2=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")

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

test_m=c("_p1","_p2","4fun")

T_cor_perm=matrix(ncol=42,nrow=100)
K_cor_perm=matrix(ncol=42,nrow=100)
for (n1 in c(1:3)){
  print (n1)
  m1=test_m[n1]
  for ( m1 in c(1:7)){
    print (m1)
    n=list1[m1]
    n_P2_p=read.csv(file=paste0(n,"/",n,"4fun_n.csv"),header=T,row.names=1)
    rownames(n_P2_p)=sapply(strsplit(rownames(n_P2_p),"\\.\\."),"[[",1)
    n_K1_p=read.table(file=paste0(n,"/",n,"_wgs_p.txt"),sep="\t")
    
    cnames=sort(intersect(colnames(n_P2_p),colnames(n_K1_p)))
    knames=union(rownames(n_P2_p),rownames(n_K1_p))
    n_K1_p=n_K1_p[,cnames]
    n_P1_p=n_P2_p[,cnames]
    n_P1_p=n_P1_p[rowSums(n_P1_p)!=0,]
    n_K1_p=n_K1_p[rowSums(n_K1_p)!=0,]
    
    
    for (j in 1:100){
      print (j)
      n_K1_p=randomizeMatrix(n_K1_p,null.model = "richness",iterations = 1000)
      
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
      
      n_T1=n_T[!is.na(n_T[,1])&!is.na(n_T[,2]),]
      T_cor_perm[j,(n1*7-7+m1)*2-1]=cor.test(n_T1[,1],n_T1[,2],method='spearman')$estimate
      T_cor_perm[j,(n1*7-7+m1)*2]=cor.test(n_T1[,1],n_T1[,2],method='spearman')$p.value
      
      n_K1=n_K[!is.na(n_K[,1])&!is.na(n_K[,2]),]
      K_cor_perm[j,(n1*7-7+m1)*2-1]=cor.test(n_K1[,1],n_K1[,2],method='spearman')$estimate
      K_cor_perm[j,(n1*7-7+m1)*2]=cor.test(n_K1[,1],n_K1[,2],method='spearman')$p.value
    }
    head(T_cor_perm)
  }
}

write.csv(T_cor_perm,file="T_cor_perm.csv")
write.csv(K_cor_perm,file="wilcx_cor_perm.csv")

K_cor_perm=read.csv(file="wilcx/wilcx_cor_perm.csv",row.names=1)
K_cor_perm_t=K_cor_perm[,seq(1,42,2)]
K_cor_perm_t=data.frame(cbind(rep(list2,3),c(rep("PICRUSt",7),rep("PICRUSt2",7),rep("Tax4Fun",7)),t(K_cor_perm_t)))
colnames(K_cor_perm_t)=c("Dataset","Method",c(1:100))
K_cor_perm_m=melt(K_cor_perm_t,id.vars=c(1,2))
K_cor_perm_m[,3]=rep("permuted",2100)
colnames(K_cor_perm_m)[3:4]=c("Type","rho")
K_cor_perm_m[,4]=as.numeric(as.character(K_cor_perm_m[,4]))
K_cor_perm_m[,1]=factor(K_cor_perm_m[,1],levels=list2)

spear_group=read.csv(file="wilcx/wilcx_inference_rho.csv",row.names=1)

spear_group[,3]=factor(spear_group[,3],levels=list2)
spear_group[,1]=as.numeric(as.character(spear_group[,1]))
spear_group[,4]=as.numeric(as.character(spear_group[,4]))
#spear_group[,1][which(spear_group[,4]>0.05)]=-0.003
#spear_group[,1][which(spear_group[,1]<0)]=-0.003
spear_group[,3]=factor(spear_group[,3],levels=c("Human_KW", "Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN"))
colnames(spear_group)=c("rho","Method","Dataset","P")


pdf("wilcx/perm_wilcx_spearman_test_p1.pdf",height=3,width=10)
plot1=ggboxplot(K_cor_perm_m[1:700,], x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.text = element_text(size = 20),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=16)))+
  labs(x = "Dataset", y = "Inference rho")
plot2=plot1+geom_point(data=spear_group[1:7,c(3,2,1)],col="red",size=2)
ggpar(plot2, ylim = c(-0.1,0.6))
dev.off()


pdf("wilcx/perm_wilcx_spearman_test_p2.pdf",height=3,width=10)
plot1=ggboxplot(K_cor_perm_m[701:1400,], x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.text = element_text(size = 20),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=16)))+
  labs(x = "Dataset", y = "Inference rho")
plot2=plot1+geom_point(data=spear_group[8:14,c(3,2,1)],col="red",size=2)
ggpar(plot2, ylim = c(-0.1,0.6))
dev.off()


pdf("wilcx/perm_wilcx_spearman_test4fun.pdf",height=3,width=10)
plot1=ggboxplot(K_cor_perm_m[1401:2100,], x = "Dataset", y = "rho", color = "blue", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),strip.text = element_text(size = 20),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),  
                                           axis.title=element_text(size=20),legend.text=element_text(size=15),legend.title=element_text(size=16)))+
  labs(x = "Dataset", y = "Inference rho")
plot2=plot1+geom_point(data=spear_group[15:21,c(3,2,1)],col="red",size=2)
ggpar(plot2, ylim = c(-0.1,0.6))
dev.off()

