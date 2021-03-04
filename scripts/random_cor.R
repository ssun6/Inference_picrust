library(picante)
library(ggplot2)
library(ggpubr)
library(reshape2)

setwd("/Users/shansun/Google Drive/picrust2/")
setwd("/users/ssun5/rho_ram/")

list1=c("china","ggs","gorilla","mouse","chicken","rhizo","defor")
list2=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")

rho_test=matrix(c("t","p"),nrow=7,ncol=12)
test_m=c("_p1","_p2","4fun")
for (n1 in c(1:3)){
  rho_cor=matrix(c("rho","dataset","group"),nrow=1,ncol=3)
  m1=test_m[n1]
  
  j=1
  for ( n in list1){
    print(n)
    
    n_P2_p=read.csv(file=paste0(n,"/",n,m1,"_n.csv"),header=T,row.names=1)
    rownames(n_P2_p)=sapply(strsplit(rownames(n_P2_p),"\\.\\."),"[[",1)
    n_K1_p=read.table(file=paste0(n,"/",n,"_wgs_p.txt"),sep="\t")
    
    cnames=sort(intersect(colnames(n_P2_p),colnames(n_K1_p)))
    knames=union(rownames(n_P2_p),rownames(n_K1_p))
    n_K1_p=n_K1_p[,cnames]
    n_P1_p=n_P2_p[,cnames]
    n_P1_p=n_P1_p[rowSums(n_P1_p)!=0,]
    n_K1_p=n_K1_p[rowSums(n_K1_p)!=0,]
    
    comm3=intersect(rownames(n_P1_p),rownames(n_K1_p))
    
    cor_n=vector()
    i=1
    for (m in 1:dim(n_K1_p)[2]){
      cor1=cor.test(n_P1_p[comm3,m],n_K1_p[comm3,m],method="spearman")
      cor_n[i]=cor1$estimate
      i=i+1
    }
    r_cor1=as.matrix(cbind(cor_n,rep(list2[j],length(cor_n)),rep("unpermuted",length(cor_n))))
    rho_cor=rbind(rho_cor,r_cor1)
    
    cor_r=vector()
    i=1
    for (r in 1:10){
      wg_perm=randomizeMatrix(n_K1_p,null.model = "richness",iterations = 1000)
      for (m in 1:dim(n_K1_p)[2]){
        cor1=cor.test(n_P1_p[comm3,m],wg_perm[comm3,m],method="spearman")
        cor_r[i]=cor1$estimate
        i=i+1
      }
    }
    
    r_cor2=as.matrix(cbind(cor_r,rep(list2[j],length(cor_r)),rep("permuted",length(cor_r))))
    rho_cor=rbind(rho_cor,r_cor2)
    
    rho_test[j,n1*4-3]=t.test(cor_n,cor_r)$estimate[1]
    rho_test[j,n1*4-2]=t.test(cor_n,cor_r)$estimate[2]
    rho_test[j,n1*4-1]=t.test(cor_n,cor_r)$statistic
    rho_test[j,n1*4]=t.test(cor_n,cor_r)$p.value    
    j=j+1
  }
  
  colnames(rho_cor)=rho_cor[1,]
  rho_cor=rho_cor[-1,]
  rho_cor=data.frame(rho_cor)
  rho_cor[,1]=as.numeric(as.character(rho_cor[,1]))
  rho_cor[,2]=factor(rho_cor[,2],levels=list2)
  
  write.csv(rho_cor,file=paste0("rho_cor",m1,".csv"))
  
  pdf_name=paste0("rho_random",m1,".pdf")
  pdf(pdf_name, width = 12, height = 5)
  gplot=ggboxplot(rho_cor, x = "dataset", y = "rho", color = "group", palette = c("blue","red"), add = "jitter",  ggtheme = theme_bw()+
              theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                    panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),
                    axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
  print(gplot)
  dev.off()
  
}
rownames(rho_test)=list2
colnames(rho_test)=c("mean_up_p1","mean_p_p1","t_p1","P_p1","mean_up_p2","mean_p_p2","t_p2","P_p2","mean_up4fun","mean_p4fun","t4fun","P4fun")
write.csv(rho_test,file="rho_random_compare.csv")

rho_gene_test=matrix(c("t","p"),nrow=7,ncol=12)
for (n1 in c(1:3)){
  rho_cor=matrix(c("rho","dataset","group","KEGG"),nrow=1,ncol=4)
  m1=test_m[n1]
  
  j=1
  for ( n in list1){
    print(n)
    
    n_P2_p=read.csv(file=paste0(n,"/",n,m1,"_n.csv"),header=T,row.names=1)
    rownames(n_P2_p)=sapply(strsplit(rownames(n_P2_p),"\\.\\."),"[[",1)
    n_K1_p=read.table(file=paste0(n,"/",n,"_wgs_p.txt"),sep="\t")
    
    cnames=sort(intersect(colnames(n_P2_p),colnames(n_K1_p)))
    knames=union(rownames(n_P2_p),rownames(n_K1_p))
    n_K1_p=n_K1_p[,cnames]
    n_P1_p=n_P2_p[,cnames]
    n_P1_p=n_P1_p[rowSums(n_P1_p)!=0,]
    n_K1_p=n_K1_p[rowSums(n_K1_p)!=0,]
    
    abd_name=rownames(n_K1_p)[apply(n_K1_p,1,function(i){sum(i!=0)})>ncol(n_K1_p)*0.9]
    comm3=intersect(rownames(n_P1_p),abd_name)
    
    cor_n=vector()
    p_n=vector()
    i=1
    for (m in comm3){
      cor1=cor.test(as.numeric(n_P1_p[m,]),as.numeric(n_K1_p[m,]),method="spearman")
      cor_n[i]=cor1$estimate
      p_n[i]=cor1$p.value
      i=i+1
    }
    r_cor1=as.matrix(cbind(cor_n,rep(list2[j],length(cor_n)),rep("unpermuted",length(cor_n)),comm3))
    rho_cor=rbind(rho_cor,r_cor1)
    
    cor_r=vector()
    i=1
    for (r in 1:10){
      wg_perm=randomizeMatrix(n_K1_p,null.model = "frequency",iterations = 1000)
      for (m in comm3){
        cor1=cor.test(as.numeric(n_P1_p[m,]),as.numeric(wg_perm[m,]),method="spearman")
        cor_r[i]=cor1$estimate
        i=i+1
      }
    }
    
    r_cor2=as.matrix(cbind(cor_r,rep(list2[j],length(cor_r)),rep("permuted",length(cor_r)),rep(comm3,100)))
    rho_cor=rbind(rho_cor,r_cor2)
    
    rho_gene_test[j,n1*4-3]=t.test(cor_n,cor_r)$estimate[1]
    rho_gene_test[j,n1*4-2]=t.test(cor_n,cor_r)$estimate[2]
    rho_gene_test[j,n1*4-1]=wilcox.test(cor_n,cor_r)$statistic
    rho_gene_test[j,n1*4]=wilcox.test(cor_n,cor_r)$p.value
    j=j+1
  }
  
  colnames(rho_cor)=rho_cor[1,]
  rho_cor=rho_cor[-1,]
  rho_cor=data.frame(rho_cor)
  rho_cor[,1]=as.numeric(as.character(rho_cor[,1]))
  rho_cor[,2]=factor(rho_cor[,2],levels=list2)
  
  write.csv(rho_cor,file=paste0("wilcx/rho_cor_genes_abd",m1,".csv"))
  
  pdf_name=paste0("wilcx/rho_random_genes",m1,".pdf")
  pdf(pdf_name, width = 12, height = 5)
  gplot=ggboxplot(rho_cor, x = "dataset", y = "rho", color = "group", palette = c("blue","red"), add = "jitter",  ggtheme = theme_bw()+
                    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                          panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=15),
                          axis.title=element_text(size=15),legend.text=element_text(size=15),legend.title=element_text(size=16)))
  print(gplot)
  dev.off()
  
}

rownames(rho_gene_test)=list2
colnames(rho_gene_test)=c("mean_up_p1","mean_p_p1","t_p1","P_p1","mean_up_p2","mean_p_p2","t_p2","P_p2","mean_up4fun","mean_p4fun","t4fun","P4fun")
write.csv(rho_gene_test,file="rho_genes_random_compare_abd.csv")
                              
#permute across genes
pdf(file="soil1_across_gene.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(as.numeric(n_P1_p["K01191",]),as.numeric(n_K1_p["K01191",]),xlab="Relative abundance from PICRUSt",ylab="Relative abundance from metagenome",pch=16,col="red",main="Unpermuted\nSpearman's rho = 0.85\n P = 2e-6")
wg_perm=randomizeMatrix(n_K1_p,null.model = "frequency",iterations = 10000000000)
plot(as.numeric(n_P1_p["K01191",]),as.numeric(wg_perm["K01191",]),xlab="Relative abundance from PICRUSt",ylab="Permuted relative abundance from metagenome",pch=16,col="blue",main="Permuted\nSpearman's rho = =0.39\n P =0.04")
dev.off()
cor.test(as.numeric(n_P1_p["K01191",]),as.numeric(n_K1_p["K01191",]),method="spearman")
# "	Spearman's rank correlation rho
# 
# data:  as.numeric(n_P1_p["K01191", ]) and as.numeric(n_K1_p["K01191", ])
# S = 490, p-value = 2e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.8504 "
cor.test(as.numeric(n_P1_p["K01191",]),as.numeric(wg_perm["K01191",]),method="spearman")

# "	Spearman's rank correlation rho
# 
# data:  as.numeric(n_P1_p["K01191", ]) and as.numeric(wg_perm["K01191", ])
# S = 4600, p-value = 0.04
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.3966 
# "

rho_p1=read.csv(file="wilcx/rho_cor_genes_abd_p1.csv",row.names=1)
rho_p2=read.csv(file="wilcx/rho_cor_genes_abd_p2.csv",row.names=1)
rho4fun=read.csv(file="wilcx/rho_cor_genes_abd4fun.csv",row.names=1)

rho_p1=cbind(rho_p1,rep("PICRUSt",nrow(rho_p1)))
rho_p2=cbind(rho_p2,rep("PICRUSt2",nrow(rho_p2)))
rho4fun=cbind(rho4fun,rep("Tax4Fun",nrow(rho4fun)))

colnames(rho_p1)[5]="method"
colnames(rho_p2)[5]="method"
colnames(rho4fun)[5]="method"

rho_perm=rbind(rho_p1,rho_p2,rho4fun)
rho_perm[,2]=factor(rho_perm[,2],levels=list2)

png("wilcx/perm_wilcx_spearman_test.png",height=500,width=500)
plot1=ggboxplot(rho_perm, x = "dataset", y = "rho", color = "group", add = "jitter", show.legend=TRUE,
                ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=5),
                                           panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=5),  
                                           axis.title=element_text(size=5),legend.text=element_text(size=5),legend.key.size = unit(1,"line"), legend.title=element_text(size=5)))
plot2=plot1+facet_wrap(~ method, ncol=1)
ggpar(plot2, ylim = c(-1,1))
dev.off()


ggplot(rho_perm, aes(dataset, rho),ggtheme = theme_bw()) +  geom_bar(aes(fill = method), position = "dodge", stat="identity")+
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
