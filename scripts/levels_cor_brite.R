b=read.table(file="/Users/Shansun/Google\ Drive/picrust2/kegg_cat.txt",sep="\t")

cat2=names(sort(table(b[,2]),decreasing = T))[c(1,2,4)]
cat3=as.character(unique(b[b[,2]%in%cat2,3]))

dataset=c("china","ggs","gorilla","mouse","chicken","rhizo","defor")
list2=c("Human_KW","Human_TY","Gorilla","Mouse","Chicken","Soil_LWM","Soil_AAN")

rsquared=matrix(nrow=length(cat3),ncol=21)
pvalue=matrix(nrow=length(cat3),ncol=21)
perc_FN=matrix(nrow=length(cat3),ncol=21)
perc_FP=matrix(nrow=length(cat3),ncol=21)
cat5len=matrix(nrow=length(cat3),ncol=21)
spearman_r=matrix(nrow=length(cat3),ncol=21)
spearman_p=matrix(nrow=length(cat3),ncol=21)
i=1
for (n1 in c("_p1","_p2","4fun")){
  for (n2 in dataset){
    filename=paste("/Users/Shansun/Google\ Drive/picrust2/wilcx/",n2,"/",n2,n1,"_K.csv",sep="")
    P_tab=read.csv(file=filename,row.names=1)
    j=0
    for (m in cat3){
      j=j+1
      cat5=as.character(b[!is.na(match(as.character(b[,3]),m)),4])
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
}

rownames(spearman_p)=cat3
rownames(spearman_r)=cat3
rownames(perc_FP)=cat3
rownames(perc_FN)=cat3
rownames(cat5len)=cat3

level1=b[match(cat3,b[,3]),2]
cat5len=cat5len[order(level1),]
perc_FN=perc_FN[order(level1),]
perc_FP=perc_FP[order(level1),]

write.csv(cat5len,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_cat5len.csv")
write.csv(perc_FN,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_perc_FN.csv")
write.csv(perc_FP,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_perc_FP.csv")
write.csv(rsquared,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_rsquared.csv")
write.csv(pvalue,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_pvalue.csv")
write.csv(spearman_r,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_rho.csv")
write.csv(spearman_p,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_spearman_p.csv")

level1=b[match(cat2,b[,3]),2]
cat5len=cat5len[order(level1),]
write.csv(cat5len,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_cat5len.csv")


spearman_sig=spearman_r
spearman_sig[is.na(spearman_sig)]=0
spearman_sig[spearman_p>=0.05]=0
spearman_sig[spearman_sig<0]=0
spearman_sig=spearman_sig[order(level1),]
write.csv(spearman_sig,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_spearman_sig.csv")


level1=b[match(cat2,b[,2]),1]
spearman_sig=spearman_sig[order(level1),]
cat3=cat2[order(level1)]
level1_s=sort(b[match(cat3,b[,2]),2])

col_bar=c("pink","coral","plum2","yellow","lightgreen","lightblue","lightgray","lightyellow")[factor(level1)]

legend_names=levels(factor(level1))
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/category_legend.pdf",height=3,width=20)
plot.new()
legend("center",legend = legend_names, fill=c("pink","coral","plum2","yellow","lightgreen","lightblue","lightgray","lightyellow"), cex=1, horiz = TRUE,bty = "n",text.width=c(0.05,0.12,0.11,-0.0,0.16,0.115,0.1))
dev.off()

#picrust1
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p1_category_rho_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(spearman_sig[,1]),names.arg=rev(rownames(spearman_sig)),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="p1")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p1_category_rho_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(spearman_sig[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="p1")
}
dev.off()

#picrust2
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p2_category_rho_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(spearman_sig[,8]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="p2")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p2_category_rho_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 8:14){
  barplot(rev(spearman_sig[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="p2")
}
dev.off()

#tax4fun
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/tax4fun_category_rho_bar1.pdf",height=10,width=7)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(spearman_sig[,15]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="tax4fun")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/tax4fun_category_rho_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 15:21){
  barplot(rev(spearman_sig[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="tax4fun")
}
dev.off()

#False positives and false negatives
level1=b[match(cat2,b[,3]),2]
perc_FP=perc_FP[order(level1),]
perc_FN=perc_FN[order(level1),]
cat3=cat2[order(level1)]
level1_s=sort(b[match(cat2,b[,3]),2])

#picrust1
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p1_category_FP_bar1.pdf",height=10,width=6)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FP[,1]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="p1")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p1_category_FP_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(perc_FP[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="p1")
}
dev.off()

#picrust2
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p2_category_FP_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FP[,8]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="p2")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p2_category_FP_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 8:14){
  barplot(rev(perc_FP[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="p2")
}
dev.off()

#tax4fun
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/tax4fun_category_FP_bar1.pdf",height=10,width=7)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FP[,15]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="tax4fun")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/tax4fun_category_FP_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 15:21){
  barplot(rev(perc_FP[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="tax4fun")
}
dev.off()

#picrust1
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p1_category_FN_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FN[,1]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="p1")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p1_category_FN_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 1:7){
  barplot(rev(perc_FN[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="p1")
}
dev.off()

#picrust2
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p2_category_FN_bar1.pdf",height=10,width=8)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FN[,8]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="p2")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/p2_category_FN_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 8:14){
  barplot(rev(perc_FN[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="p2")
}
dev.off()

#tax4fun
pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/tax4fun_category_FN_bar1.pdf",height=10,width=7)
par(mar=c(5,25,5,1),mfrow=c(1,1),las=1)
barplot(rev(perc_FN[,15]),names.arg=rev(cat3),las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.lab=2,cex.axis=1,main="tax4fun")
dev.off()

pdf("/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/tax4fun_category_FN_bar2.pdf",height=10,width=18)
par(mar=c(5,1,5,1),mfrow=c(1,7),las=1)
for (i in 15:21){
  barplot(rev(perc_FN[,i]),names.arg=NA,las=1,horiz=T,xlim=c(0,1.1),col=rev(col_bar),cex.axis=1.5,main="tax4fun")
}
dev.off()


#the average T scores of each category
t_wgs=matrix(nrow=length(cat2),ncol=21)
t_predict=matrix(nrow=length(cat2),ncol=21)
i=1
for (n1 in c("_p1","_p2","4fun")){
  for (n2 in dataset){
    filename=paste("/Users/Shansun/Google\ Drive/picrust2/",n2,"/",n2,n1,"_T.csv",sep="")
    P_tab=read.csv(file=filename,row.names=1)
    j=0
    for (m in cat2){
      j=j+1
      cat5=b[!is.na(match(as.character(b[,3]),m)),5]
      pos=na.omit(match(cat5,rownames(P_tab)))
      if (all(is.na(pos))) next
      if (length(pos)==1) next
      adjP6=P_tab[pos,]
      adjP7=as.matrix(adjP6[!is.na(adjP6[,1]) & !is.na(adjP6[,2]),])
      if (is.null(dim(adjP7))) next
      if (dim(adjP7)[1]<3) next
      t_wgs[j,i]=mean(na.omit(adjP6[,1]))
      t_predict[j,i]=mean(na.omit(adjP6[,2]))
    }
    i=i+1
  }
}

avg_tscores=cbind(t_wgs[,1:7],t_predict)
rownames(avg_tscores)=cat2
level1=b[match(cat2,b[,3]),2]
avg_tscores=avg_tscores[order(level1),]

write.csv(avg_tscores,file="/Users/Shansun/Google\ Drive/picrust2/wilcx/cats_new/avg_tscores.csv")
