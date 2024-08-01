# --------------Fig.1A----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Draw a violin chart based on the results of the difference analysis (Fig.1A)

# install.packages("vioplot")                                                   
rm(list=ls())

library(dplyr)
setwd(".")
input="./result/GSE119409_300DEG.txt"                                           #Input gene expression profile
outpdf="./result/GSE119409_30DEG_violin.pdf"                                    #Output Violin Chart
response=15                                                                     #response
non=41                                                                          #non-response

convertP<-function(p=1){
  if(p<0.001)    str="***"
  else if(p<0.01)    str="**"
  else if(p<0.05)    str="*"
  else if(p>0.05)    str="ns"
  return(str)
}

data <- read.table(input,sep = "\t",header=T,check.names=F,row.names=1)  
data=data[,which(colSums(data) > 0.2)]     

test_name ="./result/GSE119409_300DEG_ttest.txt"
data_test=read.table(test_name,sep="\t",header=T)
#Select the top ranked genes by descending the absolute value of t.
data_test = arrange(data_test,-abs(data_test$t))
data_test = data_test[data_test$p.value<0.05,]

data=data[rownames(data) %in% data_test$id,]
data=data[match(data_test$id,rownames(data)),]
data=data[1:30,]
data=t(data)

pdf(outpdf,height = 8,width = 15)
par(las=1,mar=c(8,5,3,3))                                                       
x=c(1:ncol(data))                                                               
y=c(1:ncol(data))
# 绘制框架
plot(x,y,                                          
     xlim = c(0,3*ncol(data)),ylim = c(min(data)-0.02,max(data)+0.02),
     main="Top 30 DEGs of GSE119409",
     xlab = "",
     ylab = "Gene Expression",
     cex.main=1.8,
     cex.lab=1.8,
     cex.axis=1.8,
     pch=21,
     col="white",
     xaxt="n"
     )
plot(legend("topright", fill = c("blue", "red"), legend = c("R", "NR"), title = "Sensitivity"))


library(vioplot)
#Draw a violin diagram for each immune cell cycle. Resistance is indicated in blue, sensitivity is indicated in red.

ncex=3
for(i in 1:ncol(data)){
  normalData=data[1:response,i]
  tumorData=data[(response+1):(response+non),i]
  vioplot(normalData,at=ncex*(i-1),lty=1,add=T,col="blue")
  vioplot(tumorData,at=ncex*(i-1)+1,lty=1,add=T,col="red")
  p=data_test$p.value[i]
  mx=max(c(normalData,tumorData))
  lines(c(x=ncex*(i-1)+0.1,x=ncex*(i-1)+0.8),c(mx,mx))                          
  text(x=ncex*(i-1)+0.5,y=mx+0.15,labels = convertP(p),cex=1.2)
  text(seq(1,ncex*ncol(data),ncex),1.5,xpd=NA,labels = colnames(data),cex = 1,srt=45,pos = 2) 
}

dev.off()                                                                       







