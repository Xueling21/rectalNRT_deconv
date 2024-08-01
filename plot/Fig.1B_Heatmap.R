# --------------Fig.1B----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Draw a heatmap of differentially expressed genes (Fig.1B)
heatmap_fig1 <- function() {
#rm(list=ls())
library(pheatmap)

conlabel = "response"
treatlabel = "non"  
response=15                                                                     #response
non=41                                                                          #non-response

input="./result/GSE119409_300DEG.txt"  
DA_name ="./result/GSE119409_300DEG_DA.txt"
outpdf="./result/"                                  #Output heatmap

data <- read.table(input,sep = "\t",header=T,check.names=F,row.names=1)
data_DA=read.table(DA_name,sep="\t",header=T,check.names=F,row.names=1)

geneNum=50
data_DA=data_DA[order(as.numeric(as.vector(data_DA$logFC))),]
diffGeneName=as.vector(rownames(data_DA))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep(conlabel,response),rep(treatlabel,non))
names(Type)=colnames(data)
Type=as.data.frame(Type)
colnames(Type)[1]="R_NR"

ann_colors = list(
  R_NR = c(response="#3B4992FF",non="#EE0000FF")
)

pdf(file=paste0(outpdf,"GSE119409_300DEG_heatmap.pdf"),width=15, height=15)
pheatmap(hmExp, 
         annotation=Type, 
         annotation_colors = ann_colors,
         annotation_names_row = F,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         border_color=F,
         scale="row",
         fontsize = 12,
         fontsize_row=12,
         fontsize_col=12)
dev.off()

pdf(file=paste0(outpdf,"GSE119409_300DEG_heatmap_cluster.pdf"),width=15, height=15)
pheatmap(hmExp, 
         annotation=Type, 
         annotation_colors = ann_colors,
         annotation_names_row = F,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T, 
         show_colnames = F,
         border_color=F,
         scale="row",
         fontsize = 12,
         fontsize_row=12,
         fontsize_col=12)
dev.off()
}