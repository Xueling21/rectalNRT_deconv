# --------------Fig.4----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Draw heatmap of the correlation coefficient between CS(Cell States) and csDEG(Fig.4)

rm(list=ls())
library(stringr)  
library(pheatmap)

setwd(".")
filepath=paste0("./result/GSE119409_cs_ecotyper/")
input_group="./datasets/ecotyper/"

if(!dir.exists(filepath)){
  dir.create(filepath) 
}

input_files = list.files(path = filepath,pattern = paste0(".txt$"))
input_file_group = "GSE119409_cs_minus_updown.txt"

for(i in 1:length(input_files)){
celltype = str_split(str_split(input_files[i],'[.]',simplify=T)[,1],'[_]',simplify=T)[,1]
outpdf = paste0(str_split(input_files[i],'.txt',simplify=T)[,1],".pdf")                                                            #输出的文件名
data <- read.table(paste0(filepath,input_files[i]),sep = "\t",header=T,check.names=F,row.names=1)    
data = data[,!grepl("p.value",colnames(data))]

group <- read.table(paste0(input_group,input_file_group),sep = "\t",header=T,check.names=F,row.names = 1)
colnames(group)=" Type  "
group = group[grepl(celltype,rownames(group)),,drop=FALSE]
data=data[rownames(data)%in%rownames(group),]

datamerge = data[sapply(rownames(group), function(e) {which(rownames(data) == e)}), ]
datamerge = as.matrix(datamerge) 

pdf(paste0(filepath,outpdf),height = 8,width = 8)
pheatmap(datamerge,
         annotation_row  = group, 
         color=colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols = T,
         cluster_rows = T, 
         fontsize = 15,
         fontsize_row = 10,
         fontsize_col = 15,
         border_color = NA
)
dev.off()      

}










