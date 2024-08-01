# --------------Analyze by CellChat----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Analyzing Single Cell Data of Rectal Cancer (GSE132465) Using CellChat

rm(list=ls())

library(CellChat)
library(patchwork)
library(stringr)
options(stringsAsFactors = FALSE)
outpath="./result/cellchat/"

##Part 1. Data preprocessing
load("./datasets/cellchat/GSE132465_tumor_original.RData")#TPM, no cell annotation done

data.input = data
colnames_data=str_split(colnames(data.input),'[.]',simplify=T)
colnames(data.input)=paste0(colnames_data[,2],'-',colnames_data[,3])
rownames(data.input)=str_split(rownames(data.input),'["]',simplify=T)[,2]
meta = label # a dataframe with rownames containing cell mata data
meta$Cell_type = gsub(pattern=' ', replacement='', meta$Cell_type, ignore.case = FALSE, perl = FALSE,
                      fixed = FALSE, useBytes = FALSE)
table(colnames(data.input)%in%meta$Index)
meta=meta[meta$Index%in%colnames(data.input),]
rownames(meta)=meta[,1]
meta=meta[,-1]
data.input = as.matrix(data.input[,rownames(meta)])

#Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_type")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

##Part 2: Inference of Intercellular Communication Networks

#Calculate communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
#If there are only a few cells in certain cell groups, filter out intercellular communication
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

outTab=df.net
write.table(outTab, file=paste0(outpath,"cellchat_dfnet(GSE132465_T).txt"), sep="\t", quote=F, row.names = F,col.names=T)

#Inferring intercellular communication at the level of signaling pathways
computeAveExpr(cellchat, type =  "truncatedMean", trim = 0.1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP$pathways
pathway=as.data.frame(cellchat@netP$pathways)
cellchat_sp=as.data.frame(cellchat@netP[["prob"]])
outTab=rbind(id=colnames(cellchat_sp),cellchat_sp)
write.table(outTab, file=paste0(outpath,"/3.cellchat_sp.txt"), sep="\t", quote=F, row.names = T,col.names=F)

#Inferring intercellular communication at the level of a single receptor pair
cellchat_LR=as.data.frame(cellchat@net[["prob"]])
cellchat_LR=cellchat_LR[,colSums(cellchat_LR)>0]
outTab=rbind(id=colnames(cellchat_LR),cellchat_LR)
write.table(outTab, file=paste0(outpath,"/2.cellchat_LR.txt"), sep="\t", quote=F, row.names = T,col.names=F)

#Aggregated cell cell communication networks can be obtained by calculating the number of links or summarizing communication probabilities.
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#Check the communication strength between a single cell and other cells
par(mfrow = c(1,1), xpd=TRUE)
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[1, ] <- mat[1, ]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[1])

pdf(paste0(outpath,"/3.All_Circle.pdf"),height = 12,width = 12)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

##Part 3: Visualization of Intercellular Communication Networks
pathways.show <- cellchat@netP$pathways
vertex.receiver = seq(1,3) 
netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    signaling.name="All",
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy",                   #层次图
                    vertex.weight=NULL,                     #设置绘制不同大小的顶点
                    pt.title=18,
                    title.space	=3,
                    show.legend=FALSE
)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,signaling.name="All",layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,signaling.name="All", layout = "chord")
#Heatmap
# for(i in 1:length(pathways.show)){
#   pdf(paste0(outpath,i,".",pathways.show[i],"_Heatmap.pdf"),height = 5,width = 7)
#   plot(netVisual_heatmap(cellchat, signaling = pathways.show[i], color.heatmap = "Reds"))
#   dev.off() 
# }
pdf(paste0(outpath,"Heatmap_all.pdf"),height = 5,width = 7)
for(i in 1:length(pathways.show)){
  plot(netVisual_heatmap(cellchat, signaling = pathways.show[i], color.heatmap = "Reds"))
}
dev.off()

#Visualize intercellular communication mediated by multiple ligand receptors or signaling pathways
celltype=names(table(cellchat@idents))
i=6
pdf(paste0(outpath,i,".",celltype[i],".pdf"),height =6,width = 7)
plot(netVisual_bubble(cellchat, sources.use = i, remove.isolate = FALSE))
dev.off() 
pdf(paste0(outpath,"Point_all.pdf"),height = 15,width = 8)
plot(netVisual_bubble(cellchat, sources.use = 1, remove.isolate = FALSE))
plot(netVisual_bubble(cellchat, sources.use = 2, remove.isolate = FALSE))
plot(netVisual_bubble(cellchat, sources.use = 4, remove.isolate = FALSE))
plot(netVisual_bubble(cellchat, sources.use = 5, remove.isolate = FALSE))
plot(netVisual_bubble(cellchat, sources.use = 6, remove.isolate = FALSE))
dev.off()


# Using violin/dot plot to plot the distribution of signal gene expression
# for(i in 1:length(pathways.show)){
#   pdf(paste0(outpath,i,".",pathways.show[i],"_GeneExp.pdf"),height = 9,width = 6)
#   plot(plotGeneExpression(cellchat, signaling = pathways.show[i]))
#   dev.off()
# }
pdf(paste0(outpath,"GeneExp_all.pdf"),height = 9,width = 6)
for(i in 1:length(pathways.show)){
  plot(plotGeneExpression(cellchat, signaling = pathways.show[i]))
}
dev.off()


##Part 4: Analysis of Cellular Communication Network Systems

##Calculate network centrality score
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# Using heat map visualization to calculate centrality scores, it is easy to identify the main signaling roles of cell populations
# for(i in 1:length(pathways.show)){
#   pdf(paste0(outpath,"7.Role/",i,".",pathways.show[i],"_Role.pdf"),height = 5,width = 8)
#   netAnalysis_signalingRole_network(cellchat, signaling = pathways.show[i],height = 7,width = 13, font.size = 15,font.size.title = 15)
#   dev.off()
# }
pdf(paste0(outpath,"Role_all.pdf"),height = 5,width = 8)
for(i in 1:length(pathways.show)){
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show[i],height = 7,width = 13, font.size = 15,font.size.title = 15)
}
dev.off()

# Visualize the main sender (source) and receiver (target) in two-dimensional space
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
gg1 + gg2

# Identify the signals that contribute the most to the outgoing or incoming signals of certain cell populations 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht3 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "all")
ht1 + ht2

# Identifying and visualizing the efferent communication patterns of secretory cells
library(NMF)
selectK(cellchat, pattern = "outgoing")

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

library(ggalluvial)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")


# save(cellchat,pathways.show,file = "./datasets/GSE132465_6cells_cellchat.RData")
# load("./datasets/GSE132465_6cells_cellchat.RData")


