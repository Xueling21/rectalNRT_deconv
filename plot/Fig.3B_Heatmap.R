# --------------Fig.3B----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Draw a bidirectional clustering heatmap for the results of CIBERSORTx 
# and identify the cell type, receptor, and ligand of each gene(Fig.3B)


rm(list = ls())
setwd(".")
library(stringr)
library(pheatmap)

conlabel = "response"
treatlabel = "non"

input_path = "./result/GSE119409_cs_LR/GSE119409_cs_LR_merge/"
out_path = "./result/"
filename = "GSE119409_cs_LR_merge_all.txt"

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

input = paste0(input_path, filename)
ref <- read.table(input,
                  header = T,
                  sep = "\t",
                  check.names = F)
ref[is.na(ref)] = 0

ref$Celltype = paste0(ref$Celltype, "_", ref$id)
ref[ref$R == 1, ]$L = 2
ref = ref[, c(1, 10)]
colnames(ref)[1] = "ID"
colnames(ref)[2] = "LR"
#1 represents ligand, 2 represents receptor, and 0 represents neither
ref[ref$LR == 0, ]$LR = "Neither"
ref[ref$LR == 1, ]$LR = "Ligand"
ref[ref$LR == 2, ]$LR = "Receptor"
rownames(ref) = ref[, 1]
ref = ref[, -1, drop = F]



input_path = "./datasets/"
foldername = "GSE119409_cs"

filenames = list.files(paste0(input_path, foldername), pattern = ".txt$")
celltype = str_split(str_split(filenames, '[.]', simplify = T)[, 1], '[_]', simplify = T)[, 3]

allTab = data.frame()
alltype = data.frame()
for (j in 1:length(filenames)) {
  input = paste(input_path,
                foldername,
                '/',
                filenames[j],
                sep = "",
                collapse = NULL)
  
  data <-
    read.table(
      input,
      header = T,
      sep = "\t",
      check.names = F,
      row.names = 1
    )
  data = na.omit(data)
  data = data[apply(data, 1, function(x)
    sd(x) != 0), ]
  if (nrow(data) > 1) {
    rownames(data) = paste0(celltype[j], "_", rownames(data))
    
    Type2 = c(rep(celltype[j], nrow(data)))
    names(Type2) = rownames(data)
    Type2 = as.data.frame(Type2)
    
    if (nrow(allTab) == 0) {
      allTab = data
      alltype = Type2
    } else{
      allTab = rbind(allTab, data)
      alltype = rbind(alltype, Type2)
    }
  }
}

data = allTab
Type2 = alltype

conData = data[, grep(conlabel, colnames(data))]
treatData = data[, grep(treatlabel, colnames(data))]
conNum = ncol(conData)
treatNum = ncol(treatData)

geneName = as.vector(rownames(data))
dataLength = length(geneName)
Type1 = c(rep("R", conNum), rep("NR", treatNum))
names(Type1) = colnames(data)
Type1 = as.data.frame(Type1)
colnames(Type1)[1] = "R_NR"
colnames(Type2)[1] = "CellType"

Type2 = cbind(ref, Type2)

outpdf = paste(
  out_path,
  str_split(foldername, '[.]', simplify = T)[, 1],
  '_Bicluster.pdf',
  sep = "",
  collapse = NULL
)

#Answer R blue, NR red
#Cell type B blue, E red, Ma pink, My green, S light yellow, T brown
ann_colors = list(
  R_NR = c(R = "#3B4992FF", NR = "#EE0000FF"),
  CellType = c(
    Bcells = "lightskyblue",
    Epithelialcells = "firebrick1",
    Mastcells = "plum1",
    Myeloids = "forestgreen",
    Stromalcells = "yellow",
    Tcells = "sienna"
  ),
  LR = c(
    Neither = "gray88",
    Receptor = "steelblue1",
    Ligand = "orangered"
  )
)

pdf(file = outpdf,
    width = 15,
    height = 15)
p = pheatmap(
  data,
  annotation_col = Type1,
  annotation_row = Type2,
  annotation_colors = ann_colors,
  annotation_names_row = F,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_cols = T,
  cluster_rows = T,
  show_colnames = F,
  show_rownames = F,
  scale = "row",
  fontsize = 20,
  fontsize_row = 20,
  fontsize_col = 20
)

dev.off()




