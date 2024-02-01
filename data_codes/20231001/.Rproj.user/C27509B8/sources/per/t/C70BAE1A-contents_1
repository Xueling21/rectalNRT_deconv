# --------------Differential Analysis----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Perform differential analysis on GSE119409

rm(list = ls())

library(limma)
library(pheatmap)
library(stringr)
library(dplyr)

#Input GSE119409 gene expression profile
inputFile = "./datasets/GSE119409.txt"
out_path = "./result/"

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

data = read.table(inputFile,
                  header = T,
                  sep = "\t",
                  check.names = F)
data = data[!duplicated(data[, 'id']), ]
rownames(data) = data[, 1]
data = data[, -1]
data = data[rowMeans(data) > 0, ]

conlabel = "response"
treatlabel = "non"

conData = data[, grep(conlabel, colnames(data))]                                #response
treatData = data[, grep(treatlabel, colnames(data))]                            #non-response
conNum = ncol(conData)
treatNum = ncol(treatData)
data = cbind(conData, treatData)                                                #Merge

##Difference analysis
Type = c(rep(conlabel, conNum), rep(treatlabel, treatNum))
design <- model.matrix( ~ 0 + factor(Type))
colnames(design) <- c(conlabel, treatlabel)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(response - non, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#Select the top 300 genes in ascending order of p-value.
diffSig = topTable(fit2, adjust.method = 'fdr', number = 200000)
diffSig = arrange(diffSig, -abs(diffSig$logFC))[1:1000, ]
diffSig = arrange(diffSig, diffSig$P.Value)[1:300, ]


#Output the results of difference analysis
diffSigOut = rbind(id = colnames(diffSig), diffSig)
write.table(
  diffSigOut,
  file = paste0(out_path, "GSE119409_300DEG_DA.txt")
  ,
  sep = "\t",
  quote = F,
  col.names = F
)

#Output differential gene expression profile
diffGeneExp = data[row.names(diffSig), ]
diffGeneExpOut = rbind(id = paste0(colnames(diffGeneExp)), diffGeneExp)
write.table(
  diffGeneExpOut,
  file = paste0(out_path, "GSE119409_300DEG.txt")
  ,
  sep = "\t",
  quote = F,
  col.names = F
)


##t.test
X2 = str_split(colnames(diffGeneExp), '[_]', simplify = T)[, 2]
contreatlist = data.frame(cbind(colnames(diffGeneExp), X2))
diffGeneExp = t(diffGeneExp)

#Calculate mean, variance, and t-value
avg <- aggregate(diffGeneExp, by = list(contreatlist$X2), FUN = mean)
var <- aggregate(diffGeneExp, by = list(contreatlist$X2), FUN = var)
tf <-
  t.test(diffGeneExp[, 1] ~ contreatlist$X2,
         alternative = 't',
         var.equal = T)

pf <- data.frame(t = "V1", df = "V2", p.value = "V3")

#Perform t-test on the paired samples of each gene in a loop to obtain t-values, df values, and p-values.
for (i in 1:ncol(diffGeneExp)) {
  tf <-
    t.test(diffGeneExp[, i] ~ contreatlist$X2,
           alternative = 't',
           var.equal = T)
  pf[i, ] = c(tf$statistic, tf$parameter, tf$p.value)
  rownames(pf)[i] <- colnames(diffGeneExp)[i]
}

avgt = t(avg)
colnames(avgt) <- paste("Mean", avgt[1, ], sep = " ")
avgt = avgt[-1, ]
result = cbind(pf, avgt)
vart = t(var)
colnames(vart) <- paste("Var", vart[1, ], sep = " ")
vart = vart[-1, ]
result = cbind(result, vart)

outTab = rbind(id = colnames(result), result)
write.table(
  outTab,
  file = paste(
    "./result/GSE119409_300DEG_ttest.txt",
    sep = "",
    collapse = NULL
  )
  ,
  sep = "\t",
  quote = F,
  col.names = F
)










