# --------------Fig.3A----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Draw a violin chart based on the results of the CIEBRSORTx (Fig.3A)

# Sort the results of CIEBRSORTx by the absolute value of t, and then select the top 12 to draw a violin chart
rm(list = ls())
setwd(".")

library(stringr)
library(ggpubr)
library(dplyr)
conlabel = "response"
treatlabel = "non"

data_filepath = "./datasets/"
out_filepath = "./result/"
input_file = "GSE119409_cs_LR_merge.txt"

convertP <- function(p = 1) {
  if (p < 0.001)
    str = "***"
  else if (p < 0.01)
    str = "**"
  else if (p < 0.05)
    str = "*"
  else if (p > 0.05)
    str = "ns"
  return(str)
}

test_name = paste0(data_filepath, input_file)
data_test = read.table(test_name, sep = "\t", header = T)
data_test = arrange(data_test, -abs(data_test$t))
data_test = data_test[data_test$p.value < 0.05, ]
out_path = out_filepath

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

data_violin = data.frame()

ngene = nrow(data_test)
if (ngene > 30)
  ngen = 30
for (j in 1:ngene) {
  data_name = paste0(data_filepath,
                     "GSE119409_cs/GSE119409_cs_",
                     data_test[j, ]$Celltype,
                     ".txt")
  cellid = data_test[j, ]$Celltype
  geneid = data_test[j, ]$id
  
  datam = read.table(data_name, sep = "\t", header = T)
  rownames(datam) = datam[, 1]
  datam = datam[, -1]
  datam = t(datam)
  datam = as.data.frame(datam[, geneid])
  colnames(datam)[1] = paste0(cellid, "_", geneid)
  
  if (j == 1) {
    data_violin = datam
  } else{
    data_violin = cbind(data_violin, datam)
  }
}
data_violin = round(log2(data_violin + 1), 3)

response = 15
non = 41
outpdf = paste0(out_path, "GSE119409_30csDEG_violin.pdf")

pdf(outpdf, height = 8, width = 13)
par(las = 1, mar = c(8, 5, 3, 3))
x = c(1:ncol(data_violin))
y = c(1:ncol(data_violin))
plot(
  x,
  y,
  xlim = c(0, 3 * ncol(data_violin)),
  ylim = c(min(data_violin) - 0.02, max(data_violin) + 0.02),
  main = paste0("Top DEGs of GSE119409 hrCTD(s-mode)"),
  xlab = "",
  ylab = "Gene Expression",
  cex.main = 1.8,
  cex.lab = 1.8,
  cex.axis = 1.8,
  pch = 21,
  col = "white",
  xaxt = "n"
)
plot(legend(
  "topright",
  fill = c("blue", "red"),
  legend = c("R", "NR"),
  title = "Sensitivity"
))

library(vioplot)
ncex = 3
for (i in 1:ncol(data_violin)) {
  normalData = data_violin[1:response, i]
  tumorData = data_violin[(response + 1):(response + non), i]
  vioplot(
    normalData,
    at = ncex * (i - 1),
    lty = 1,
    add = T,
    col = "blue"
  )
  vioplot(
    tumorData,
    at = ncex * (i - 1) + 1,
    lty = 1,
    add = T,
    col = "red"
  )
  wilcoxTest = wilcox.test(normalData, tumorData)
  p = data_test$p.value[i]
  mx = max(c(normalData, tumorData))
  lines(c(x = ncex * (i - 1) + 0.1, x = ncex * (i - 1) + 1), c(mx, mx))
  text(
    x = ncex * (i - 1) + 0.5,
    y = mx + 0.2,
    labels = convertP(p),
    cex = 1.5
  )
  text(
    seq(1, ncex * ncol(data_violin), ncex),
    -0.1,
    xpd = NA,
    labels = colnames(data_violin),
    cex = 1.1,
    srt = 45,
    pos = 2
  )
}

dev.off()




