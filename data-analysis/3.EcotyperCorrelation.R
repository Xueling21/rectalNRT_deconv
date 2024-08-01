# --------------Analyze by Ecotyper----------------
# Prepred by Xiao Sun, Xueling Li and Min Zhu
# Contact Xueling Li via email: xlli@cmpt.ac.cn or xuelingli16@foxmail.com
# Analyze the results of Ecotyper and CIBERSORTx

# Calculate the correlation coefficient between CS(Cell States) and csDEG
rm(list = ls())
library(stringr)

setwd(".")

mode_cor = "pearson"
input_filepath = "./datasets/ecotyper/"
out_filepath = paste0("./result/GSE119409_cs_ecotyper/")

if (!dir.exists(out_filepath)) {
  dir.create(out_filepath)
}

input_file = "GSE119409_cs_minus.txt"

input1 = paste0(input_filepath, input_file)
data1 <-
  read.table(
    input1,
    header = T,
    sep = "\t",
    check.names = F,
    row.names = 1
  )
celltype = names(table(str_split(rownames(data1), '[_]', simplify = T)[, 1]))
data1 = as.data.frame(t(data1))

for (num_file_state in 1:length(celltype)) {
  input2 = paste0(
    input_filepath,
    "/Ecotyper_Cell_States/",
    celltype[num_file_state],
    "_Cell_State_Abundance.txt"
  )
  data2 <-
    read.table(
      input2,
      header = T,
      sep = "\t",
      check.names = F,
      row.names = 1
    )
  
  merge_cor = data.frame()
  merge_p = data.frame()
  
  for (i in 1:ncol(data1)) {
    for (j in 1:ncol(data2)) {
      #The method can be "spearman," "pearson," and "kendall," 
      #corresponding to the calculation and testing of three correlation coefficients.
      cor = cor.test(x = data1[, i],
                     y = data2[, j],
                     method = mode_cor)
      merge_cor[i, j] = cor[["estimate"]][["cor"]]         #pearson
      # merge_cor[i,j]=cor[["estimate"]][["rho"]]       #spearman
      # merge_cor[i,j]=cor[["estimate"]][["tau"]]       #kendall
      merge_p[i, j] = cor[["p.value"]]
    }
  }
  
  rownames(merge_cor) = colnames(data1)
  colnames(merge_cor) = colnames(data2)
  colnames(merge_p) = paste0("p.value_", colnames(data2))
  rownames(merge_p) = colnames(data1)
  
  merge_all = data.frame()
  merge_all = as.data.frame(merge_cor[, 1, drop = FALSE])
  for (i in 2:ncol(merge_cor)) {
    merge_all = cbind(merge_all, merge_p[, i - 1, drop = FALSE])
    merge_all = cbind(merge_all, merge_cor[, i, drop = FALSE])
  }
  merge_all = cbind(merge_all, merge_p[, ncol(merge_cor), drop = FALSE])
  
  
  outTab = rbind(id = colnames(merge_all), merge_all)
  write.table(
    outTab,
    file = paste0(out_filepath, celltype[num_file_state], '_cor_', mode_cor, '.txt'),
    sep = "\t",
    quote = F,
    row.names = T,
    col.names = F
  )
}
