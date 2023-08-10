library(dplyr)
library(Seurat)
library(HGNChelper)
library(data.table)
library(jsonlite)
library(openxlsx)

source('./granatum_sdk.R') # Uses GranatumX SDK
source("./gene_sets_prepare.R")
source("./sctype_score_.R")

threshold <- as.numeric(gn_get_arg('threshold'))
tissue <- gn_get_arg('tissue')
jsonAddition <- gn_get_arg('jsonAddition')
assay <- gn_get_import('assay')

in_matrix_with_gids <- as.data.table(assay$matrix)
rownames(in_matrix_with_gids) <- assay$geneIds
colnames(in_matrix_with_gids) <- assay$sampleIds
print("Check table head")
print(in_matrix_with_gids[1:5, 1:5])
cds <- CreateSeuratObject(counts = in_matrix_with_gids, project = "Proj", assay = "RNA", min.cells = 3, min.features = 1, row.names = assay$geneIds)
print("Check rownames")
print(head(rownames(cds$RNA)))
print("Check colnames")
print(head(colnames(cds$RNA)))
cds <- FindVariableFeatures(cds, selection.method = "vst")
cds <- ScaleData(cds, features = rownames(cds))
cds <- RunPCA(cds, features = VariableFeatures(object = cds))
fload <- "./data/ScTypeDB_from_web.xlsx"
if (nchar(jsonAddition) > 0) {
  json_data <- fromJSON(jsonAddition)
  tissue <- json_data$tissueType[1]
  table_data <- read.xlsx(fload, 1)
  totdata <- rbind(json_data, table_data)
  fload <- "./data/ScTypeDB_from_web_2.xlsx"
  write.xlsx(totdata, fload, row.names = FALSE, showNA = FALSE)
} 
print("=============== Check this ===================")
print(cds)
gs_list <- gene_sets_prepare(fload, tissue)
print("Running scoring now")

scores <- t(sctype_score(scRNAseqData=cds[["RNA"]]@scale.data, scaled=TRUE, gs=gs_list$gs_positive, gs2=gs_list$gs_negative))

#scores['Unknown'] <- threshold

print("Score results")
head(scores)

celltype <- {}
if (threshold > -1000) {
  celltype[assay$sampleIds] <- append(colnames(scores), "Unknown")[max.col(cbind(scores, threshold),ties.method="first")
]
} else {
  celltype[assay$sampleIds] <- colnames(scores)[max.col(scores, ties.method="first")]
}
celltype[assay$sampleIds] <- append(colnames(scores), "Unknown")[max.col(cbind(scores, threshold),ties.method="first")]
#print(celltype)
#print(unbox(celltype))
#gn_export_statically(unbox(as.data.frame(celltype)), 'cellTypeAssignment')
#print(t(as.data.frame(celltype)))
#flush(stdout())
#Sys.sleep(10)


gn_export_statically(unbox(as.data.frame(t(as.data.frame(celltype)))), 'cellTypeAssignment')

#gn_add_result(
#  output[['results']][['numbers_of_genes']],
#  'markdown',
#  list(label = 'Number of genes') # the "description" argument is no longer used
#)

#
# FOR TESTING
# write_json(output, "temp", simplifyVector=T)
#

# write_json(output, stdout(), simplifyVector=T)
gn_commit()
