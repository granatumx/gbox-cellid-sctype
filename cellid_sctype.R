library(dplyr)
library(Seurat)
library(HGCNhelper)

source('./granatum_sdk.R') # Uses GranatumX SDK
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

tissue <- gn_get_arg('tissue')
assay <- gn_get_import('assay')

in_matrix_with_gids <- assay$matrix
rownames(in_matrix_with_gids) <- assay$geneIds
colnames(in_matrix_with_gids) <- assay$sampleIds
cds <- newCellDataSet(in_matrix_with_gids, lowerDetectionLimit = 0, expressionFamily = negbinomial.size())
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", tissue)
print("Running scoring now")

scores <- sctype_score(scRNAseqData=cds[["RNA"]]@scale.data, scaled=TRUE, gs=gs_list$gs_positive, gs2=gs_list$gs_negative)

celltype <- colnames(scores)[max.col(scores,ties.method="first")]
celltype$sampleIds <- assay$sampleIds
gn_export_statically(celltype, 'cellTypeAssignment')

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
