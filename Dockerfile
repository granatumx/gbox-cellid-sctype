## Emacs, make this -*- mode: sh; -*-
 
# [tkwolf 20180712] Using specific dated version
FROM granatumx/gbox-r-sdk:1.0.0

RUN R -e 'install.packages("HGNChelper")'
RUN R -e 'install.packages("Seurat")'
RUN R -e 'install.packages("openxlsx")'

COPY . .

RUN wget "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"
RUN wget "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R"
#RUN wget "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
