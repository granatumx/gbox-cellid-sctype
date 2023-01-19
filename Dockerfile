## Emacs, make this -*- mode: sh; -*-
 
# [tkwolf 20180712] Using specific dated version
FROM granatumx/gbox-r-sdk:1.0.0

RUN R -e 'install.packages("HGNChelper")'
RUN R -e 'install.packages("Seurat")'
RUN R -e 'install.packages("openxlsx")'

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
