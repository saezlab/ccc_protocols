library(remotes)
library(devtools)
library(BiocManager)

# devtools::install_github("renozao/repotools@0f76e52253c08063084074061f8ca2c05e8a4818", quiet = T, upgrade = F)
# library(repotools)
# cellchat
# nmf.url<-'https://cran.r-project.org/src/contrib/Archive/NMF/NMF_0.23.0.tar.gz'
#repotools::install.pkgs(nmf.url, repos = NULL, force = F, quiet = T)
# repotools::install.pkgs('NMF', force = F, quiet = T) # how to specify version?
# remotes::install_github('sqjin/CellChat@c982c7132386454b3611371a29c4593b3c4a01cc', upgrade=F)
# liana
remotes::install_github('saezlab/liana@ab70b34066f68df60e9ed0d0ce72b0d00f871b7e', upgrade=F)# most recent

remotes::install_github("mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c", upgrade=F)

# decoupleR
BiocManager::install("saezlab/decoupleR@c17d635e0720c86f2386c39ad7dea8614df393f1", update = FALSE)