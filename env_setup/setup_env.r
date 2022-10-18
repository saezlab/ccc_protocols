library(remotes)
library(devtools)

devtools::install_github("renozao/repotools@0f76e52253c08063084074061f8ca2c05e8a4818", quiet = T, upgrade = F)
library(repotools)

# nmf.url<-'https://cran.r-project.org/src/contrib/Archive/NMF/NMF_0.23.0.tar.gz'
#repotools::install.pkgs(nmf.url, repos = NULL, force = F, quiet = T)
repotools::install.pkgs('NMF', force = F, quiet = T) # how to specify version?

remotes::install_github('sqjin/CellChat@c982c7132386454b3611371a29c4593b3c4a01cc', upgrade=F)
remotes::install_github('saezlab/liana@f7ec81b18eb63bc7a7f941e8e2195139f6018b9c', upgrade=F)
remotes::install_github('theislab/kBET@f35171dfb04c7951b8a09ac778faf7424c4b6bc0', upgrade = F)