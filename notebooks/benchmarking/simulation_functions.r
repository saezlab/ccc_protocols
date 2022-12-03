suppressPackageStartupMessages({
    library(splatter)
    
    library(scater)
    library(scran)
    library(bluster)
    
    library(stringr)
    library(ggpubr, quietly = T)
    
    library(plyr, quietly = T)
    library(reshape2, quietly = T)
    
    library(StabEco, quietly = T)
    library(igraph, quietly = T)
    
    library(liana, quietly = T)
    library(tibble, quietly = T)
    c2c <- reticulate::import(module = "cell2cell", as="c2c")
    
    library(ggplot2, quietly = T)

})

seed <- 888
set.seed(seed)
n.cores <- 20

env.name<-'ccc_protocols'
data.path<-'/data3/hratch/ccc_protocols/'

# generate a scale-free, undirected, bipartite PPI network
# emulate from c2c_sim (https://github.com/hmbaghdassarian/c2c_sim) based on Simulate.LR_network method
# with same parameters as in first tensor-cell2cell paper (https://github.com/hmbaghdassarian/tc2c_analyses_1/tree/master/notebooks/time_simulation)
generate.lr.ppi<-function(lr.genes, seed = 888){
    alpha<-2
    degrees<-3
    edges<-NULL

    set.seed(seed)
    B = StabEco::BiGraph$new(n1=length(lr.genes)/2, beta=alpha, k=degrees, m=edges, 
                             type = 'bipartite_sf', directed = T, is_adj = F)# simulate
    G = B$get_graph()  # adjacency matrix

    if (!isSymmetric(G)){stop('Not a bipartite network')}

    node_groups = list(ligand = 1:B$n1, receptor = (B$n1+1):(B$n1+B$n2))

    G[lower.tri(G, diag = TRUE)] <- NA # symmetric
    G<-reshape2::melt(G) #adjacency list
    G = G[!is.na(G$value), ] # remove symmetric bidirectionality
    G = G[G$value != 0, ] # remove disconnected nodes

    colnames(G)<-c('Ligand', 'Receptor', 'Interaction.Strength')
    rownames(G)<-NULL
    G<-G[, c('Ligand', 'Receptor')]

    # map to simulated dataset gene names
    # the nodes must be a subset of those that are interacting
    if ((length(setdiff(G$Ligand, node_groups$ligand)) > 0) | (length(setdiff(G$Receptor, node_groups$receptor)) > 0)){
        stop('Something went wrong in ligand & receptor assignment')
    }
    ligand.map<-setNames(lr.genes[1:B$n1], node_groups$ligand)
    receptor.map<-setNames(lr.genes[(B$n1+1):(B$n1+B$n2)], node_groups$receptor)
    lr.map<-c(ligand.map, receptor.map)
    G[['Ligand']]<-ligand.map[as.character(G$Ligand)]
    G[['Receptor']]<-receptor.map[as.character(G$Receptor)]

    if (length(intersect(G$Ligand, G$Receptor)) != 0){stop('Not a bipartite network')}
    
    # format for LIANA
    colnames(G)<-c('source_genesymbol', 'target_genesymbol')
    G<-as_tibble(G)
    
    return(G)
    
}

qc.data<-function(sce){
    # taken from PMID: 34949812
    
    # QC of cells
    sce <- scater::addPerCellQC(sce) # typical QC as in batch correction paper
    discard <- lapply(unique(colData(sce)$Batch), function(batch) {
        in_batch <- colData(sce)$Batch == batch
        scater::quickPerCellQC(colData(sce)[in_batch, ], nmads = 2)$discard
    })
    discard <- unlist(discard)
    colData(sce)$Discard <- discard
    sce <- sce[, !discard]

    # QC of genes
    sce <- scater::addPerFeatureQC(sce)
    is_exprs <- rowData(sce)$detected >= 0.01
    sce <- sce[is_exprs, ]
    
    return(sce)
}

split.by.context<-function(sim, context_lab = 'Batch'){
    sim.bc<-list()
    contexts<-sort(unique(sim[[context_lab]]))
    for (context in contexts){
        bc<-rownames(colData(sim)[(colData(sim)[[context_lab]] == context),])
        sim.bc[[context]]<-sim[, bc]
    }
    return(sim.bc)
}

score.communication.sce<-function(sce, lr.ppi, 
                                  seed = 888, pos = T, n.cores = 15, expr_prop = 0.1, assay_type = 'logcounts'){
    communication.scores<-liana_wrap(sce = sce, 
                                       method = c('natmi', 'sca'), 
                                       idents_col = 'Group', 
                                       assay.type = assay_type,
                                       expr_prop = expr_prop, # liana default
                                       seed = seed,
                                       parallelize = T, 
                                       workers = n.cores, 
                                       permutation.params = list(nperms = 1), # since we don't use p-values
                                       resource = 'custom',
                                       external_resource = lr.ppi
                                      )

    # filter for columns of interest and format
    communication.scores[['natmi']] <- communication.scores$natmi[,c('source', 'target', 'ligand', 'receptor', 'prod_weight')]
    communication.scores[['sca']] <- communication.scores$sca[,c('source', 'target', 'ligand', 'receptor', 'LRscore')]

    colnames(communication.scores$natmi) <- c('source', 'target', 'ligand', 'receptor', 'score')
    colnames(communication.scores$sca) <- c('source', 'target', 'ligand', 'receptor', 'score')
    
    if (pos){
        if (min(communication.scores$natmi$score) < 0){stop('Unexpected negative score')}
        if (min(communication.scores$sca$score) < 0){stop('Unexpected negative score')}
        }
    
    return(communication.scores)
}

# ... into score.communication.sce
score.communication<-function(sim.list, lr.ppi, do.contexts=NULL, ...){
    if (!is.null(do.contexts)){
        sim.list <- sim.list[do.contexts]
    }
    
    #### score communication
    suppressMessages({
        suppressWarnings({
            score.list<-lapply(sim.list, FUN = function(sce) score.communication.sce(sce = sce, lr.ppi = lr.ppi, ...))
            names(score.list)<-names(sim.list)
        })
    })
        
    # separate into the two scoring methods
    natmi.scores<-list()
    sca.scores<-list()
    for (context in names(sim.list)){
        natmi.scores[[context]]<-score.list[[context]]$natmi
        sca.scores[[context]]<-score.list[[context]]$sca
    }
    
    all.scores <- list(natmi = natmi.scores, sca = sca.scores)
    return(all.scores)
}