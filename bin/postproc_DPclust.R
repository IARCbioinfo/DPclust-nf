# postprocessing of DPclust results
DPclustsres.files = list.files(".",pattern = "bestConsensusResults.RData",recursive = T,all.files = F,full.names = T)

subclrect <- data.frame (xmin=0, xmax=1, ymin=0, ymax=1)

dataset_fracs_list=vector("list",length(DPclustsres.files))

for(i in 1:length(DPclustsres.files)){
  load(DPclustsres.files[i])
  # find cluster names
  masterfile = read_tsv(list.files(".",pattern="_DPmasterfile.txt"))
  dataset_fracs = bind_cols(chr=dataset$chromosome[,1] , pos=dataset$position[,1], subclonal.fractions=dataset$subclonal.fraction , best.assignment.likelihoods=clustering$best.assignment.likelihoods, all.assignment.likelihoods=clustering$all.assignment.likelihoods  , Cluster=clustering$best.node.assignments, DP=dataset$WTCount+dataset$mutCount, 
                            allSameCN = sapply(1:nrow(dataset$totalCopyNumber),function(i) all(dataset$totalCopyNumber[i,]==dataset$totalCopyNumber[i,1])&all(dataset$copyNumberAdjustment[i,]==dataset$copyNumberAdjustment[i,1]) ) )
  colnames(dataset_fracs$subclonal.fractions) = masterfile$subsample
  colnames(dataset_fracs$DP) = masterfile$subsample
  colnames(dataset_fracs$all.assignment.likelihoods) = clustering$cluster.locations[,1]
  
  # merge clusters with uncertain clustering
  dataset_fracs$Cluster[dataset_fracs$best.assignment.likelihoods<0.95 | is.na(dataset_fracs$best.assignment.likelihoods)] = "Uncertain"
  clust_prop = table(dataset_fracs$Cluster[dataset_fracs$Cluster!="Uncertain"])/sum(dataset_fracs$Cluster!="Uncertain") # proportion of high-confidence mutations in each cluster
  # remove clusters with less than 2.5% of mutations assigned with high confidence, and recompute likelihoods
  all.assignment.likelihoods.tmp = dataset_fracs$all.assignment.likelihoods[,colnames(dataset_fracs$all.assignment.likelihoods) %in%names(clust_prop[clust_prop>=0.025])]
  dataset_fracs$all.assignment.likelihoods = sweep(all.assignment.likelihoods.tmp,1,rowSums(all.assignment.likelihoods.tmp),"/")
  dataset_fracs$Cluster = apply(dataset_fracs$all.assignment.likelihoods,1,function(x){ res=colnames(dataset_fracs$all.assignment.likelihoods)[which.max(x)];if(length(res)==0){res=NA};return(res)} )
  dataset_fracs$best.assignment.likelihoods = apply(dataset_fracs$all.assignment.likelihoods,1,function(x){ res=max(x);if(length(res)==0){res=NA};return(res)} )
  dataset_fracs$Cluster[dataset_fracs$best.assignment.likelihoods<0.95 | is.na(dataset_fracs$best.assignment.likelihoods)] = "Uncertain"
  
  # set to 1 CCF>1
  dataset_fracs$subclonal.fractions[dataset_fracs$subclonal.fractions>1] = 1
  
  # merge clonal clusters
  clonalclusts = clustering$cluster.locations[,1][which(clustering$cluster.locations[,2]>0.95)]
  ## set cluster numbering to 1 for all clonal clusters
  dataset_fracs$Cluster[dataset_fracs$Cluster%in%clonalclusts] = 1
  
  # rename clusters to have contiguous numbers
  dataset_fracs$Cluster = as.numeric(factor(dataset_fracs$Cluster))
  dataset_fracs$Cluster[dataset_fracs$Cluster==max(dataset_fracs$Cluster)] = "Uncertain"
  
  # find alterations belonging to a clonal cluster
  dataset_fracs$Clonal = "Uncertain"
  if(length(clonalclusts)==0){#when no clonal cluster detected, all alterations are considered subclonal
    dataset_fracs$Clonal = FALSE
  }else{
    if(length(clonalclusts)>1){#when there are multiple clonal clusters, sum their likelihoods
      # merge likelihoods
      all.assignment.likelihoods.tmp = dataset_fracs$all.assignment.likelihoods
      all.assignment.likelihoods.tmp[,clonalclusts[1]] = rowSums(all.assignment.likelihoods.tmp[,clonalclusts])
      dataset_fracs$all.assignment.likelihoods = all.assignment.likelihoods.tmp[,-clonalclusts[-1]]
      
      dataset_fracs$Clonal[dataset_fracs$all.assignment.likelihoods[,clonalclusts[1]]>=0.95 ] = TRUE
      if(ncol(dataset_fracs$all.assignment.likelihoods)>2){# if there more than 1 subclonal cluster, sum their likelihoods
        dataset_fracs$Clonal[rowSums(dataset_fracs$all.assignment.likelihoods[,-clonalclusts[1]])>=0.95 ] = FALSE
      }else{# when there is a single subclonal cluster, use its likelihood
        dataset_fracs$Clonal[dataset_fracs$all.assignment.likelihoods[,-clonalclusts[1]]>=0.95 ] = FALSE
      }
    }else{#when there is a single clonal cluster, use its likelihood
      dataset_fracs$Clonal[dataset_fracs$all.assignment.likelihoods[,clonalclusts]>=0.95 ] = TRUE
      if(ncol(dataset_fracs$all.assignment.likelihoods)>2){# if there more than 1 subclonal cluster, sum their likelihoods
        dataset_fracs$Clonal[rowSums(dataset_fracs$all.assignment.likelihoods[,-clonalclusts])>=0.95 ] = FALSE
      }else{# when there is a single subclonal cluster, use its likelihood
        dataset_fracs$Clonal[dataset_fracs$all.assignment.likelihoods[,-clonalclusts]>=0.95 ] = FALSE
      }
    }
  }
  dataset_fracs_list[[i]] = dataset_fracs
}

## CCF plots 
for(i in 1:length(DPclustsres.files)){
  # find driver alterations
  # plot
  Fig1 <-  ggplot(dataset_fracs_list[[i]], aes(x=subclonal.fractions[,1],fill=Cluster)) + geom_histogram()+
    theme_classic()   + xlab("CCF") + ylab("Number of alterations")
  
  ggsave( filename = "Cluster_CCF.pdf", Fig1, width = 3.5,height = 3.5)
  print(FigS5A)
}

save(dataset_fracs_list,file = "dataset_fracs_list.Rdata")

