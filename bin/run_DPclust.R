library(DPClust)

masterfiles = list.files("./","DPmasterfile.txt",full.names = T)

dir.create("DPclust")

for(i in 1:length(masterfiles)){
  opt = list(run_sample= 1,
             data_path="./",
             outputdir="DPclust/",
             input=masterfiles[i],
             keep_temp_files=TRUE,
             analysis_type="nd_dp", 
             iterations=5000, burnin=1000,
             mut_assignment_type=1,
             min_muts_cluster=-1,min_frac_muts_cluster=0.02,num_muts_sample=50000,bin_size=NULL,seed=123, 
             assign_sampled_muts=TRUE )
  
  source("dpclust_pipeline.R")
}
