###################################################
# Script to preprocess mutation data for DPclust. #
# Adapted from DPclust3p github repo              #
###################################################

# create useful function
PurpleToBattenberg <- function (samplename, outfile.folder, copynumber.purple.file, copynumber.purple.summary){
  d = read_tsv(copynumber.purple.file ) %>% filter(sample==samplename) %>% 
    mutate(minorAlleleCopyNumber = as.numeric(minorAlleleCopyNumber),majorAlleleCopyNumber = as.numeric(majorAlleleCopyNumber))
  # keep only clonal CN, i.d. with fractional part less than 0.2 or greater than 0.8
  d= d %>% filter(minorAlleleCopyNumber %%1<=0.2 | minorAlleleCopyNumber %%1>=0.8,
                  majorAlleleCopyNumber %%1<=0.2 | majorAlleleCopyNumber %%1>=0.8)
  
  subclones = d[, c(1:3)] # chr, start, end
  colnames(subclones) = c("chr", "startpos", "endpos")
  subclones$normal_total = 2
  if(sex=="male") subclones$normal_total[subclones$chr%in%c("chrX","chrY")] = 1
  subclones$nMaj1_A = round(d$majorAlleleCopyNumber)
  subclones$nMin1_A = round(d$minorAlleleCopyNumber)
  subclones$frac1_A = 1
  subclones$nMaj2_A = NA
  subclones$nMin2_A = NA
  subclones$frac2_A = NA
  
  subclones = subclones %>% filter(!is.na(frac1_A))# %>% mutate(chr=str_remove(chr,"chr"))
  
  write_tsv(subclones, paste(outfile.folder, samplename, "_subclones.txt", sep = ""))
  
  samplestats = read_tsv(copynumber.purple.summary) %>% filter(tumor_id==samplename)
  purity_ploidy = array(NA, c(3, 5))
  colnames(purity_ploidy) = c("rho", "psi", "ploidy", "distance", "is.best")
  rownames(purity_ploidy) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
  purity_ploidy["FRAC_GENOME", "rho"] = as.numeric(samplestats$purity )
  purity_ploidy["FRAC_GENOME", "psi"] = as.numeric(samplestats$ploidy)
  purity_ploidy["FRAC_GENOME", "ploidy"] = as.numeric(samplestats$ploidy)
  purity_ploidy["FRAC_GENOME", "is.best"] = TRUE
  write.table(purity_ploidy, paste(outfile.folder, samplename, "_rho_and_psi.txt", sep = ""), quote = F, sep = "\t")
}

# get options
library(optparse)
if(!exists("opt")){
  option_list = list(
    make_option(c("-n", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
    #make_option(c("-a", "--alignfolder"), type="character", default=".", help="Folder with BAM/CRAM files to extract mutation allele counts from, there must be an index with the same name, but with an extra .bai or .crai extension", metavar="character"),
    make_option(c("-v", "--vcffolder"), type="character", default=NULL, help="Folder with VCF files with somatic mutation data", metavar="character"),
    make_option(c("-c", "--CNV"), type="character", default=".", help="Path to CNV segmentation file [default= %default]", metavar="character"),
    make_option(c("-t", "--tumor"), type="character", default=".", help="Path to tumor BAM/CRAM file [default= %default]", metavar="character"),
    make_option(c("-s","--CNVsummary"), type="character", default=NULL, help="Facets samplestatistics output file", metavar="character"),
    #make_option(c("-o", "--output"), type="character", default=".", help="Output directory", metavar="character"),
    make_option(c("-f","--fai"), type="character", default=NULL, help="Reference genome index", metavar="character"),
    make_option(c("-g","--ign"), type="character", default=NULL, help="File with a list of contigs to ignore", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

## Load libraries
library(dpclust3p)
library(DPClust)
library(tidyverse)

samplename  = opt$samplename
# find cram file associated with samplename
#crams = opt$tumor #read_tsv(opt$TNfile) %>% filter(str_detect(sampleID,samplename)) %>% pull(tumor)
align_files = opt$tumor #grep(list.files(path=opt$alignfolder,pattern = paste(crams,collapse = "|") ),pattern = "ai$",value = T,invert = T)
samplenames = opt$samplename #read_tsv(opt$TNfile) %>% filter(str_detect(sampleID,samplename)) %>% pull(sampleID) #sapply(align_files, function(x) str_split(x,pattern = "_")[[1]][1] )
vcf_file    = opt$vcffolder
copynumber.purple.files = opt$CNV
CNVsummary = read_tsv(opt$CNVsummary)
sex = tolower(unique(CNVsummary$gender[CNVsummary$tumor_id==samplename])) # Get sex from CNVsummary 
output_dir = "./"
fai_file   = opt$fai
ign_file   = opt$ign

# check if files exist
.checkfile = function(infile) {
  if (!file.exists(infile)) {
    stop(paste("File", infile, "does not exist", sep=" "))
  }
}

.checkfile(vcf_file)
.checkfile(fai_file)
.checkfile(ign_file)

# Define the final output file
dpoutput_files = file.path(output_dir, paste(samplenames, "_allDirichletProcessInfo.txt", sep=""))

# Define various temp files
loci_file = file.path(output_dir, paste(samplename, "_loci.txt", sep=""))
allelecounts_files = file.path(output_dir, paste0(samplenames, "_alleleFrequencies.txt", sep=""))

# Dump loci - this function can take multiple vcf files when multiple samples from same donor
vcf2loci(vcf_file=vcf_file, fai_file=fai_file, ign_file=ign_file, outfile=loci_file)

for(i in 1:length(allelecounts_files)){
  print(i)
  # Fetch allele counts
  cmd = paste("~/bin/alleleCounter", "-b", paste0(opt$alignfolder,align_files[i]), "-o", allelecounts_files[i], 
              "-l",loci_file, "-m", 20, "-q", 35, sep = " ")
  system(cmd, wait = T)
  print("done alleleCounter")
  # Transform ASCAT NGS output into Battenberg output
  PurpleToBattenberg(samplenames[i], output_dir, copynumber.purple.files, opt$CNVsummary )
  print("done conversion to battenberg")
  
  # Create dpIn file
  print(c(loci_file=loci_file, allele_frequencies_file=allelecounts_files[i], 
           cellularity_file=paste(output_dir, samplenames[i], "_rho_and_psi.txt", sep = ""), 
           subclone_file=paste(output_dir, samplenames[i], "_subclones.txt", sep = ""), 
           gender=sex, SNP.phase.file="NA", mut.phase.file="NA", output_file=dpoutput_files[i]))
  
  runGetDirichletProcessInfo(loci_file=loci_file, allele_frequencies_file=allelecounts_files[i], 
                             cellularity_file=paste(output_dir, samplenames[i], "_rho_and_psi.txt", sep = ""), 
                             subclone_file=paste(output_dir, samplenames[i], "_subclones.txt", sep = ""), 
                             gender=sex, SNP.phase.file="NA", mut.phase.file="NA", output_file=dpoutput_files[i])
  print("done run getdirichletprocessinfo")
}

# remove the "chr" prefix from the DPclust files to avoid incompatibilities with the package
system("for file in `ls *allDirichletProcessInfo.txt`; do sed -i '2,$s/chr//g' $file; done")


#remove alterations in low-coverage regions (<60X) 
input_files = list.files("./",pattern = "allDirichletProcessInfo.txt",full.names = T)

for(i in 1:length(input_files)){
  input.tmp.files = input_files[i]
  input.tmp = lapply( input.tmp.files , read_tsv)
  ALT.tmp = sapply(input.tmp,function(x) x$mut.count )
  DP.tmp = sapply(input.tmp,function(x) x$WT.count + x$mut.count )
  CNT.tmp = sapply(input.tmp,function(x) x$subclonal.CN )
  CNmaj1.tmp = sapply(input.tmp,function(x) x$nMaj1 )
  CNmin1.tmp = sapply(input.tmp,function(x) x$nMin1 )
  
  minDP = rep(60,length(input.tmp.files))
  tokeep = which( apply(DP.tmp,1,function(x)all(x>=minDP)) &
                    rowSums(is.na(CNmaj1.tmp))==0 & rowSums(ALT.tmp)>0  )
  
  for(j in 1:length(input.tmp)) write_tsv(input.tmp[[j]][tokeep,],input.tmp.files[j])
}

#create DPclust master files with info on file location for each sample
rp_files = list.files("./",pattern = "rho_and_psi",full.names = T)

opt = list(samplenames=paste0(samplenames,collapse = ",") , 
             donornames = paste(rep(str_remove(samplenames,"T$|_TU[0-9]$"), length( samplenames )),collapse = "," ),
             rho_and_psi=paste(rp_files,collapse = ","),
             sex=sex,
             output=paste0("", samplenames, "_DPmasterfile.txt") ) 

print(opt)

source("preproc_dpclust_master_file.R")


