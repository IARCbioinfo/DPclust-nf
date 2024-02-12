library(MutationTimeR)
library(tidyverse)

# get options
library(optparse)
if(!exists("opt")){
  option_list = list(
    make_option(c("-n", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
    make_option(c("-v", "--vcf"), type="character", default=NULL, help="Path to VCF files with somatic mutation data", metavar="character"),
    make_option(c("-c", "--CNV"), type="character", default=".", help="Path to CNV segmentation file [default= %default]", metavar="character"),
    make_option(c("-s","--CNVsummary"), type="character", default=NULL, help="Facets samplestatistics output file", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# input data
opt = list()
opt$vcf = "/data/lungNENomics/work/organoids/WGS/variant_calling/release3_25102021/somatic/split_per_sample/reannotated_11102022/LNET6T_final_deannotated.hg38_multianno.vcf.gz"
opt$CNV = "/data/lungNENomics/work/alcalan/WGS/CNV/release1_PURPLE_11052022/lungNENomics.all.purple.cnv.somatic.10052022.tsv"
opt$CNVsummary = "/data/lungNENomics/work/alcalan/WGS/evolution/purple_summary.tsv"
opt$samplename = "LNET6T"  

svp <- ScanVcfParam(info=c("Func.ensGene","Gene.ensGene","ExonicFunc.ensGene","AAChange.ensGene"),geno=c("GT","AD","AF","DP"))
vcf <- readVcf(opt$vcf,"hg38", svp)

GT = unlist(geno(vcf)$GT)
AD = array(unlist(geno(vcf)$AD),dim = c(nrow(geno(vcf)$AD),1,2))
AF = unlist(geno(vcf)$AF)
DP = unlist(geno(vcf)$DP)
geno(vcf)$AD = AD

#vcf <- MutationTimeR::readVcf(opt$vcf) # Point mutations, needs `geno` entries `AD` and `DP` or `info` columns t_alt_count t_ref_count.
CNVs = read_tsv(opt$CNV) %>% filter(sample==opt$samplename)
CNVsummary = read_tsv(opt$CNVsummary) %>% filter(tumor_id==opt$samplename)

bb <- GRanges(seqnames=CNVs$chromosome , ranges = IRanges(start=CNVs$start,end=CNVs$end), 
              major_cn=CNVs$majorAlleleCopyNumber.corrected.integer , minor_cn=CNVs$majorAlleleCopyNumber.corrected.integer, 
              clonal_frequency=CNVsummary$purity) # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment

clusters <- data.frame(cluster= , n_ssms=, proportion=) # Optional data.frame with subclonal cluster locations (VAF proportion) and size (number of variants n_ssms)

#
mt <- mutationTime(vcf, bb, n.boot=10) #clusters=clusters, 

vcf2 <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

plotSample(vcf2,bb)
