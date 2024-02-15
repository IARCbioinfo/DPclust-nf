library(MutationTimeR)
library(tidyverse)
library(patchwork)

# get options
library(optparse)
if(!exists("opt")){
  option_list = list(
    make_option(c("-n", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
    make_option(c("-v", "--vcf"), type="character", default=NULL, help="Path to VCF files with somatic mutation data", metavar="character"),
    make_option(c("-c", "--CNV"), type="character", default=".", help="Path to CNV segmentation file [default= %default]", metavar="character"),
    make_option(c("-s","--CNVsummary"), type="character", default=NULL, help="PURPLE sample statistics output file", metavar="character"),
    make_option(c("-k","--clusters"), type="character", default=NULL, help="DPclust generated clusters file", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

svp <- ScanVcfParam(info=c("Func.ensGene","Gene.ensGene","ExonicFunc.ensGene","AAChange.ensGene"),geno=c("GT","AD","AF","DP"))
vcf <- readVcf(opt$vcf,"hg38", svp)

GT = unlist(geno(vcf)$GT)
AD = t(apply(geno(vcf)$AD, 1, function(x) x[[1]]))
AD2 = array(AD,dim=c(nrow(AD),1,2))
AF = unlist(geno(vcf)$AF)
DP = unlist(geno(vcf)$DP)
geno(vcf)$AD = AD2

CNVs = read_tsv(opt$CNV) %>% filter(sample==opt$samplename)
CNVsummary = read_tsv(opt$CNVsummary) %>% filter(tumor_id==opt$samplename)

bb <- GRanges(seqnames=CNVs$chromosome , ranges = IRanges(start=CNVs$start,end=CNVs$end), 
              major_cn=CNVs$majorAlleleCopyNumber.corrected.integer , minor_cn=CNVs$majorAlleleCopyNumber.corrected.integer, 
              clonal_frequency=CNVsummary$purity) # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment

clusters <- read.table(opt$clusters,header = T) # Optional data.frame with subclonal cluster locations (VAF proportion) and size (number of variants n_ssms)

#


if(sum(vcf%over% bb)>20){
if(nrow(clusters)==0){#no subclones identified
  mt <- mutationTime(vcf, bb, n.boot=10)
}else{
  mt0 <- mutationTime(vcf, bb, n.boot=10)
  mt <- mutationTime(vcf, bb, n.boot=10,clusters=clusters)
  # visualise alterations and CCF
  gg0 = ggplot(data_frame(VAF = AD2[,,2]/DP[,1], CLS=mt0$V$CLS), aes(x=VAF,fill=CLS) ) + geom_histogram(bins = 100) + theme_classic() +ggtitle("Automatic subclone detection")
  gg1 = ggplot(data_frame(VAF = AD2[,,2]/DP[,1], CLS=mt$V$CLS), aes(x=VAF,fill=CLS) ) + geom_histogram(bins = 100) + theme_classic() +ggtitle("DPclust subclone detection")
  ggsave(filename = paste0(opt$samplename,"_VAFdistribution_timing.pdf"),gg0+gg1,width = 3.5*2,height = 3.5)
}

vcf2 <- addMutTime(vcf, mt$V)
# remove chr
seqlevels(rowRanges(vcf2)) = str_remove(seqlevels( rowRanges(vcf2)),"chr")

vcf20 <- addMutTime(vcf, mt0$V)
# remove chr
seqlevels(rowRanges(vcf20)) = str_remove(seqlevels( rowRanges(vcf20)),"chr")

bb2 = bb
seqlevels(bb2) = str_remove(seqlevels( bb2),"chr")
mcols(bb2) <- cbind(mcols(bb2),mt$T)

bb20 = bb
seqlevels(bb20) = str_remove(seqlevels( bb20),"chr")
mcols(bb20) <- cbind(mcols(bb20),mt0$T)

# get hg38 ref lengths
regions = GRanges(seqnames= c(1:22,"X","Y"), 
                  ranges=IRanges(start = 1, end = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,
                                                    133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,
                                                    64444167,46709983,50818468,156040895,57227415)))

pdf(paste0(opt$samplename,"_mutationtimeR.pdf"))
plotSample(vcf2,bb2,regions = regions,title = opt$samplename)
dev.off()

pdf(paste0(opt$samplename,"_mutationtimeR_autosubclone.pdf"))
plotSample(vcf20,bb20,regions = regions,title = opt$samplename)
dev.off()

write_tsv(as.data.frame(mt$V),file = paste0(opt$samplename,"_SmallVariant_timing.tsv"))
write_tsv(as.data.frame(bb2),file = paste0(opt$samplename,"_CNV_timing.tsv"))

write_tsv(as.data.frame(mt0$V),file = paste0(opt$samplename,"_SmallVariant_timing_autosubclone.tsv"))
write_tsv(as.data.frame(bb20),file = paste0(opt$samplename,"_CNV_timing_autosubclone.tsv"))

# for high-coverage regions only
mt.deep <- mutationTime(vcf[DP >=60,], bb, n.boot=10,clusters=clusters)
vcf2.deep <- addMutTime(vcf[DP >=60,], mt.deep$V)
# remove chr
seqlevels(rowRanges(vcf2.deep)) = str_remove(seqlevels( rowRanges(vcf2.deep)),"chr")

bb2.deep = bb
seqlevels(bb2.deep) = str_remove(seqlevels( bb2.deep),"chr")
mcols(bb2.deep) <- cbind(mcols(bb2.deep),mt.deep$T)

pdf(paste0(opt$samplename,"_mutationtimeR_60Xcov.pdf"))
plotSample(vcf2.deep,bb2,regions = regions,title = opt$samplename)
dev.off()

write_tsv(as.data.frame(mt.deep$V),file = paste0(opt$samplename,"_SmallVariant_timing_60Xcov.tsv"))
write_tsv(as.data.frame(bb2.deep),file = paste0(opt$samplename,"_CNV_timing_60Xcov.tsv"))


}else{
  write_tsv(data.frame(comment="too_few_variants"),file = paste0(opt$samplename,"_SmallVariant_timing.tsv"))
  write_tsv(data.frame(comment="too_few_variants"),file = paste0(opt$samplename,"_CNV_timing.tsv"))
  write_tsv(data.frame(comment="too_few_variants"),file = paste0(opt$samplename,"_SmallVariant_timing_60Xcov.tsv"))
  write_tsv(data.frame(comment="too_few_variants"),file = paste0(opt$samplename,"_CNV_timing_60Xcov.tsv"))
}
