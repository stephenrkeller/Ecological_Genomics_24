# load libraries
library(vcfR)
library(ggplot2)

# set wd to PopulationGenomics folder on gpfs1 class drive

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

list.files("variants/")
list.files("reference/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")
#dna2 = dna$CM058040.1

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")
#gff2 = grep("CM058040.1",gff)

vcf

head(vcf)

# create the chromR object class
chr1 <- create.chromR(name="Centaurea chr1", vcf=vcf, seq=dna, ann=gff, verbose=F)
#doesn't work: chr2 <- create.chromR(name="Centaurea chr2", vcf=vcf, seq=dna2, ann=gff2, verbose=F)

#head(chr1@var.info)
#chr1
#quantile(chr1@var.info$MQ)
#quantile(chr1@var.info$DP)

# example plot
plot(chr1)
chromoqc(chr1, xlim=c(1e1, 1.1e8))

# example with exporting to pdf -- note path requirements
# also note zoom ability via xlim=c() option
pdf(file="~/courses/Ecological_Genomics_24/population_genomics/figures/chromoPlot_chr1.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()


# Mask poor quality variants
chr1_masked <- masker(chr1,
                      min_QUAL=50,
                      min_DP=1000,
                      max_DP=10000,
                      min_MQ=30)

plot(chr1_masked)
chromoqc(chr1_masked, xlim=c(1e1, 1.1e8))

# Now process the chromR object with proc.chromR
# The default window size = 1000bp, can set with win.size
chr1_proc <- proc.chromR(chr1_masked, win.size = 1e5)
plot(chr1_proc)
chromoqc(chr1_proc, xlim=c(1e1, 1.1e8))

#head(chr1_masked) # see an overview of each list element


#How to filter a vcf file for minDP
dp <- extract.gt(vcf, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)
quantile(dp)

dp[dp==0] <- NA

quantile(dp, na.rm=T)

#pdf(file="~/courses/Ecological_Genomics_24/population_genomics/figures/chromoPlot_chr1.pdf")
heatmap.bp(dp[1:1000,], rlabels=F, clabels=F)
#dev.off()

vcf@gt[,-1][is.na(dp)==TRUE] <- NA

indMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
indMiss <- indMiss/nrow(vcf@gt)
hist(indMiss)
quantile(indMiss)

siteMiss <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) })
siteMiss <- siteMiss/ncol(vcf@gt)
hist(siteMiss)

vcf2 = vcf

vcf2@gt = vcf2@gt[,c(TRUE, indMiss < 0.7)]

vcf2

siteMiss <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) })
siteMiss <- siteMiss/ncol(vcf@gt)
hist(siteMiss)

vcf2 = vcf2[siteMiss < 0.25,]

vcf2

heatmap.bp(dp2[1:1000,], rlabels=F, clabels=F)

#the third is the count and fourth the frequency.
vcf2_maf <- as.data.frame(maf(vcf2, element=2)) 
# element=2 returns info on the minor allele
hist(vcf2_maf$Frequency, breaks=50)

quantile(vcf2_maf$Frequency)

vcf3 = vcf2

vcf3_maf <- maf(vcf2, element=2)

vcf3 = vcf3[which(vcf3_maf$Frequency>0.05),]


write.vcf(vcf, "~/test.vcf.gz")


# indDP <- colMeans(dp2)
# hist(dp2)
# hist(indDP)






myINFO = c("INDEL","DP","AC","AN","MQ") 
#AC=count of Alt alleles (for doing --mac)
#AN=count of total (REF+ALT) alleles in called genotypes (after masking??)

vcf_tdy <- vcfR2tidy(vcf, 
                     info_fields=myINFO,
                     dot_is_NA=TRUE)

vcf_tdy$meta
vcf_tdy$fix
vcf_tdy$gt


print(vcf_tdy$gt, n=50)  # view `n` rows in the tibble

# Can extract fields from the `gt` using either vcfR or chom
indDP = extract.gt(chr1, 
                   element="DP", 
                   mask=TRUE, 
                   as.numeric=TRUE, 
                   convertNA=TRUE) 

quantile(indDP) # get quantiles
hist(indDP, xlim=c(0,20),breaks=100) # plot histogram
mean(indDP) # mean 20X


# The first column is the total number of alleles, 
# the second is the number of NA genotypes, 
#the third is the count and fourth the frequency.
vcf_maf <- as.data.frame(maf(vcf, element=2)) 
# element=2 returns info on the minor allele
hist(vcf_maf$Frequency, breaks=50)

vcf_tdy_gt = extract_gt_tidy(vcf,
                             format_fields = c("DP","GT", "PL"),
                             format_types = TRUE,
                             dot_is_NA = TRUE,
                             verbose = TRUE)


PL = tidyr::separate_wider_delim(vcf_tdy$gt, cols="gt_PL", ",", names=c("PL1","PL2","PL3"))

pl <- extract.gt(vcf, element = "PL")
pl_1 <- masplit(pl, sort=1, decreasing=0)
pl_2 <- masplit(pl, record = 2)
pl_3 <- masplit(pl, record = 3)


#######################################

# Playing with the annotation file
SEARCHTERM="FRI"
chr1_masked_processed@ann[grep(SEARCHTERM, chr1_processed@ann$V9),4:5] # access annotation info

head(chr1_masked_processed@var.info) # access variant table info

head(chr1_masked_processed@win.info) # access info on windows


# what if we want to ignore the windows with no SNPs?
chr1_masked_processed@win.info$variants[chr1_masked_processed@win.info$variants==0] <- NA

hist(chr1_masked_processed@win.info$variants)
summary(chr1_masked_processed@win.info$variants)

ggplot(as.data.frame(chr1_masked_processed@win.info), aes(x=start, y=variants)) +
  geom_point(size=2, shape=2, color="red")

#chr1_masked_processed@win.info$genic # what is this?


# We can write out the VCF with masking
write.vcf(chr1_masked_processed, 
          "~/courses/Ecological_Genomics_24/population_genomics/outputs/myfiltered.vcf.gz",
          mask=TRUE)

vcf_field_names(vcf) # gets the INFO names and def's so you can choose what you want



############  NOT WORKING #################

GQ = median(as.numeric(PL[,4:6]))  # not working...

# Looking at ploidy info in allele dosage frequencies
gt <- extract.gt(vcf)
hets <- is_het(gt)
# Censor non-heterozygous positions.
is.na(vcf@gt[,-1][!hets]) <- TRUE

# Extract allele depths.
ad <- extract.gt(vcf, element = "AD")
ad1 <- masplit(ad, record = 1)
ad2 <- masplit(ad, record = 2)
freq1 <- ad1/(ad1+ad2)
freq2 <- ad2/(ad1+ad2)
myPeaks1 <- freq_peak(freq1, getPOS(vcf))  # breaks on this step :()
is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
myPeaks2 <- freq_peak(freq2, getPOS(vcf), lhs = FALSE)
is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE
myPeaks1



# calculates allele frequencies per genotypes for assessing ploidy assumptions
AD_frequency(ad, delim = ",", allele = 1L, sum_type = 0L, decreasing = 1L) 