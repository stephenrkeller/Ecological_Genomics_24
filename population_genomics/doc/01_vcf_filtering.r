# load libraries
library(vcfR)
#library(ggplot2)

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
# chr1_masked <- masker(chr1,
#                       min_QUAL=50,
#                       min_DP=1000,
#                       max_DP=10000,
#                       min_MQ=30)
# 
# plot(chr1_masked)
# chromoqc(chr1_masked, xlim=c(1e1, 1.1e8))
# 
# # Now process the chromR object with proc.chromR
# # The default window size = 1000bp, can set with win.size
# chr1_proc <- proc.chromR(chr1_masked, win.size = 1e5)
# plot(chr1_proc)
# chromoqc(chr1_proc, xlim=c(1e1, 1.1e8))

#head(chr1_masked) # see an overview of each list element

####### Left off here last time #######

#How to filter a vcf file for minDP
DP <- extract.gt(vcf, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)
quantile(DP)

DP[DP==0] <- NA. #set the 0 depth genos to `NA`

quantile(DP, na.rm=T) # much better

dim(DP) # ensure loci are in rows; samples in columns

#Let's look at mean DP per individual...
# avgDP = colMeans(DP, na.rm=T)
# summary(avgDP)
# hist(avgDP, breaks=50) 
# mean is ~24X, range from 2.8-171X.  
# Pretty good! Ideeally want avg of 15-20X/ind

#What about missingness? We can use the heatmap function:

#pdf(file="~/courses/Ecological_Genomics_24/population_genomics/figures/chromoPlot_chr1.pdf")
heatmap.bp(DP[1:5000,], rlabels=F, clabels=F)
#dev.off()

# set individual genotypes with DP<X to `NA`
# not needed with SNPfiltR??
#vcf@gt[,-1][is.na(DP)==TRUE] <- NA 

vcf # check to see % missing data -- 25.3% if DP==0 <- NA

# Now that we see the data attribtues, let's start filtering 
library(SNPfiltR)

meta <- read.csv("metadata/meta4vcf.txt, header=T")
head(meta)
meta2=meta[,c(1,4)]
names(meta2) = c("id","pop")

# Look at mean depth per ind
hard_filter(vcf)
vcf.filt <- hard_filter(vcf, 
                        depth=3) #What's reasonable while keeping variants?

# Look at allele balance (note these are autotetratploid...0.25, 0.5, 0.75)
vcf.filt <- filter_allele_balance(vcf.filt,
                                  min.ratio = 0.15,
                                  max.ratio=0.85)

vcf.filt <- max_depth(vcf.filt, 
                      maxdepth=60) # generally set filter to 2X mean depth 

# start with no cutoff for exploratory, 
# then add cuttoff=0.8 or 0.75 (N=36 sammples)

vcf.filt.indMiss <- missing_by_sample(vcf.filt, 
                  popmap = meta,
                  cutoff=0.8) 
#subset popmap to only include retained individuals
meta <- meta[meta$id %in% colnames(vcf.filt.indMiss@gt),]

# gets rid of monomrophic or multi-allelic sites
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss) 


# use PCA to verify that missingngess is not driving clustering
# maybe try just 1 or 2 threholds at first...note that lower is more permissive
# didn't work for me on the VACC-OOD!
#library(adgenet)
#missPCA <- assess_missing_data_pca(vcfR=vcf.filt.indMiss, 
                              #popmap = meta, 
                              #thresholds = 0.5, 
                              #clustering = FALSE)

#Filter out by SNP missingness -- higher cuttoff is more stringent
vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, 
                                      cutoff=0.5)

vcf.filt.indSNPMiss <- min_mac(vcf.filt.indSNPMiss,
                               min.mac=2)

#assess clustering without MAC cutoff
#miss<-assess_missing_data_tsne(vcf.filt.indSNPMiss, 
                              # popmap=meta, 
                               #clustering = FALSE)

DP2 <- extract.gt(vcf.filt.indSNPMiss, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)


heatmap.bp(DP2[1:5000,], rlabels=F, clabels=F)

vcfR::write.vcf(vcf.filt.indSNPMiss, 
                "~/myrepo/outputs/Centaurea_finalfiltered.vcf.gz")

# Can also thin for LD:

vcf.filt.indSNPMiss.thin <- distance_thin(vcf.filt.indSNPMiss,
                                          min.distance=500)

vcfR::write.vcf(vcf.filt.indSNPMiss.thin, 
                "~/myrepo/outputs/Centaurea_finalfiltered_thinned.vcf.gz")

mydiff <- vcfR::genetic_diff(vcf.filt.indSNPMiss,
                             pops=as.factor(meta$region),
                             method="nei")

#if there are still problematic samples, drop them using the following syntax
#vcfR <- vcfR[,colnames(vcfR@gt) != "mysampleID" & colnames(vcfR@gt) != "mysampleID2"]





#### Below is optional, depending if SNPfiltR poackage works

#optional if class wants to filter:
#DP[DP<3 <- NA]

vcf@gt[,-1][is.na(DP)==TRUE] <- NA # set individual genotypes with DP<X to `NA`

DP.filtered <- extract.gt(vcf, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)

quantile(DP.filtered, na.rm=T) # good -- our DP filtering step worked (min=3)

#Now, calculate how many missing (NA) genos per individual
indMiss <- apply(DP, 
                 MARGIN = 2, 
                 function(x){ sum(is.na(x)) })

indMissProp <- indMiss/nrow(vcf@gt) #make a proportion

hist(indMissProp)
summary(indMissProp)

siteMiss <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) })
siteMiss <- siteMiss/ncol(vcf@gt)
hist(siteMiss)

vcf2 = vcf 

vcf2@gt = vcf2@gt[,c(TRUE, indMiss < 0.7)]

vcf2

#vcf = vcf[siteMiss < 0.25,]

dp2 <- extract.gt(vcf2, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)
dp2[dp2==0] <- NA
quantile(dp2, na.rm=T)

siteMiss <- apply(dp2, MARGIN = 1, function(x){ sum(is.na(x)) })
siteMiss <- siteMiss/ncol(vcf2@gt)
hist(siteMiss)
quantile(siteMiss)

heatmap.bp(dp2[1:1000,], rlabels=F, clabels=F)

vcf2 = vcf2[siteMiss < 0.25,]

#the third is the count and fourth the frequency.
#vcf_maf <- as.data.frame(maf(vcf, element=2)) 
# element=2 returns info on the minor allele
#hist(vcf_maf$Frequency, breaks=50)

write.vcf(vcf2, "~/test.vcf.gz")


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