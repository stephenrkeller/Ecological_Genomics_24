library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)


#vcf.thin <- read.pcadapt("/gpfs1/home/s/r/srkeller/courses/Ecological_Genomics_24/population_genomics/outputs/vcf_final.filtered.thinned.lfmm",type="lfmm")

vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type="vcf")

# Alt if above doesn't work b/c read vcf is deprecated in pcadapt
#vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.lfmm",
#                    type="lfmm")


pcadapt.pca = pcadapt(vcf,
                      K=2,
                      min.maf=0.01,
                      LD.clumping = list(size=500, thr=0.2))

vcfR <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),] 

#poplist.names = meta2$region

plot(pcadapt.pca, option="scores", 
     pop=meta2$region,
     i=1,j=2)

str(pcadapt.pca)

plot(pcadapt.pca, option="stat.distribution") + 
  xlim(0,50) + ylim(0,0.1)

plot(pcadapt.pca, option="qqplot")

#plot(pcadapt.pca, option="manhattan")

hist(pcadapt.pca$pvalues, breaks=50, col="gold", 
     xlab="p-values")

padj = p.adjust(pcadapt.pca$pvalues, method="BH")
hist(padj, breaks=50)

outliers = which(padj < 0.05) # FDR
length(outliers)


# which PCs are the outliers associated with?
outliers.pc <- get.pc(pcadapt.pca,outliers)
table(outliers.pc$PC)

#plot(pcadapt.pca$loadings[,2]) # for PC1 SNP loadings

 
# Bring in VCF file to get info on Chrom/POS
View(head(vcfR@fix))

vcfR.fix <- as.data.frame(vcfR@fix[,1:2])
chr.main <- unique(vcfR.fix$CHROM)[1:8]

# We then use the 'seq' function to number the chromosomes from 1 to 8
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

# And finally merge these numerically named chromosomes with the full diversity results:
Pval <- pcadapt.pca$pvalues
Qval <- p.adjust(pcadapt.pca$pvalues, method="BH")

pcadapt.MHplot <- cbind(vcfR.fix,Pval,Qval)

pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, 
                            join_by(chr.main==CHROM))


# create a new "SNP" column that concatenates the chromosome ID with the bp position
pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP=paste0(chr.main,"_",POS))

# Make these fields numbers instead of characters:
pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$Pval = as.numeric(pcadapt.MHplot$Pval)
#pcadapt.MHplot$Qval = as.numeric(pcadapt.MHplot$Qval)


# get rid of NAs
pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(Pval)

outliers=pcadapt.MHplot$SNP[which(pcadapt.MHplot$Qval<0.005)]


# and finally....plot!!!
pdf("figures/pcadapt_K2_Manhattan.pdf", height=5,width=10)
manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="Pval",
          col=c("blue4","orange3"),
          logp=T,
          ylab="-log10 P-value",
          genomewideline = F,
          suggestiveline = 20,
          highlight = outliers,
          main="PCAdapt genome scan for selection (K=2)")
dev.off()

# do for component-wise PCAdapt run:
# component-wise runs for each pc separately
pcadapt.pca2 = pcadapt(vcf, 
                       K=2, 
                       method="componentwise",
                       min.maf=0.01,
                       LD.clumping = list(size=500, thr=0.2))
summary(pcadapt.pca2)
str(pcadapt.pca2)

Pval2 <- pcadapt.pca2$pvalues

pcadapt.MHplot2 <- cbind(vcfR.fix,Pval2)

pcadapt.MHplot2 <- left_join(chrnum, pcadapt.MHplot2, 
                            join_by(chr.main==CHROM))


# create a new "SNP" column that concatenates the chromosome ID with the bp position
pcadapt.MHplot2 <- pcadapt.MHplot2 %>%
  mutate(SNP=paste0(chr.main,"_",POS))

# Make these 2 fields numbers instead of characters:
pcadapt.MHplot2$V2 = as.numeric(pcadapt.MHplot2$V2)
pcadapt.MHplot2$POS = as.numeric(pcadapt.MHplot2$POS)
pcadapt.MHplot2$pPC1 = as.numeric(pcadapt.MHplot2[,4])
pcadapt.MHplot2$pPC2 = as.numeric(pcadapt.MHplot2[,5])

# drop NAs
pcadapt.MHplot2 <- pcadapt.MHplot2 %>% drop_na(pPC1)

# and finally....plot!!!
pdf("figures/pcadapt_PC1_2_Manhattan.pdf", height=10,width=10)
par(mfrow=c(2,1))
manhattan(pcadapt.MHplot2,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4","orange3"),
          logp=T,
          ylab="-log10 p-value",
          genomewideline = F,
          suggestiveline = 13,
          main="PCAdapt genome scan for selection (PC1)")

manhattan(pcadapt.MHplot2,
          chr="V2",
          bp="POS",
          p="pPC2",
          col=c("blue4","orange3"),
          logp=T,
          ylab="-log10 p-value",
          genomewideline = F,
          suggestiveline = 13,
          main="PCAdapt genome scan for selection (PC2)")

dev.off()

# SNP2GO annotation and enrichment testing
library(SNP2GO)
library(GenomicRanges)

# extract outliers and Non-outliers from vcf.fix

outliers=Qval<0.05

vcfR.fix.outliers <- vcfR.fix[outliers==TRUE,]
vcfR.fix.outliers$POS <- as.numeric(vcfR.fix.outliers$POS)
vcfR.fix.outliers <- vcfR.fix.outliers %>% drop_na()

candSNPs <- GRanges(seqnames=vcfR.fix.outliers[,1],
                    ranges=IRanges(vcfR.fix.outliers[,2],
                                   vcfR.fix.outliers[,2]))

gff=read.table("/gpfs1/cl/pbio3990/PopulationGenomics/reference/GCA.gff", sep="\t",quote=)
gff2 <- gff[which(gff$V3=="gene"),c(4,5)]

vcfR.fix.Nonoutliers <- vcfR.fix[-outliers,]
vcfR.fix.Nonoutliers$POS <- as.numeric(vcfR.fix.Nonoutliers$POS)
noncandSNPs <- GRanges(seqnames=vcfR.fix.Nonoutliers[,1],
                    ranges=IRanges(vcfR.fix.Nonoutliers[,2],
                                   vcfR.fix.Nonoutliers[,2]))

pcadapt_snp2go <- snp2go(gff="/gpfs1/cl/pbio3990/PopulationGenomics/reference/GCA.gff",
                         candidateSNPs = candSNPs,
                         noncandidateSNPs = noncandSNPs,
                         FDR=0.05,
                         runs=100000,
                         extension=50000)