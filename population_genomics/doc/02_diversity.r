library(vcfR)
        
########

# 02_diversity.r

vcf.div <- readvcfR("pathto myfiletered.vcf.gz")

meta3 <- read.csv("pathtomymeta_filtered.csv")

# make sure same dimensions!
vcf.div
str(meta3)


mydiff <- genetic_diff(vcf.div,
                             pops=as.factor(meta3$region),
                             method="nei")

View(mydiff)

# Use tidy to extract chromosomes and stats

unique(mydiff$CHROM)

mydiff.chr1 = mydiff %>% 
  select(CHROM,POS,Ht,Gst) %>%
  filter(CHROM=="CM058040.1")

plot(mydiff.chr1$POS,mydiff.chr1$Ht,
     pch=21, col="lightblue",
     xlab="Genomic position (bp)",
     ylab="Total genetic diversity (Ht)")



# knitr::kable(round(colMeans(myDiff[,c(3:9,16)], na.rm = TRUE), digits = 3))

###  In PopGenomeR

library(PopGenome)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

dat <- readVCF("variants/Centaurea_filtered.vcf.gz", 
               numcols=20000,
               tid="CM058040.1",
               frompos=1,
               topos=1.1e8,
               include.unknown=T,
               approx=T,
               gffpath="reference/GCA.gff")

get.sum.data(dat)

dat.slide = sliding.window.transform(dat, 10000, 10000, type=2)

length(dat.slide@region.names)
out = neutrality.stats(dat.slide)
get.neutrality(out)
hist(out@Tajima.D)
hist(out@theta_Tajima/10000)
hist(out@theta_Watterson/10000)
hist(out@n.segregating.sites)

div = diversity.stats(dat.slide)
get.diversity(div)
plot(div@nuc.diversity.within/div@n.sites, 
     pch=21, 
     col="blue", 
     cex=0.7,
     xlab="Chromosome 1 window (10 kb)",
     ylab="Nucleotide diversity (per-site)")

sfs = detail.stats(dat,site.spectrum = T)
get.detail(sfs)
allelefrqs = sfs@region.stats@minor.allele.freqs[[1]]
hist(allelefrqs,breaks=50)

dat@n.biallelic.sites

genes = splitting.data(dat, subsites="gene")

length(genes@region.names)

genediv = diversity.stats(genes)
hist(genediv@nuc.diversity.within[genediv@nuc.diversity.within!=0]/genediv@n.sites[genediv@nuc.diversity.within!=0],xlim=c(0,0.0005),breaks=500)

gene.neut = neutrality.stats(genes)

plot(gene.neut@n.segregating.sites, 
     pch=21, 
     col="blue", 
     cex=0.7,
     xlab="Chromosome 1 window (10 kb)",
     ylab="Nucleotide diversity (per-site)")


genediv2 = neutrality.stats(genes)
plot(genediv2@SLIDE.POS,genediv2@Tajima.D)

plot(genediv2@theta_Watterson[is.na(genediv2@n.segregating.sites==0)==F]/genediv2@n.sites[is.na(genediv2@n.segregating.sites==0)==F])

