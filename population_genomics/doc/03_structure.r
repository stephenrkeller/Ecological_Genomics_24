library(vcfR)
library(LEA)
library(SNPfiltR)
library(tidyverse)

X11.options(type="cairo")
#options(bitmapType = "cairo")

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")


list.files("~/projects/eco_genomics/population_genomics/outputs/")
list.files("~/courses/Ecological_Genomics_24/population_genomics/outputs/")

setwd("~/courses/Ecological_Genomics_24/population_genomics")


vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
 
vcf.thin = distance_thin(vcf, min.distance=500)

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta) # meta has 629 inds

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]),] 

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

# Can't store these data in GH repo, so make a temp folder to store...
system("mkdir ~/temp_dat")

#Only need to do once -- save to outside your repo to avoid github max files size issues
system("gunzip -c ~/courses/Ecological_Genomics_24/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz >~/temp_dat/vcf_final.filtered.thinned.vcf")

mygenothinned <- vcf2geno(input.file="/gpfs1/home/s/r/srkeller/temp_dat/vcf_final.filtered.thinned.vcf", 
                   output.file="outputs/vcf_final.filtered.thinned2.geno")

# This will make a PCA project with the same name as the input file with extension .pca
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

show(CentPCA)

# Tracy Widon test -- basically out to PC 37 or so is significant...
tw = tracy.widom(CentPCA)
tw$percentage[1:5]
  
pdf("figures/PCA_3.6KSNPs.pdf")
# Make the PCA plot
plot(CentPCA$projections,
     pch=21,lwd=1.5,
     col=as.factor(meta2$region),
     xlab="Genetic PC1 (2.3%)",
     ylab="Genetic PC2 (1.1%)",
     main="PCA of Centaurea 3.6K SNPs")
# add the legend in
legend("bottomright", 
       legend=as.factor(unique(meta2$region)),
       fill=as.factor(unique(meta2$region)))
# Turn the pdf writer off
dev.off()


# Try smnf clustering

CentAdmix = snmf("outputs/vcf_final.filtered.thinned.geno", 
               K=1:5, 
               entropy=T, 
               repetitions=3,
               project="new")

plot(CentAdmix, col="blue4",cex=1.5,pch=19)

myK=4
  
CE = cross.entropy(CentAdmix, K=myK)

best=which.min(CE)

myKQ = Q(CentAdmix, K=myK, run=best)
myKQmeta = cbind(myKQ, meta2)

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(levels=region, .by_group=TRUE)

my.colors <- c("blue4","gold","tomato","lightblue","olivedrab")

pdf("figures/Admixture_thinned_K4.pdf", height=4, width=10)
barplot(as.matrix(t(myKQmeta[,1:myK])), 
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab="Geographic regions",
        ylab="Ancestry proportions",
        main = paste0("Ancestry matrix K=",myK),
        cex.names=0.7)
axis(1, 
     at=1:length(myKQmeta$region),
     labels=myKQmeta$region,
     tick=F,
     las=3)
dev.off()


# multi-K barplot
pdf("figures/Admixture_thinned_multK.pdf", height=10, width=10)
par(mfrow=c(3,1))

for(i in seq(2,4,1)) {
  myK=i
  myKQ = Q(CentAdmix, K=myK, run=best)
  barplot(as.matrix(t(myKQ[,1:myK])), 
        border=NA,
        space=0,
        col=my.colors[1:myK],
        ylab="Ancestry proportions",
        main = paste0("Ancestry matrix K=",myK),
        cex.names=0.7)
axis(1, 
     at=1:length(meta2$region),
     labels=meta2$region,
     tick=F,
     las=3)
}

dev.off()

###### using LEA::barchart function:

fig <- barchart(CentAdmix,
         K=myK, run=best,
         sort.by.Q = F,
         col=my.colors,
         lab=meta2$id,
         border=NA,
         space=0,
         ylab="Ancestry proportions",
         main = "K=4 ancestry matrix")
axis(1, 
     at=1:length(meta2$region),
     labels=meta2$region,
     las=3,
     srt=45,
     cex.axis=.7)
