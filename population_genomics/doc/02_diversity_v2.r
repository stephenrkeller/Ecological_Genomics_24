# Diversity estimation in Centaurea

library(vcfR)
library(tidyverse)
library(qqman)

#options(bitmapType='cairo')
X11.options(type="cairo")

# read in analysis-ready vcf file saved from our filtering steps
vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

# read in our metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

# always a good idea to check that the sample
# sizes are the same, given our filtering, etc.
length(colnames(vcf@gt[,-1])) # 529 samples
dim(meta) # 629 samples!  we need to subset the meta like we did last time...

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
dim(meta2) # 529 -- much better!

# Now, we can use the `genetic_diff` function in vcfR
# to get a lot of info on diversity in terms of heterozygosity (Hs for each pop, Ht for total sample)
# and in terms of Fst.  Let's do it...

vcf.div <- genetic_diff(vcf,
                        pops=as.factor(meta2$continent),
                        method="nei")

# Use your tools to take a look at the contents (View(), head(), str(), etc...)

#vcf.div.subset <- vcf.div[,c(1:5,8,11)]
#vcf.div.subset <- vcf.div[,]

chr <- unique(vcf.div$CHROM)
chr.main = chr[1:8]
chrnum = as.data.frame(cbind(chr.main, seq(1,length(chr.main),1)))

vcf.div.MHplot <- left_join(chrnum,
                            vcf.div,
                            join_by(chr.main==CHROM))


# set up the plotting df using tidyr:
# no negative Gst values, need a "SNP" column (basically chr + pos)

vcf.div.MHplot <- vcf.div.MHplot  %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main,"_",POS)) 

# check the var types:
str(vcf.div.MHplot)
levels(as.factor(vcf.div.MHplot$chr.main))

vcf.div.MHplot$V2 <- as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS <- as.numeric(vcf.div.MHplot$POS)


# main = "Manhattan Plot",

pdf("~/courses/Ecological_Genomics_24/population_genomics/figures/ManhattanPlot_Gst_byCont.pdf")
manhattan(vcf.div.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          cex = 1, 
          col = c("blue4", "orange3"), 
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999),
          chrlabs=as.character(c(1:8)),
          logp=F,
          ylab="Gst (EU vs. NA)",
          ylim=c(0,0.3))
dev.off()

outliers = vcf.div.MHplot %>%
  select(chr.main, POS, SNP, Gst) %>%
  filter(Gst > quantile(Gst, 0.999))

write.csv(outliers, 
          "~/courses/Ecological_Genomics_24/population_genomics/outputs/Gst_byCont_outliers.csv",
          quote=F,
          row.names = F)

vcf.div.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4,5)) %>% 
  ggplot(aes(POS, value, color = name)) + 
  geom_point(alpha=0.5) +
  facet_wrap(vars(chr.main), scales = "free_y", nrow = 4, strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(title="Genome-wide diversity", colour="Continent", x="Chromosomal position (bp)", y="Gene diversity (Hs)")

ggsave("Chromplot_Diversity_byCont.pdf", 
       path="~/courses/Ecological_Genomics_24/population_genomics/figures/")

ggplot(vcf.div.MHplot, aes(POS, Gst)) +
  geom_point() +
  facet_wrap(vars(chr.main), scales = "free_y", nrow = 4, strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  labs(title="Genome-wide differentiation", x="Chromosomal position (bp)",y="Fst (among continents)")

ggsave("Chromplot_Gst_Cont.pdf", 
       path="~/courses/Ecological_Genomics_24/population_genomics/figures/")


vcf.div.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4,5)) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title="Genome-wide diversity", fill="Continent", x="Gene diversity (Hs)",y="Frequency")

ggsave("Histogram_Diversity_byCont.pdf", 
       path="~/courses/Ecological_Genomics_24/population_genomics/figures/")

################################################################
##########  using a different grouping for the pop division ####
################################################################

vcf.div2 <- genetic_diff(vcf,
                         pops=as.factor(meta2$region),
                         method="nei")

# Use your tools to take a look at the contents (View(), head(), str(), etc...)

vcf.div2.subset <- vcf.div2[,c(1:9,16,19)]

chr <- unique(vcf.div2.subset$CHROM)
chr.main = chr[1:8]
chrnum = as.data.frame(cbind(chr.main, seq(1,length(chr.main),1)))

vcf.div2.sub.MHplot <- left_join(chrnum,
                                 vcf.div2.subset,
                                 join_by(chr.main==CHROM))

vcf.div2.sub.MHplot <- vcf.div2.sub.MHplot  %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr,"_",POS)) 

vcf.div2.sub.MHplot$POS <- as.numeric(vcf.div2.sub.MHplot$POS)
vcf.div2.sub.MHplot$V2 <- as.numeric(vcf.div2.sub.MHplot$V2)

# main = "Manhattan Plot",

pdf("~/courses/Ecological_Genomics_24/population_genomics/figures/ManhattanPlot_Gst_byRegions.pdf")
manhattan(vcf.div2.sub.MHplot,
          chr="V2",
          bp="POS",
          p="Gst",
          cex = 1, 
          cex.axis = 1,
          col = c("blue4", "orange3"), 
          suggestiveline = quantile(vcf.div2.sub.MHplot$Gst, 0.999), 
          chrlabs=as.character(c(1:8)),
          logp=F,
          ylab="Gst (among regions)",
          ylim=c(0,0.5))
dev.off()


outliers = vcf.div2.sub.MHplot %>%
  select(chr.main, POS, SNP, Gst) %>%
  filter(Gst > quantile(Gst, 0.999))

vcf.div2.sub.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(POS, value, color = name)) + 
  geom_point(alpha=0.5) +
  facet_wrap(vars(chr.main), scales = "free_y", nrow = 4, strip.position = "top") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(title="Genome-wide diversity", colour="Region", x="Chromosomal position (bp)", y="Gene diversity (Hs)")


vcf.div2.sub.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title="Genome-wide diversity", fill="Region", x="Gene diversity (Hs)",y="Frequency") +
  xlim(-0.01,0.2) + ylim(0,1500)


vcf.div2.sub.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  filter(value<0.5) %>%
  ggplot(aes(x = name, y=value, fill = name)) +
  geom_boxplot() +
  labs(x="Region", y="Gene diversity (Hs)",fill="Region")


vcf.div2.sub.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  summarise(avg=mean(value))


###################################################
################  In PopGenomeR ##################
###################################################

library(PopGenomeR)

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

genediv2 = neutrality.stats(genes)
plot(genediv2@SLIDE.POS,genediv2@Tajima.D)

plot(genediv2@theta_Watterson[is.na(genediv2@n.segregating.sites==0)==F]/genediv2@n.sites[is.na(genediv2@n.segregating.sites==0)==F])

