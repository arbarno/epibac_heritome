## Methylation analysis of Acropora spawning nanopore data ##

## Plotting CpG O/E from Acropora reference ----

# adapted from https://github.com/jldimond/Coral-CpG/blob/master/analyses/scripts/CpG_Density.R

acro_cpg <-read.delim("acropora_ref_cpg_bias.tsv", header=TRUE)

#Fitting mixture model with mixtools normalmixEM
library(mixtools) # v2.0.0.1
library(ggplot2) # v3.5.1
library(dplyr) # v1.1.4

cpg_data <- acro_cpg %>%
  filter(!grepl("^tRNA", gene)) %>%
  filter(cpg_bias >= 0.001 & cpg_bias <= 1.5) %>% #Cutting off high and low values (high value just cuts off 10 genes)
  select(gene, cpg_bias)
# 23946 genes passed filtration

cpg_data$cpg_bias <- as.numeric(cpg_data$cpg_bias)

set.seed(1)
mix.model <- normalmixEM(cpg_data$cpg_bias, k = 2)

mix.model[["mu"]] #means
## [1] 0.4369254 0.8418351
mix.model[["sigma"]] #standard deviations
## [1] 0.1351937 0.1688484
mix.model[["lambda"]] #amplitudes
## [1] 0.3252951 0.6747049

plot(mix.model, which = 2, col2 = c("red", "blue"), xlab2 = "CpG O/E", ylab2 = "Density", main2 = "")
abline(v = 0.588, lty = 2, col = "grey30", lwd = 1.5)

#Finds intersection point of two component model
intersect <- function(m1, s1, m2, s2, prop1, prop2){
  B <- (m1/s1^2 - m2/s2^2)
  A <- 0.5*(1/s2^2 - 1/s1^2)
  C <- 0.5*(m2^2/s2^2 - m1^2/s1^2) - log((s1/s2)*(prop2/prop1))
  (-B + c(1,-1)*sqrt(B^2 - 4*A*C))/(2*A)
}

model_intersect <- intersect(0.4369254, 0.1351937, 
                             0.8418351, 0.1688484, 
                             0.3252951, 0.6747049)
# Intersect is at 0.588

# get median CpG O/E value of genes
median(cpg_data$cpg_bias, na.rm = TRUE) # [1] 0.7372

# number of genes in high CpG O/E component
sum(cpg_data$cpg_bias > 0.588) # [1] 16138
# number of genes in low CpG O/E component
sum(cpg_data$cpg_bias < 0.588) # [1] 7808

#Test fit of single component model and compare
null.model <- MASS::fitdistr(cpg_data$CpG_bias, "normal")

loglik.null <- null.model$loglik
loglik.mix <- mix.model$loglik

## AIC (Akaike Information Criterion) + BIC (Bayesian Information Criterion)
# AIC and BIC for the null model
n <- length(cpg_data$CpG_bias)  # Number of data points
aic_null <- -2 * loglik.null + 2 * 2  # 2 parameters: mean and sd
bic_null <- -2 * loglik.null + log(n) * 2

# AIC and BIC for the mixture model
aic_mix <- -2 * loglik.mix + 2 * 5  # 5 parameters: mu1, mu2, sigma1, sigma2, lambda
bic_mix <- -2 * loglik.mix + log(n) * 5

# higher log-likelihood indicates better fit
loglik.null #[1] -523.2005
loglik.mix #[1] 261.18

# lower AIC and BIC values indicate better model performance
aic_null #[1] 1050.401
bic_null #[1] 1066.568
aic_mix #[1] -512.3599
bic_mix #[1] -471.9421
#

## Comparing methylation percents across of low/high CpG O/E ----

library(dplyr) # v1.1.4
library(ggplot2) # v3.5.1

# read in median methylation tsv file
mean.meths <- read.table("merged_mean_meths.tsv", sep = '\t', header = TRUE)
gene_universe <- read.delim("gene_universe.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# retain genes that passed filtering steps
mean.meths.f <- mean.meths %>%
  filter(!grepl("^tRNA", gene)) %>% # 23960 genes
  filter(meth_pos >= 3) %>%
  select(-meth_pos) %>%filter(gene %in% gene_universe) # 10883 genes

# Get aggregated percents for each gene (19831 total genes)
merged <- mean.meths.f %>%
  select(gene, merged)

# Read the list of genes from the first column of the tab-separated file
gene_list <- read.delim("cpg_all_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
high_genes <- read.delim("cpg_high_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
low_genes  <- read.delim("cpg_low_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# Filter just genes present in gene_list (23946 genes in gene list)
alls <- mean.meths.f %>%
  filter(gene %in% gene_list) %>%
  mutate(high_low = case_when(
    gene %in% high_genes ~ "high_cpg", # 16138 genes
    gene %in% low_genes  ~ "low_cpg", # 7808 genes
    TRUE ~ NA_character_  # if not present in either list
  ))

write.table(alls$gene, file = "filtered_gene_universe.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# alls results in 10877 genes
high <- alls %>%
  filter(high_low == 'high_cpg') # 8791 genes
write.table(high$gene, file = "filtered_highcpg_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
low <- alls %>%
  filter(high_low == 'low_cpg') # 2086 genes
write.table(low$gene, file = "filtered_lowcpg_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Plot boxplot of methylation of CpGs
box_plot <- ggplot(alls, aes(x = factor(high_low, levels = c("low_cpg", "high_cpg")), y = merged)) +
  geom_boxplot(aes(fill=high_low), colour= "black", lwd=0.2, fatten = 0.5, outliers = FALSE) +
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  labs(y="Methylation (%)") +
  scale_x_discrete(labels=c(expression("Low CpG"[O/E], "High CpG"[O/E]))) +
  scale_fill_manual(values = c("low_cpg" = "blue", "high_cpg" = "red")) +  
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12, family = "sans", colour = "black"),axis.title=element_text(size=12, family = "sans", colour = "black")) +
  theme(legend.position="none")
box_plot

ggsave("methylation_highlow_boxplot.pdf", plot = box_plot, width = 4, height = 3)

# Statistics comparing the high/low CpG O/E
hist(alls$merged) # not normal, many zeros
wilcox.test(merged ~ high_low, data = alls, conf.int = T) # W = 4658576, p-value < 2.2e-16
t.test(x = alls$merged) # t = 42.991, df = 10876, p-value < 2.2e-16
#

## topGO consider universe ----

# change this folder to point to your own "go_annot" folder
library(topGO) # v2.58.0

# remember to change the folder name to point to the folder containing your genes of interest lists
folder_of_interest = "topGO/"

# exclude files with "universe" in it
mult_files = grep(list.files(folder_of_interest), pattern="*genes.txt", value = T)
universe_file = grep(list.files(folder_of_interest), pattern="*universe.txt", value =T)

for (m in mult_files) {
  annot_filename = 'acropora_ref_goIDs.tsv'
  gene_id_to_go = readMappings(file=annot_filename)

  # shrink list of all GO terms down to the correct universe
  universe_genes = scan(paste0(folder_of_interest, universe_file), character(0), sep="\n")
  
  gene_id_to_go = gene_id_to_go[universe_genes]
  gene_id_to_go = gene_id_to_go[gene_id_to_go != 'no_hit']
  gene_names = names(gene_id_to_go)
  
  for (go_category in c('bp', 'cc', 'mf')) {
    print(paste("Current file:", m))
    genes_of_interest_filename = paste0(folder_of_interest, m)
    genes_of_interest = scan(genes_of_interest_filename, character(0), sep="\n")
    genes_of_interest <- genes_of_interest[genes_of_interest %in% names(gene_id_to_go)]
    
    genelist = factor(as.integer(gene_names %in% genes_of_interest))
    names(genelist) = gene_names
    
    GOdata = try(new("topGOdata", ontology=toupper(go_category), allGenes=genelist, gene2GO=gene_id_to_go, annotationFun=annFUN.gene2GO))
    
    # handle error
    if (class(GOdata) == "try-error") {
      print (paste0("Error for file", m, "!"))
      next
    }
    
    # weight01 is the default algorithm used in Alexa et al. (2006)
    weight01.fisher <- runTest(GOdata, statistic = "fisher")
    
    # generate a results table (for only the top 1000 GO terms)
    #   topNodes: highest 1000 GO terms shown
    #   numChar: truncates GO term descriptions at 1000 chars (basically, disables truncation)
    if (length(genes_of_interest) < 500) {
      results_table = GenTable(GOdata, P_value=weight01.fisher, orderBy="P_value", topNodes=100, numChar=1000)
    } else {
      results_table = GenTable(GOdata, P_value=weight01.fisher, orderBy="P_value", topNodes=250, numChar=1000)
    }
    
    # write it out into a file for python post-processing
    output_filename = paste0("topGO/", go_category, "_", m)
    write.table(results_table, file=output_filename, quote=FALSE, sep='\t')
  }
}
#

## NanoMethViz to visualize methylation levels across a ROI ----

# adapted from https://www.bioconductor.org/packages/release/bioc/vignettes/NanoMethViz/inst/doc/UsersGuide.html

library(rrvgo) # v1.20.0
library(org.Ce.eg.db) # v3.21.0

# Import the go terms data
l.go_analysis <- read.delim("bp_filtered_lowcpg_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Perform the similarity matrix (using C elegans, from H Putnam)
l.simMatrix <- calculateSimMatrix(l.go_analysis$GO.ID,
                                  orgdb="org.Ce.eg.db",
                                  ont="BP",
                                  method="Rel")


# Set scores, minus log-transform the pvalues first because for this package, higher is better
l.scores <- setNames(-log(l.go_analysis$P_value), l.go_analysis$GO.ID)

# Reduce GO terms to higher order terms, grouping a threshold and weighting by pvalue
# higher thresholds result in fewer groups
l.reducedTerms <- reduceSimMatrix(l.simMatrix,
                                  l.scores,
                                  threshold=0.90,
                                  orgdb="org.Ce.eg.db")

pdf("LowCpG_GOSemSim.pdf", width = 8, height = 4)
treemapPlot(l.reducedTerms)
dev.off()

# Import the go terms data
h.go_analysis <- read.delim("bp_filtered_highcpg_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Perform the similarity matrix (using C elegans, from H Putnam)
h.simMatrix <- calculateSimMatrix(h.go_analysis$GO.ID,
                                  orgdb="org.Ce.eg.db",
                                  ont="BP",
                                  method="Rel")


# Set scores, minus log-transform the pvalues first because for this package, higher is better
h.scores <- setNames(-log(h.go_analysis$P_value), h.go_analysis$GO.ID)

# Reduce GO terms to higher order terms, grouping a threshold and weighting by pvalue
# higher thresholds result in fewer groups
h.reducedTerms <- reduceSimMatrix(h.simMatrix,
                                  h.scores,
                                  threshold=0.90,
                                  orgdb="org.Ce.eg.db")

pdf("HighCpG_GOSemSim.pdf", width = 8, height = 4)
treemapPlot(h.reducedTerms)
dev.off()
#

## Boxplots of methylation levels for low & high CpG O/E ----

library(dplyr) # v1.1.4
library(ggplot2) # v3.5.1
library(tidyr) # v1.3.1
library(tibble) # v3.2.1
library(data.table) # v1.17.0

## Low CpG O/E genes

# read in percent methylation tsv file
filt_pcts <- fread("all_filt_pct_context.tsv.gz", sep = '\t', header = TRUE)
gene_list <- read.delim("cpg_all_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
high_genes <- read.delim("cpg_high_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
low_genes  <- read.delim("cpg_low_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# Apply filter functions to just the low CpG O/E genes
lows <- filt_pcts %>%
  filter(gene %in% low_genes) %>% # 35745 positions
  select(gene, nb22:sc5) %>%
  pivot_longer(-'gene', names_to = 'sample', values_to = 'meth_pct') %>%
  mutate(group = case_when(startsWith(sample, 'nb') ~ 'bmc_parent', startsWith(sample, 'nc') ~ 'placebo_parent', startsWith(sample, 'sb') ~ 'bmc_sperm', startsWith(sample, 'sc') ~ 'placebo_sperm')) %>%
  select(sample, group, everything() )%>%
  group_by(sample, group) %>%
  summarise(means = mean(meth_pct))

# Plot boxplot of methylation of low CpGs
lowcpg_plot <- ggplot(lows, aes(x = factor(group, levels = c("placebo_parent", "bmc_parent", "placebo_sperm", "bmc_sperm")), y = means)) +
  geom_boxplot(aes(fill=group), colour= "black", lwd=0.5) +
  labs(y="Methylation (%)") +
  scale_x_discrete(labels=c(expression("Placebo Parent", "Probiotic Parent", "Placebo Sperm", "Probiotic Sperm"))) +
  scale_fill_manual(values=c("blue", "blue", "blue", "blue")) +  
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black")) +
  theme(legend.position="none")
lowcpg_plot

ggsave("lowcpg_methylation_plot.pdf", plot = lowcpg_plot, width = 8, height = 6)

# stats
hist(lows$means)
shapiro.test(lows$means) # W = 0.7228, p-value = 0.00141
kruskal.test(lows$means, lows$group) 
#Kruskal-Wallis chi-squared = 9.359, df = 3, p-value = 0.02488
FSA::dunnTest(lows$means ~ lows$group, method="bh")
#                       Comparison          Z     P.unadj      P.adj
# 1         bmc_parent - bmc_sperm -2.6042372 0.009207901 0.05524741
# 2    bmc_parent - placebo_parent -1.0190493 0.308179547 0.36981546
# 3     bmc_sperm - placebo_parent  1.5851878 0.112923661 0.22584732
# 4     bmc_parent - placebo_sperm -2.4910095 0.012738072 0.03821422
# 5      bmc_sperm - placebo_sperm  0.1132277 0.909850033 0.90985003
# 6 placebo_parent - placebo_sperm -1.4719601 0.141031641 0.21154746

#--

## high CpG O/E genes

# Apply filter functions to just the high CpG O/E genes
highs <- filt_pcts %>%
  filter(gene %in% high_genes) %>% # 410127 positions
  select(gene, nb22:sc5) %>%
  pivot_longer(-'gene', names_to = 'sample', values_to = 'meth_pct') %>%
  mutate(group = case_when(startsWith(sample, 'nb') ~ 'bmc_parent', startsWith(sample, 'nc') ~ 'placebo_parent', startsWith(sample, 'sb') ~ 'bmc_sperm', startsWith(sample, 'sc') ~ 'placebo_sperm')) %>%
  select(sample, group, everything()) %>%
  group_by(sample, group) %>%
  summarise(means = mean(meth_pct))

# Plot boxplot of methylation of high CpGs
highcpg_plot <- ggplot(highs, aes(x = factor(group, levels = c("placebo_parent", "bmc_parent", "placebo_sperm", "bmc_sperm")), y = means)) +
  geom_boxplot(aes(fill=group), colour= "black", lwd=0.5) +
  labs(y="Methylation (%)") +
  scale_x_discrete(labels=c(expression("Placebo Parent", "Probiotic Parent", "Placebo Sperm", "Probiotic Sperm"))) +
  scale_fill_manual(values=c("red", "red", "red", "red")) +  
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black")) +
  theme(legend.position="none")
highcpg_plot

ggsave("highcpg_methylation_plot.pdf", plot = highcpg_plot, width = 8, height = 6)

# stats
hist(highs$means)
shapiro.test(highs$means) # W = 0.78344, p-value = 0.006111
kruskal.test(highs$means,highs$group) 
#Kruskal-Wallis chi-squared = 8.7436, df = 3, p-value = 0.0329
FSA::dunnTest(highs$means ~ highs$group, method="bh")
#                       Comparison          Z    P.unadj      P.adj
# 1         bmc_parent - bmc_sperm -2.4910095 0.01273807 0.07642843
# 2    bmc_parent - placebo_parent -0.5661385 0.57129962 0.68555955
# 3     bmc_sperm - placebo_parent  1.9248710 0.05424550 0.10849101
# 4     bmc_parent - placebo_sperm -2.1513264 0.03145045 0.09435135
# 5      bmc_sperm - placebo_sperm  0.3396831 0.73409518 0.73409518
# 6 placebo_parent - placebo_sperm -1.5851878 0.11292366 0.16938549

#--

## total genes methylation

# Apply filter functions to just the CpG O/E genes
tots <- filt_pcts %>%
  filter(gene %in% gene_list) %>% # 445870 positions of a total 446080 positions
  select(gene, nb22:sc5) %>%
  pivot_longer(-'gene', names_to = 'sample', values_to = 'meth_pct') %>%
  mutate(group = case_when(startsWith(sample, 'nb') ~ 'bmc_parent', startsWith(sample, 'nc') ~ 'placebo_parent', startsWith(sample, 'sb') ~ 'bmc_sperm', startsWith(sample, 'sc') ~ 'placebo_sperm')) %>%
  select(sample, group, everything()) %>%
  group_by(sample, group) %>%
  summarise(means = mean(meth_pct))

# Plot boxplot of methylation of all genes
meth_plot <- ggplot(tots, aes(x = factor(group, levels = c("placebo_parent", "bmc_parent", "placebo_sperm", "bmc_sperm")), y = means)) +
  geom_boxplot(aes(fill=group), colour= "black", lwd=0.5) +
  labs(y="Methylation (%)") +
  scale_x_discrete(labels=c(expression("Placebo Parent", "Probiotic Parent", "Placebo Sperm", "Probiotic Sperm"))) +
  scale_fill_manual(values=c("grey", "grey", "grey", "grey")) +  
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor=element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black")) +
  theme(legend.position="none")
meth_plot

ggsave("methylation_plot.pdf", plot = meth_plot, width = 8, height = 6)

# stats
hist(tots$means)
shapiro.test(tots$means) # W = 0.75136, p-value = 0.002761
kruskal.test(tots$means,tots$group) 
#Kruskal-Wallis chi-squared = 9.359, df = 3, p-value = 0.02488
FSA::dunnTest(tots$means ~ tots$group, method="bh")
#                                 Comparison          Z     P.unadj      P.adj
# 1                   bmc_parent - bmc_sperm -2.6042372 0.009207901 0.05524741
# 2              bmc_parent - placebo_parent -1.0190493 0.308179547 0.36981546
# 3               bmc_sperm - placebo_parent  1.5851878 0.112923661 0.22584732
# 4               bmc_parent - placebo_sperm -2.4910095 0.012738072 0.03821422
# 5                bmc_sperm - placebo_sperm  0.1132277 0.909850033 0.90985003
# 6           placebo_parent - placebo_sperm -1.4719601 0.141031641 0.21154746
#

## Genome context of methylation ----

library(dplyr) # v1.1.4
library(data.table) # v1.17.0

df<-fread("merged_bed_annotated.bed.gz", sep = '\t', header=FALSE, fill = 14)
head(df)

df$V5 <- as.numeric(df$V5) # CpG methylated
df$V6 <- as.numeric(df$V6) # unmethylated

total_meth <- sum(df$V5, na.rm = TRUE) # 43546365
total_unmeth <- sum(df$V6, na.rm = TRUE) # 593856473
# total percent = 6.83%

intergenic <- df %>%
  filter(V7 == "intergenic") %>%
  mutate(intergenic_meth = sum(V5), na.rm = TRUE) %>% # intergenic_meth = 34717180
  mutate(intergenic_unmeth = sum(V6), na.rm = TRUE) %>% # intergenic_unmeth = 468446942
  select(V7, intergenic_meth, intergenic_unmeth) %>%
  slice(1)
# intergenic percent = 6.90%

exons <- df %>%
  filter(startsWith(V11, "Exon")) %>%
  mutate(exon_meth = sum(V5), na.rm = TRUE) %>% # inter_meth = 4678535
  mutate(exon_unmeth = sum(V6), na.rm = TRUE) %>% # inter_unmeth = 51882513
  select(V11, exon_meth, exon_unmeth) %>%
  slice(1)
# exon percent = 8.27%

introns <- df %>%
  filter(startsWith(V11, "Intron")) %>%
  mutate(intron_meth = sum(V5), na.rm = TRUE) %>% # intron_meth = 3994722
  mutate(intron_unmeth = sum(V6), na.rm = TRUE) %>% # intron_unmeth = 69751983
  select(V11, intron_meth, intron_unmeth) %>%
  slice(1)
# intron percent = 5.73%
#

## Genome context of samples ----

library(dplyr) # v1.1.4
library(tidyr) # v1.3.1
library(tibble) # v3.2.1
library(ggplot2) # v3.5.1

# read in percent methylation tsv file
samp_context <- read.table("sample_genome_context.tsv", sep = '\t', header = TRUE)

context_data <- samp_context %>%
  pivot_longer(-c(sample, generation, inoculation, group), names_to = 'context', values_to = 'percents')

# Reorder levels
context$group <- factor(context$group, levels = c("placebo_parent", "bmc_parent", "placebo_sperm", "bmc_sperm"))
context$context <- factor(context$context, levels = c('Exonic', 'Intronic', 'Intergenic', 'TotalCpGs'))

# plot the boxplot
p<-ggplot(context, aes(x= context, y=percents))+
  geom_boxplot(aes(fill = group), colour= "black", linewidth=0.3, fatten = 0.3) +
  scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40"), labels = c("Placebo Parent", "BMC Parent", "Placebo Sperm", "BMC Sperm"), name = NULL) +
  geom_point(aes(group = group), size=0.6, position = position_dodge(width=0.75)) +
  facet_wrap(vars(context), scales= "free_x", nrow = 1, ncol = 4) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 0, vjust = 0, hjust = 0.5),
        strip.text = element_blank(), strip.background = element_blank(),
        legend.position.inside = c(0.95, 0.95), legend.justification = c("right", "top"),
        axis.text=element_text(size=12, family = "sans", colour = "black"),
        axis.title=element_text(size=12, family = "sans", colour = "black"),
        legend.text = element_text(size=12, family = "sans", colour = "black"),
        legend.title = element_text(size=12, family = "sans", colour = "black")) +
  labs(y="Methylation (%)", x="Genomic Context")
p

ggsave("sample_genomic_context_boxplots.pdf", plot = p, width = 8, height = 3)

# stats
exons <- context_data %>%
  filter(context == 'Exonic')
hist(exons$percents)
shapiro.test(exons$percents) # W = 0.96177, p-value = 0.8087
t<-aov(percents ~ generation + inoculation, data = exons)
summary(t)
#                       Df Sum Sq Mean Sq F value Pr(>F)  
# generation              1  0.069   0.069   0.068 0.8008  
# inoculation             1  4.175   4.175   4.123 0.0768 .
# generation:inoculation  1  1.043   1.043   1.030 0.3398  
# Residuals               8  8.101   1.013 

introns <- context_data %>%
  filter(context == 'Intronic')
hist(introns$percents)
shapiro.test(introns$percents) # W = 0.9026, p-value = 0.1714
t<-aov(percents ~ generation * inoculation, data = introns)
summary(t)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# generation              1 0.1162 0.11616   1.383  0.273
# inoculation             1 0.0653 0.06530   0.777  0.404
# generation:inoculation  1 0.0307 0.03066   0.365  0.563
# Residuals               8 0.6721 0.08401

inter <- context_data %>%
  filter(context == 'Intergenic')
hist(inter$percents)
shapiro.test(inter$percents) # W = 0.77079, p-value = 0.004444
rcompanion::scheirerRayHare(percents ~ generation * inoculation, data = inter) 
#                       Df  Sum Sq      H p.value
# generation              1 108.000 8.3077 0.00395
# inoculation             1   3.000 0.2308 0.63095
# generation:inoculation  1   5.333 0.4103 0.52184
# Residuals               8  26.667  
wilcox.test(percents ~ generation, data = inter)
#W = 0, p-value = 0.002165

totals <- context %>%
  filter(context == 'TotalCpGs')
hist(totals$percents)
shapiro.test(totals$percents) # W = 0.83045, p-value = 0.02123
rcompanion::scheirerRayHare(percents ~ generation * inoculation, data = totals) 
#                       Df  Sum Sq      H p.value
# generation              1 108.000 8.3077 0.00395
# inoculation             1   5.333 0.4103 0.52184
# generation:inoculation  1   3.000 0.2308 0.63095
# Residuals               8  26.667 
wilcox.test(percents ~ generation, data = totals)
#W = 0, p-value = 0.002165
#

## Correlation matrix ----

library(dplyr) # v1.1.4
library(tibble) # v3.2.1
library(data.table) # v1.17.0
library(pheatmap) # v1.0.12
library(ggplot2) # v3.5.1

# read in mean meths tsv fil to run glms
filt_pcts <- fread("all_filt_pct_context.tsv.gz", sep = '\t', header = TRUE)

# Filter out the data for the correlation plot
df <- filt_pcts %>%
  mutate(scaf_pos = paste(.$scaffold, .$pos, sep = "_")) %>%
  select(scaf_pos, nb22:sc5) %>%
  column_to_rownames('scaf_pos')

hist(filt_pcts$merged) # not normal, many zeros
corrk <- pcaPP::cor.fk(df)

# set the self-correlations as NA
diag(corrk) <- NA

# read in kendall tau table
taus <- fread("kendall_tau_table.tsv", sep = '\t', header = TRUE)

taus.f <- taus %>%
  filter(group != "sperm")
hist(taus.f$coefficient)
shapiro.test(taus.f$coefficient) # W = 0.85214, p-value = 0.004625
wilcox.test(coefficient ~ group, data = taus.f, conf.int = T) # W = 90, p-value = 3.686e-05
t.test(x = taus.f$coefficient) # t = 195.86, df = 20, p-value < 2.2e-16

# plot the correlation as a heatmap
pheatmap(corrk,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(20)),
         breaks = seq(0.60, 0.8, length.out = 20 + 1),
         clustering_method = "complete",
         display_numbers = T,
         legend_breaks   = c(0.6, 0.7, 0.8),
         fontsize = 12,
         fontsize_row = 10,
         fontsize_col = 10,
         filename = "sample_correlation_plot.pdf",
         width = 6,
         height = 6
         )
#

## PCA analysis of methylation ----

library(dplyr) # v1.1.4
library(data.table) # v1.17.0
library(FactoMineR) # v2.11
library(factoextra) # v1.0.7

# read in percent methylation tsv file
filt_pcts <- fread("all_filt_pct_context.tsv.gz", sep = '\t', header = TRUE)
gene_list <- read.delim("cpg_all_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
high_genes <- read.delim("cpg_high_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
low_genes  <- read.delim("cpg_low_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# Apply filter functions to just the CpG O/E genes
tots <- filt_pcts %>%
  mutate(scaf_pos = paste(.$scaffold, .$pos, sep = "_")) %>%
  select(scaf_pos, nb22:sc5) %>%
  column_to_rownames('scaf_pos') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  mutate(group = case_when(startsWith(sample, 'nb') ~ 'bmc_parent', startsWith(sample, 'nc') ~ 'placebo_parent', startsWith(sample, 'sb') ~ 'bmc_sperm', startsWith(sample, 'sc') ~ 'placebo_sperm')) %>%
  select(sample, group, everything())

# Run the PCA
res.pca <- PCA(tots, quali.sup=c(1,2), scale = TRUE, graph = FALSE)

fviz_pca_ind(res.pca,  geom = "point", pointsize = 5,
             habillage="group", label = 'none',
             palette = c("#5D3A9B","#5D3A9B", "#E66100", "#E66100"),
             mean.point = FALSE,
             invisible=c("quali.sup"),
             alpha.ind = 1,
             addEllipses = FALSE) +
  scale_shape_manual(values=c(19,17,19,17)) +
  ggtitle("Total") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),
        axis.title=element_text(size=14, family = "sans", colour = "black")) +
  theme(legend.text = element_text(size=12, family = "sans", colour = "black"), 
        legend.title = element_text(size=14, family = "sans", colour = "black"))
  
ggsave("PCA_scaf_positions_all.pdf", width=8, height=6)
#


## Normality testing ----

library(dplyr) # v1.1.4

# read in median methylation tsv file
mean.meths <- read.table("all_mean_meths.tsv", sep = '\t', header = TRUE)
gene_universe <- read.delim("gene_universe.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# reformat the mean methylation data to run the glms
mean.meths.f <- mean.meths %>%
  filter(!grepl("^tRNA", gene)) %>% # 23960 genes
  select(-meth_pos, -merged) %>%
  filter(gene %in% gene_universe) %>% # 12492 genes
  column_to_rownames('gene') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'adult', startsWith(sample, 's') ~ 'sperm'),
         inoculation = case_when(startsWith(sample, 'nb') | startsWith(sample, 'sb') ~ 'bmc',
                                 startsWith(sample, 'nc') | startsWith(sample, 'sc') ~ 'placebo')) %>%
  select(sample, generation, inoculation, everything()) %>%
  pivot_longer(-c(sample, generation, inoculation), names_to = 'gene', values_to = 'mean_meth')
  
mean.meths.f$generation <- as.factor(mean.meths.f$generation)
mean.meths.f$inoculation <- as.factor(mean.meths.f$inoculation)

# perform the shapiro test for all the conditions and count how many are not normally distributed
lev.gen <- mean.meths.f %>%
  group_by(gene) %>%
  levene_test(mean_meth ~ generation) %>%
  ungroup()  %>% 
  mutate(bh = p.adjust(p, method = 'BH')) %>%
  filter(bh < 0.05) # 209 genes of 12492 total are not-normal

lev.inoc <- mean.meths.f %>%
  group_by(gene) %>%
  levene_test(mean_meth ~ inoculation) %>%
  ungroup()  %>% 
  mutate(bh = p.adjust(p, method = 'BH')) %>%
  filter(bh < 0.05) # 44 genes of 12492 total are not-normal

shap.adult <- mean.meths.f %>%
  filter(grepl('adult', generation)) %>%
  group_by(gene)  %>%
  filter(n_distinct(mean_meth) > 1) %>%
  do(tidy(shapiro.test(.$mean_meth)))  %>% 
  ungroup()  %>%
  mutate(bh = p.adjust(p.value, method = 'BH')) %>%
  filter(bh < 0.05) # 3221 genes of 12492 total are not-normal

shap.sperm <- mean.meths.f %>%
  filter(grepl('sperm', generation)) %>%
  group_by(gene)  %>%
  filter(n_distinct(mean_meth) > 1) %>%
  do(tidy(shapiro.test(.$mean_meth)))  %>% 
  ungroup()  %>%
  mutate(bh = p.adjust(p.value, method = 'BH')) %>%
  filter(bh < 0.05) # 4730 genes of 12492 total are not-normal

shap.placebo <- mean.meths.f %>%
  filter(grepl('placebo', inoculation)) %>%
  group_by(gene)  %>%
  filter(n_distinct(mean_meth) > 1) %>%
  do(tidy(shapiro.test(.$mean_meth)))  %>% 
  ungroup()  %>%
  mutate(bh = p.adjust(p.value, method = 'BH')) %>%
  filter(bh < 0.05) # 3744 genes of 12492 total are not-normal

shap.bmc <- mean.meths.f %>%
  filter(grepl('bmc', inoculation)) %>%
  group_by(gene)  %>%
  filter(n_distinct(mean_meth) > 1) %>%
  do(tidy(shapiro.test(.$mean_meth))) %>%
  ungroup()  %>%
  mutate(bh = p.adjust(p.value, method = 'BH')) %>%
  filter(bh < 0.05) # 3768 genes of 12492 total are not-normal
#

## GLM ----

library(dplyr) # v1.1.4
library(tibble) # v3.2.1
library(tidyr) # v1.3.1
library(purrr) # v1.0.4
library(broom) # v1.0.7

# read in mean meths tsv fil to run glms
mean.meths <- read.table("all_mean_meths.tsv", sep = '\t', header = TRUE)
gene_universe <- read.delim("gene_universe.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# reformat the mean methylation data to run the glms
mean.meths.f <- mean.meths %>%
  filter(!grepl("^tRNA", gene)) %>% # 23960 genes
  filter(meth_pos >= 3) %>%
  select(-meth_pos) %>%
  filter(gene %in% gene_universe) %>% # 10883 genes
  column_to_rownames('gene') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'adult', startsWith(sample, 's') ~ 'sperm'),
         inoculation = case_when(startsWith(sample, 'nb') | startsWith(sample, 'sb') ~ 'bmc',
                                 startsWith(sample, 'nc') | startsWith(sample, 'sc') ~ 'placebo')) %>%
  select(sample, generation, inoculation, everything()) %>%
  pivot_longer(-c(sample, generation, inoculation), names_to = 'gene', values_to = 'mean_meth')

mean.meths.f$generation <- factor(mean.meths.f$generation, levels = c('adult', 'sperm'))
mean.meths.f$inoculation <- factor(mean.meths.f$inoculation, levels = c('placebo', 'bmc'))

# run the glm on each gene
hope <- mean.meths.f %>%
  nest(data = -gene) %>%
  mutate(
    model = map(data, ~ {
      l <- glm(mean_meth ~ generation + inoculation + generation:inoculation, data = ., family = gaussian())
      s <- tryCatch(step(l, trace = 0), error = function(e) e)
      if (inherits(s, "error") || !isTRUE(s$converged)) {
        tibble(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA)
      } else {
        tidy(s)
      }
    })
  ) %>%
  unnest_wider(model) %>%
  select(-data)

# Filter the results to keep only the relevant comparisons
res <- hope %>%
  unnest(c(term, estimate, std.error, statistic, p.value)) %>%
  filter(!if_all(2:ncol(.), is.na))

# Pivot the table for better readability
res_wide <- pivot_wider(res, names_from = term, values_from = c(estimate, std.error, statistic, p.value))
res_wide <- arrange(res_wide, gene)

# Adjust for multiple testing (e.g., Benjamini-Hochberg) for each term
res_gen <- res_wide %>%
  filter(!is.na(`p.value_generationsperm`)) %>%
  mutate(bh = p.adjust(`p.value_generationsperm`, method = 'BH')) %>%
  filter(bh < 0.05)

write.table(res_gen$gene, "glm_gen_converged_bh.txt", sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
write.table(res_gen, "glm_gen.tsv", sep="\t", row.names=FALSE)

#--

res_bmc <- res_wide %>%
  filter(!is.na(`p.value_inoculationbmc`)) %>%
  mutate(bh = p.adjust(`p.value_inoculationbmc`, method = 'BH')) %>%
  filter(bh < 0.05)

write.table(res_bmc$gene, "glm_bmc_converged_bh.txt", sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
write.table(res_bmc, "glm_bmc.tsv", sep="\t", row.names=FALSE)

#--

res_interact <- res_wide %>%
  filter(!is.na(`p.value_generationsperm:inoculationbmc`)) %>%
  mutate(bh = p.adjust(`p.value_generationsperm:inoculationbmc`, method = 'BH')) %>%
  filter(bh < 0.05)

write.table(res_interact, "glm_gen_bmc_converged_bh.txt", sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
write.table(res_interact, "glm_interact.tsv", sep="\t", row.names=FALSE)
#

## CpG of GLMs ----

library(dplyr) # v1.1.4
library(purrr) # v1.0.4
library(ggplot2) # v3.5.1

# read in mean meths tsv fil to run glms
acro_cpg <- read.delim("acropora_ref_cpg_bias.tsv", header=TRUE)
glm_gen <- read.delim("glm_gen_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm_bmc <- read.delim("glm_bmc_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# filter the cpg data to only include the relevant genes
df1 <- acro_cpg %>%
  filter(!grepl("^tRNA", gene)) %>%
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>%
  select(gene, CpG.bias)
df1$CpG.bias <- as.numeric(df1$CpG.bias)

df2 <- acro_cpg %>%
  filter(!grepl("^tRNA", gene)) %>%
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>%
  filter(gene %in% glm_gen) %>%
  select(gene, CpG.bias)
df2$CpG.bias <- as.numeric(df2$CpG.bias)

df3 <- acro_cpg %>%
  filter(!grepl("^tRNA", gene)) %>%
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>%
  filter(gene %in% glm_bmc) %>%
  select(gene, CpG.bias)
df3$CpG.bias <- as.numeric(df3$CpG.bias)

# make list of dataframes to overlay
my_list <- list(CpGbias = df1, glm_gen = df2, glm_bmc = df3)

# stack them all into a tibble
long_df <- imap_dfr(my_list, ~ tibble(CpG.bias = .x$CpG.bias, source = .y))
long_df$source <- factor(long_df$source, levels = c('glm_bmc', 'glm_gen', 'CpGbias'))

hist(long_df$CpG.bias) # not normal, many zeros
t<-aov(long_df$CpG.bias ~ long_df$source)
summary(t)
TukeyHSD(t)

# plot
ridge <- ggplot(long_df, aes(x = CpG.bias, y = source, fill = source)) +
  ggridges::geom_density_ridges(scale = 1.5, rel_min_height = 0.001, alpha = 0.5) +
  scale_fill_manual(values=c("#7648E3", "#2A6D61", "#9C9C9C")) +
  scale_x_continuous(breaks = seq(0, 1.5, by = 0.5), limits = c(0, 1.5)) +
  geom_vline(xintercept = 0.588, color = "grey30", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
        axis.text = element_text(size=14, family = "sans", colour = "black"),
        axis.ticks.x = element_line(colour="black"),
        axis.title = element_text(size=14, family = "sans", colour = "black"),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
        legend.position = 'none') +
  scale_y_discrete(labels = c("GLM inoculation", "GLM generation", "Total genes"),
                   expand = expansion(add = c(0.2)),
                   position = 'right') +
  labs(x = expression("CpG"[O/E]), y = NULL)

ggsave("cpg_genes_ridgeplot.pdf", plot = ridge, width=8, height=6)
#

## Differential abundances of genes ----

library(dplyr) # v1.1.4
library(tibble) # v3.2.1
library(tidyr) # v1.3.1
library(broom) # v1.0.7
library(pheatmap) # v1.0.12
library(ggplot2) # v3.5.1

# read in mean meths tsv fil to run glms
mean.meths <- read.table("all_mean_meths.tsv", sep = '\t', header = TRUE)
gene_universe <- read.delim("gene_universe.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm_gen <- read.delim("glm_gen_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm_bmc <- read.delim("glm_bmc_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# reformat the mean methylation data to run the glms
mean.meths.f <- mean.meths %>%
  filter(!grepl("^tRNA", gene)) %>% # 23960 genes
  filter(meth_pos >= 3) %>%
  select(-meth_pos, -merged, -sperm, -parent) %>%
  filter(gene %in% gene_universe) %>% # 10883 genes
  column_to_rownames('gene') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'adult', startsWith(sample, 's') ~ 'sperm'),
         inoculation = case_when(startsWith(sample, 'nb') | startsWith(sample, 'sb') ~ 'bmc',
                                 startsWith(sample, 'nc') | startsWith(sample, 'sc') ~ 'placebo')) %>%
  select(sample, generation, inoculation, everything()) %>%
  pivot_longer(-c(sample, generation, inoculation), names_to = 'gene', values_to = 'mean_meth')

mean.meths.f$generation <- factor(mean.meths.f$generation, levels = c('adult', 'sperm'))
mean.meths.f$inoculation <- factor(mean.meths.f$inoculation, levels = c('placebo', 'bmc'))

# run the heatmap for the glm_gen and glm_bmc
gen_meths <- mean.meths.f %>% 
  filter(gene %in% glm_gen) %>%
  group_by(gene, generation) %>%
  summarise(avg = mean(mean_meth)) %>%
  mutate(diff = max(avg) - min(avg)) %>%
  filter(diff >= 10) %>%
  distinct(gene) %>%
  inner_join(., mean.meths.f, by ='gene') %>%
  select(sample, generation, inoculation, gene, mean_meth) %>%
  pivot_wider(id_cols = c('sample', 'generation', 'inoculation'), names_from = 'gene', values_from = 'mean_meth') %>%
  column_to_rownames('sample')

# plot the heatmap
heatr <- t(gen_meths[,3:ncol(gen_meths)]) # just the data, adjust the group as needed
heatr <- t(scale(t(heatr))) # scale before pheatmap
heatr.gs <- select(gen_meths, generation, inoculation) # groups

heatr.gs$generation <- factor(heatr.gs$generation, levels = c('adult', 'sperm'))
heatr.gs$inoculation <- factor(heatr.gs$inoculation, levels = c('placebo', 'bmc'))

pheatmap(heatr,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, 'RdYlBu'))(256)),
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2', 
         cutree_rows = 1, cutree_cols = 2,
         show_rownames = FALSE,
         angle_col = 45,
         border_color= NA,
         annotation_col = heatr.gs,
         filename = "/Users/barnoar/Documents/acropora_spawning_23/figures/heatmap_gen_glm_10.pdf",
         width = 6,
         height = 6
         )

#--

# run the heatmap for the glm_bmc
bmc_meths <- mean.meths.f %>% 
  filter(gene %in% glm_bmc) %>%
  group_by(gene, inoculation) %>%
  summarise(avg = mean(mean_meth)) %>%
  mutate(diff = max(avg) - min(avg)) %>%
  filter(diff >= 10) %>%
  distinct(gene) %>%
  inner_join(., mean.meths.f, by ='gene') %>%
  select(sample, generation, inoculation, gene, mean_meth) %>%
  pivot_wider(id_cols = c('sample', 'generation', 'inoculation'), names_from = 'gene', values_from = 'mean_meth') %>%
  column_to_rownames('sample')

# plot the heatmap
heatr <- t(bmc_meths[,3:ncol(bmc_meths)]) # just the data, adjust the group as needed
heatr <- t(scale(t(heatr))) # scale before pheatmap
heatr.gs <- select(bmc_meths, generation, inoculation) # groups

heatr.gs$generation <- factor(heatr.gs$generation, levels = c('adult', 'sperm'))
heatr.gs$inoculation <- factor(heatr.gs$inoculation, levels = c('placebo', 'bmc'))

pheatmap(heatr,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, 'RdYlBu'))(256)),
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2', 
         cutree_rows = 1, cutree_cols = 2,
         show_rownames = FALSE,
         angle_col = 45,
         border_color= NA,
         annotation_col = heatr.gs,
         filename = "/Users/barnoar/Documents/acropora_spawning_23/figures/heatmap_bmc_glm_10.pdf",
         width = 6,
         height = 3
         )
#


## NanoMethViz to visualize methylation levels across a ROI ----

# https://www.bioconductor.org/packages/release/bioc/vignettes/NanoMethViz/inst/doc/UsersGuide.html

library(NanoMethViz) # v3.4.0
library(dplyr) # v1.1.4
library(data.table) # v1.17.0
library(tidyr) # v1.3.1
library(Rsamtools) # v2.25.0
library(patchwork) # v1.3.0

# Import the gff data to fit the format of NanoMethViz
anno <- rtracklayer::import("acropora_ref_assembly.gff.gz")

genes <- anno %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  dplyr::rename(chr = "seqnames", gene_id = "ID") %>%
  select("gene_id", "chr", "strand", "start", "end") %>%
  filter(!grepl("^tRNA", gene_id)) %>%
  mutate("transcript_id" = gene_id, "symbol" = gene_id)

# Define the samples annotations
sample_anno <- read.table("sample_anno.tsv", sep = '\t', header = TRUE)

# Build the NanoMethResult
nmr <- NanoMethResult(
  methy = "combined_data.sorted.tsv.gz",
  samples = sample_anno,
  exons = genes
)

# Plot gene g2113, but put the window prop further left to include g2112
gs <- plot_gene(nmr, "g2113", window_prop = c(1.655,0.12), smoothing_window = 500) #
plots <- gs$patches$plots
top_plot <- plots[[1]] + 
  scale_y_continuous(limits = c(0.5, 0.75))
gs_plot <- top_plot / plots[[2]]
print(gs_plot)
#
