## Methylation analysis of Acropora spawning nanopore data ##

## Create universal gene universe file with minimum CpG >= 3 in all samples ----

# read in files with gene methylation data
gene.counts <- read.table("gene_level_counts.tsv", sep = '\t', header = TRUE)

library(dplyr) # v1.2.1

filtered.genes <- gene.counts %>%
  group_by(gene) %>%
  summarise(min_n_cpgs = min(n_cpgs), n_samples_present = n(), .groups = "drop") %>%
  filter(min_n_cpgs >= 3, n_samples_present == 12, !grepl("intergenic|no_info", gene)) %>%
  distinct(gene)

write.table(filtered.genes$gene, file = "gene_universe_conserved_3CpGs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 10763 genes

## Plotting CpG O/E from Acropora reference ----

# adapted from https://github.com/jldimond/Coral-CpG/blob/master/analyses/scripts/CpG_Density.R

acro_cpg <-read.delim("acropora_ref_cpg_bias.tsv", header=TRUE)

#Fitting mixture model with mixtools normalmixEM
library(mixtools) # v2.0.0.1
library(ggplot2) # v4.0.3
library(dplyr) # v1.2.1

cpg_data <- acro_cpg %>%
  filter(!grepl("^tRNA", gene)) %>%
  filter(cpg_bias >= 0.001 & cpg_bias <= 1.5) %>% #Cutting off high and low values (high value just cuts off 8 genes)
  select(gene, cpg_bias)
# 23946 genes passed filtration

cpg_data$cpg_bias <- as.numeric(cpg_data$cpg_bias)

set.seed(1)
mix.model <- normalmixEM(cpg_data$cpg_bias, k = 2)
# number of iterations= 232

mix.model[["mu"]] #means
## [1] 0.4369254 0.8418351
mix.model[["sigma"]] #standard deviations
## [1] 0.1351937 0.1688484
mix.model[["lambda"]] #amplitudes
## [1] 0.3252951 0.6747049

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

plot(mix.model, which = 2, col2 = c("red", "blue"), xlab2 = "CpG O/E", ylab2 = "Density", main2 = "")
abline(v = 0.588, lty = 2, col = "grey30", lwd = 1.5)

# get median CpG O/E value of genes
median(cpg_data$cpg_bias, na.rm = TRUE) # [1] 0.7372
# number of genes in high CpG O/E component
sum(cpg_data$cpg_bias > 0.588) # [1] 16138
high <- cpg_data %>%
  filter(cpg_bias > 0.588)
write.table(high$gene, file = "highcpg_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# number of genes in low CpG O/E component
sum(cpg_data$cpg_bias < 0.588) # [1] 7808
low <- cpg_data %>%
  filter(cpg_bias > 0.588)
write.table(low$gene, file = "lowcpg_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Test fit of single component model and compare
null.model <- MASS::fitdistr(cpg_data$cpg_bias, "normal")

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

library(dplyr) # v1.2.1
library(ggplot2) # v4.0.3

# read in weighted mean methylation tsv file
weighted.avgs <- read.table("gene_level_counts.tsv", sep = '\t', header = TRUE)
gene.universe <- read.delim("gene_universe_conserved_3CpGs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# retain genes that passed filtering steps
merged <- weighted.avgs %>%
  filter(gene %in% gene.universe) %>% # 10763 genes passing threshold
  group_by(gene) %>%
  summarize(gene_avg = (sum(methylated_count)/sum(total_count) * 100), .groups = "drop") # adds all methylated counts/total counts

# Read the list of genes from the first column of the tab-separated file
high <- read.delim("highcpg_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
low  <- read.delim("lowcpg_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# Filter just genes present in gene_list (23946 genes in gene list)
alls <- merged %>%
  filter(gene %in% gene.universe) %>%
  mutate(high_low = case_when(
    gene %in% high ~ "high_cpg", # 8760 genes
    gene %in% low  ~ "low_cpg", # 1997 genes
    TRUE ~ NA_character_  # if not present in either list
  )) %>%
  filter(!is.na(high_low)) # only 6 genes removed
table(alls$high_low) # to get the number of genes in each group

# Plot boxplot of methylation of high/low CpGs
box_plot <- ggplot(alls, aes(x = factor(high_low, levels = c("low_cpg", "high_cpg")), y = gene_avg)) +
  geom_boxplot(aes(fill = high_low), colour = "black", linewidth = 0.2, fatten = 0.5, outlier.shape = NA) +
  theme_classic() +
  labs(x = NULL, y = "Methylation (%)") +
  scale_x_discrete(labels = c(expression("Low CpG"[O/E]), expression("High CpG"[O/E]))) +
  scale_fill_manual(values = c("low_cpg" = "blue", "high_cpg" = "red")) +  
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text = element_text(size = 12, family = "sans", colour = "black"), axis.title = element_text(size = 12, family = "sans", colour = "black")) +
  theme(legend.position="none")
box_plot

ggsave("methylation_highlow_boxplot.pdf", plot = box_plot, width = 4, height = 3)

# Statistics comparing the high/low CpG O/E
hist(alls$gene_avg) # not normal, many zeros
wilcox.test(gene_avg ~ high_low, data = alls, conf.int = T) # W = 4635151, p-value < 2.2e-16
#

## topGO consider universe ----

# change this folder to point to your own "go_annot" folder
library(topGO) # v2.62.0

# remember to change the folder name to point to the folder containing your genes of interest lists
folder_of_interest = "cpgs/"

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
    output_filename = paste0("cpgs/", go_category, "_", m)
    write.table(results_table, file=output_filename, quote=FALSE, sep='\t')
  }
}
#

## Reduce and visualize GO terms with rrvgo ----

# https://github.com/ssayols/rrvgo
# http://revigo.irb.hr/

library(rrvgo) # v1.20.0
library(org.Ce.eg.db) # v3.21.0

# Import the go terms data
l.go_analysis <- read.delim("bp_lowcpg_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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
h.go_analysis <- read.delim("bp_highcpg_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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

## Genome context of methylation ----

library(dplyr) # v1.2.1
library(data.table) # v1.18.4

bed <-fread("merged_clair3_bed_annotated.bed", sep = '\t', header=FALSE, fill = 14)
head(bed)

bed$V5 <- as.numeric(bed$V5) # CpG methylated
bed$V6 <- as.numeric(bed$V6) # unmethylated

total_meth <- sum(bed$V5, na.rm = TRUE) # 39097652
total_unmeth <- sum(bed$V6, na.rm = TRUE) # 584146003
# total percent = 6.27%

intergenic <- bed %>%
  filter(V7 == "intergenic") %>%
  mutate(intergenic_meth = sum(V5), na.rm = TRUE) %>% # intergenic_meth = 31205812
  mutate(intergenic_unmeth = sum(V6), na.rm = TRUE) %>% # intergenic_unmeth = 459285313
  dplyr::select(V7, intergenic_meth, intergenic_unmeth) %>%
  .[!duplicated(.), ]
# intergenic percent = 6.36%

exons <- bed %>%
  filter(startsWith(V11, "Exon")) %>%
  mutate(exon_meth = sum(V5), na.rm = TRUE) %>% # exon_meth = 4236771
  mutate(exon_unmeth = sum(V6), na.rm = TRUE) %>% # exon_unmeth = 53141674
  dplyr::select(V11, exon_meth, exon_unmeth) %>%
  .[!duplicated(.$exon_meth), ]
# exon percent = 7.38%

introns <- bed %>%
  filter(startsWith(V11, "Intron")) %>%
  mutate(intron_meth = sum(V5), na.rm = TRUE) %>% # intron_meth = 3512946
  mutate(intron_unmeth = sum(V6), na.rm = TRUE) %>% # intron_unmeth = 68063739
  dplyr::select(V11, intron_meth, intron_unmeth) %>%
  .[!duplicated(.$intron_meth), ]
# intron percent = 4.91%
#

## Genome context of samples ----

library(dplyr) # v1.2.1
library(tidyr) # v1.3.1
library(tibble) # v3.2.1
library(ggplot2) # v4.0.3

# read in the sample context tsv using the weighted averages
contexts <- read.table("context_level_counts.tsv", sep = '\t', header = TRUE) %>%
  mutate(group = case_when(startsWith(sample, 'nb') ~ 'bmc_parent', startsWith(sample, 'nc') ~ 'placebo_parent', 
                           startsWith(sample, 'sb') ~ 'bmc_sperm', startsWith(sample, 'sc') ~ 'placebo_sperm')) %>%
  filter(context != "other") %>%
  dplyr::select(sample, group, context, weighted_average)

# Reorder levels
contexts$group <- factor(contexts$group, levels = c("placebo_parent", "bmc_parent", "placebo_sperm", "bmc_sperm"))
contexts$context <- factor(contexts$context, levels = c('exon', 'intron', 'intergenic', 'total'))

# plot the boxplot
p<-ggplot(contexts, aes(x= context, y=weighted_average))+
  geom_boxplot(aes(fill = group), colour= "black", linewidth=0.3, fatten = 0.3) +
  scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40"), labels = c("Placebo Parent", "BMC Parent", "Placebo Sperm", "BMC Sperm"), name = NULL) +
  geom_point(aes(group = group), size=0.6, position = position_dodge(width=0.75)) +
  facet_wrap(vars(context), scales= "free_x", nrow = 1, ncol = 4) +
  theme_classic() +
  scale_y_continuous(breaks = seq(4, 11 , by = 1), limits = c(4,11)) +
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
exons <- contexts %>%
  filter(context == 'exon') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'parent', startsWith(sample, 's') ~ 'sperm')) %>%
  mutate(inoculation = case_when(startsWith(sample, 'nb')|startsWith(sample, 'sb') ~ 'bmc', startsWith(sample, 'nc')|startsWith(sample, 'sc') ~ 'placebo'))
hist(exons$weighted_average)
shapiro.test(exons$weighted_average) # W = 0.92945, p-value = 0.3743
t<-aov(weighted_average ~ generation * inoculation, data = exons)
summary(t)
#                         DfSum Sq  Mean Sq    F value      Pr(>F)  
# generation              1  0.008    0.008      0.007      0.9358  
# inoculation             1  4.700    4.700      4.180      0.0751 .
# generation:inoculation  1  1.089    1.089      0.969      0.3538  
# Residuals               8  8.994    1.124      

introns <- contexts %>%
  filter(context == 'intron') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'parent', startsWith(sample, 's') ~ 'sperm')) %>%
  mutate(inoculation = case_when(startsWith(sample, 'nb')|startsWith(sample, 'sb') ~ 'bmc', startsWith(sample, 'nc')|startsWith(sample, 'sc') ~ 'placebo'))
hist(introns$weighted_average)
shapiro.test(introns$weighted_average) # W = 0.91985, p-value = 0.2847
t<-aov(weighted_average ~ generation * inoculation, data = introns)
summary(t)
#                         Df Sum Sq Mean Sq F value Pr(>F)
# generation              1 0.0384 0.03843   0.471  0.512
# inoculation             1 0.0629 0.06294   0.772  0.405
# generation:inoculation  1 0.0283 0.02833   0.347  0.572
# Residuals               8 0.6523 0.08154

inter <- contexts %>%
  filter(context == 'intergenic') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'parent', startsWith(sample, 's') ~ 'sperm')) %>%
  mutate(inoculation = case_when(startsWith(sample, 'nb')|startsWith(sample, 'sb') ~ 'bmc', startsWith(sample, 'nc')|startsWith(sample, 'sc') ~ 'placebo'))
hist(inter$weighted_average)
shapiro.test(inter$weighted_average) # W = 0.77902, p-value = 0.005463 ***not normal
rcompanion::scheirerRayHare(weighted_average ~ generation * inoculation, data = inter) 
#                         Df      Sum Sq            H         p.value
# generation               1     108.000       8.3077         0.00395 ** generation significant
# inoculation              1       1.333       0.1026         0.74877
# generation:inoculation   1       3.000       0.2308         0.63095
# Residuals                8      30.667
wilcox.test(weighted_average ~ generation, data = inter, conf.int = T) # W = 0, p-value = 0.002165

totals <- contexts %>%
  filter(context == 'total') %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'parent', startsWith(sample, 's') ~ 'sperm')) %>%
  mutate(inoculation = case_when(startsWith(sample, 'nb')|startsWith(sample, 'sb') ~ 'bmc', startsWith(sample, 'nc')|startsWith(sample, 'sc') ~ 'placebo'))
hist(totals$weighted_average)
shapiro.test(totals$weighted_average) # W = 0.84437, p-value = 0.0313 ***not normal
rcompanion::scheirerRayHare(weighted_average ~ generation * inoculation, data = totals) 
#                       Df  Sum Sq      H p.value
# generation              1 108.000 8.3077 0.00395 ** generation significant
# inoculation             1   5.333 0.4103 0.52184
# generation:inoculation  1   3.000 0.2308 0.63095
# Residuals               8  26.667 
wilcox.test(weighted_average ~ generation, data = inter, conf.int = T) # W = 0, p-value = 0.002165
#

## Correlation matrix ----

library(dplyr) # v1.2.1
library(tibble) # v3.3.1
library(data.table) # v1.18.4
library(pheatmap) # v1.0.13
library(ggplot2) # v4.0.3
library(reshape2) #v1.4.5

# read in mean meths tsv fil to run glms
filt.pcts <- fread("all_clair3_pct_context.tsv.gz", sep = '\t', header = TRUE)

# Filter out the data for the correlation plot
df <- filt.pcts %>%
  mutate(scaf_pos = paste(.$scaffold, .$pos, sep = "_")) %>%
  dplyr::select(scaf_pos, nb22:sc5) %>%
  column_to_rownames('scaf_pos')

hist(filt.pcts$merged) # not normal, many zeros
corrk <- pcaPP::cor.fk(df)
#fast estimation of kendall's tau rank correlation coefficient
#corrp <- cor(df, method = "kendall", use = "everything")

# set the self-correlations as NA
diag(corrk) <- NA

write.table(corrk, "kendall_tau_table.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

# make tau table
corrk_long <- melt(corrk, varnames = c("sample1", "sample2"), value.name = "tau")
corrk_long <- corrk_long[!is.na(corrk_long$tau), ] # drops the NA diagonal
corrk_long <- corrk_long[as.character(corrk_long$sample1) < as.character(corrk_long$sample2), ]

write.table(corrk_long, "kendall_correlation_pairs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# read in kendall tau table
taus <- fread("kendall_correlation_pairs.tsv", sep = '\t', header = TRUE)

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
dev.off()
#

## PCA analysis of methylation ----

library(dplyr) # v1.2.1
library(data.table) # v1.18.4
library(factoextra) # v2.2.0
library(tibble) # v3.3.1

# read in percent methylation tsv file
filt.pcts <- fread("all_clair3_pct_context.tsv.gz", sep = '\t', header = TRUE)

# Apply filter functions to just the CpG O/E genes
tots <- filt.pcts %>%
  mutate(scaf_pos = paste(.$scaffold, .$pos, sep = "_")) %>%
  select(scaf_pos, nb22:sc5) %>%
  column_to_rownames('scaf_pos') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  mutate(group = case_when(startsWith(sample, 'nb') | startsWith(sample, 'nc') ~ 'parent',
                           startsWith(sample, 'sb') | startsWith(sample, 'sc') ~ 'sperm')) %>%
  mutate(pairs = case_when(endsWith(sample, 'c3') ~ 'c3', endsWith(sample, 'c4') ~ 'c4', endsWith(sample, 'c5') ~ 'c5', 
                           endsWith(sample, 'b22') ~ 'b22', endsWith(sample, 'b26') ~ 'b26', endsWith(sample, 'b27') ~ 'b27')) %>%
  select(group, pairs, everything()) %>%
  select(-sample)

num <- as.matrix(tots[setdiff(names(tots), c("group", "pairs"))])
num <- qlogis(pmin(pmax(num / 100, 1e-6), 1 - 1e-6))
num <- num[, apply(num, 2, sd) > 0, drop = FALSE]

# Run the PCA
res.pca <- prcomp(num, center = TRUE, scale. = TRUE)

# Get the PCA coordinates
ind_df <- data.frame(PC1 = res.pca$x[, 1], PC2 = res.pca$x[, 2], group = tots$group, pairs = tots$pairs)
pc_var <- (res.pca$sdev^2) / sum(res.pca$sdev^2) * 100

# Add the grouping variables
ind_df$group <- tots$group
ind_df$pairs <- tots$pairs

# Plot with ggplot2
pca_plot <- ggplot(ind_df, aes(x = PC1, y = PC2, color = pairs, shape = group)) +
  geom_point(size = 5, alpha = 1) +
  scale_colour_manual(name = "pairs",
                      values = c(b22 = "#332288", b26 = "#117733", b27 = "#44AA99",
                                 c3 = "#88CCEE",  c4 = "#DDCC77",  c5 = "#882255")) +
  scale_shape_manual(name = "group", values = c(parent = 19, sperm = 1)) +
  labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
       y = paste0("PC2 (", round(pc_var[2], 1), "%)")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text        = element_text(size=14, family="sans", colour="black"),
        axis.title       = element_text(size=14, family="sans", colour="black"),
        legend.title     = element_text(size=14, family="sans", colour="black"),
        legend.text      = element_text(size=12, family="sans", colour="black")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey10") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey10")

ggsave("PCA_scaf_positions_all.pdf", width=8, height=6)
#

## GLM ----

library(dplyr) # v1.1.4
library(tibble) # v3.2.1
library(tidyr) # v1.3.1
library(purrr) # v1.0.4
library(broom) # v1.0.7

# read in mean meths tsv fil to run glms
gene.counts <- read.table("gene_level_counts.tsv", sep = '\t', header = TRUE)
gene.universe <- read.delim("gene_universe_conserved_3CpGs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# reformat the mean methylation data to run the glms
gene.counts.f <- gene.counts %>%
  filter(gene %in% gene.universe) %>%
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'adult', startsWith(sample, 's') ~ 'sperm'),
         inoculation = case_when(startsWith(sample, 'nb') | startsWith(sample, 'sb') ~ 'bmc',
                                 startsWith(sample, 'nc') | startsWith(sample, 'sc') ~ 'placebo')) %>%
  select(sample, generation, inoculation, everything())

gene.counts.f$generation <- factor(gene.counts.f$generation, levels = c('adult', 'sperm'))
gene.counts.f$inoculation <- factor(gene.counts.f$inoculation, levels = c('placebo', 'bmc'))

# run the glm on each gene
glm.results <- gene.counts.f %>%
  nest(data = -gene) %>%
  mutate(
    model = map(data, ~ {
      d <- .
      d$unmeth_count <- d$total_count - d$methylated_count
      
      # Primary model: additive only (interaction term resulted in 0 significant genes)
      fit <- tryCatch(
        glm(cbind(methylated_count, unmeth_count) ~ generation + inoculation,
            data = d, family = binomial()),
        error = function(e) e
      )
      
      if (inherits(fit, "error")) {
        return(tibble(term = NA, estimate = NA, std.error = NA,
                      statistic = NA, p.value = NA, dispersion = NA))
      }
      
      disp <- sum(residuals(fit, type = "pearson")^2) / fit$df.residual
      
      if (!is.na(disp) && disp > 1.5) {
        fit <- glm(cbind(methylated_count, unmeth_count) ~ generation + inoculation,
                   data = d, family = quasibinomial())
      }
      
      tidy(fit) %>% mutate(dispersion = disp)
    })
  ) %>%
  unnest_wider(model) %>%
  select(-data)

# Reshape results
res <- glm.results %>%
  unnest(c(term, estimate, std.error, statistic, p.value, dispersion)) %>%
  filter(!if_all(c(estimate, std.error, statistic, p.value), is.na))

res.wide <- pivot_wider(res, names_from = term,
                        values_from = c(estimate, std.error, statistic, p.value))
res.wide <- arrange(res.wide, gene)

# How many genes were fit with quasibinomial (overdispersed) vs binomial
n_overdispersed <- res %>% distinct(gene, dispersion) %>% filter(dispersion > 1.5) %>% nrow() #4051 of 10763

# Adjust for multiple testing (e.g., Benjamini-Hochberg) for each term
glm.gen <- res.wide %>%
  filter(!is.na(`p.value_generationsperm`)) %>%
  mutate(bh = p.adjust(`p.value_generationsperm`, method = 'BH')) %>%
  filter(bh < 0.05)

write.table(glm.gen$gene, "glm_gen_converged_bh.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(glm.gen, "glm_gen.tsv", sep = "\t", row.names = FALSE)

#--

glm.bmc <- res.wide %>%
  filter(!is.na(`p.value_inoculationbmc`)) %>%
  mutate(bh = p.adjust(`p.value_inoculationbmc`, method = 'BH')) %>%
  filter(bh < 0.05)

write.table(glm.bmc$gene, "glm_bmc_converged_bh.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(glm.bmc, "glm_bmc.tsv", sep = "\t", row.names = FALSE)
#

## CpG of GLMs ridgeplot ----

library(dplyr) # v1.2.1
library(purrr) # v1.2.2
library(ggplot2) # v4.0.3

# read in mean meths tsv fil to run glms
agland<-read.delim("acropora_cpg_bias.tsv", header=TRUE)
gene.universe <- read.delim("gene_universe_conserved_3CpGs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm.gen <- read.delim("glm_gen_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
sperm.hyper <- read.delim("sperm_hypermethylated_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
sperm.hypo <- read.delim("sperm_hypomethylated_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm.bmc <- read.delim("glm_bmc_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# filter the cpg data
cpg.genes <- agland %>%
  filter(gene %in% gene.universe) %>% # 10763 genes
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>% # same filter as before (6 genes removed)
  select(gene, CpG.bias)
cpg.genes$CpG.bias <- as.numeric(cpg.genes$CpG.bias) # median = 0.8161

gen.genes <- agland %>%
  filter(gene %in% glm.gen) %>% # 868 genes
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>% # same filter as before (0 genes removed)
  select(gene, CpG.bias)
gen.genes$CpG.bias <- as.numeric(gen.genes$CpG.bias) # median = 0.62125

hyper.cpg <- agland %>%
  filter(gene %in% sperm.hyper) %>% # 4 genes
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>% # same filter as before (0 genes removed)
  select(gene, CpG.bias)
hyper.cpg$CpG.bias <- as.numeric(hyper.cpg$CpG.bias) # median = 0.5263

hypo.cpg <- agland %>%
  filter(gene %in% sperm.hypo) %>% # 4 genes
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>% # same filter as before (0 genes removed)
  select(gene, CpG.bias)
hypo.cpg$CpG.bias <- as.numeric(hypo.cpg$CpG.bias) # median = 0.8797

bmc.genes <- agland %>%
  filter(gene %in% glm.bmc) %>% # 4 genes
  filter(CpG.bias >= 0.001 & CpG.bias <= 1.5) %>% # same filter as before (0 genes removed)
  select(gene, CpG.bias)
bmc.genes$CpG.bias <- as.numeric(bmc.genes$CpG.bias) # median = 0.70965

# make list of dataframes to overlay
cpg.cats <- list(totals = cpg.genes, glm_gen = gen.genes, glm_bmc = bmc.genes)

# stack them all into a tibble
long.cpg.cats <- imap_dfr(cpg.cats, ~ tibble(CpG.bias = .x$CpG.bias, source = .y))
long.cpg.cats$source <- factor(long.cpg.cats$source, levels = c('glm_bmc', 'glm_gen', 'totals'))

# test variances across the groups
car::leveneTest(CpG.bias ~ source, data = long.cpg.cats) # unequal variances, use nonparametric
#           Df F value    Pr(>F)    
# group     2   16.42 7.563e-08 ***
kruskal.test(CpG.bias ~ source, data = long.cpg.cats) # chi-squared = 326.91, df = 2, p-value < 2.2e-16
FSA::dunnTest(CpG.bias ~ source, data = long.cpg.cats, method = "bh")

# plot
ggplot(long.cpg.cats, aes(x = CpG.bias, y = source, fill = source)) +
  ggridges::geom_density_ridges(scale = 0.8, rel_min_height = 0.001, alpha = 0.5) +
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
  scale_y_discrete(labels = c("GLM - inoculation", "GLM - life stage", "Total genes"),
                   expand = expansion(add = c(0.1)),
                   position = 'right') +
  labs(x = expression("CpG"[O/E]), y = NULL)

ggsave("cpg_genes_ridgeplot.pdf", width=8, height=6, dpi = 300)
#

## Differential abundances of genes ----

library(dplyr) # v1.2.1
library(tibble) # v3.3.1
library(tidyr) # v1.3.2
library(broom) # v1.0.13
library(pheatmap) # v1.0.13
library(ggplot2) # v4.0.3

# read in mean meths tsv fil to run glms
weighted.avgs <- read.table("gene_level_counts.tsv", sep = '\t', header = TRUE)
gene.universe <- read.delim("gene_universe_conserved_3CpGs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm.gen <- read.delim("glm_gen_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()
glm.bmc <- read.delim("glm_bmc_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# reformat the mean methylation data to run the glms
weighted.avgs.f <- weighted.avgs %>%
  filter(gene %in% gene.universe) %>% # 10763
  mutate(generation = case_when(startsWith(sample, 'n') ~ 'adult', startsWith(sample, 's') ~ 'sperm'),
         inoculation = case_when(startsWith(sample, 'nb') | startsWith(sample, 'sb') ~ 'bmc',
                                 startsWith(sample, 'nc') | startsWith(sample, 'sc') ~ 'placebo')) %>%
  select(sample, generation, inoculation, gene, weighted_average)

weighted.avgs.f$generation <- factor(weighted.avgs.f$generation, levels = c('adult', 'sperm'))
weighted.avgs.f$inoculation <- factor(weighted.avgs.f$inoculation, levels = c('placebo', 'bmc'))

# run the heatmap for the glm_gen and glm_bmc
gen.meths <- weighted.avgs.f %>% 
  filter(gene %in% glm.gen) %>%
  pivot_wider(id_cols = c('sample', 'generation', 'inoculation'), names_from = 'gene', values_from = 'weighted_average') %>%
  column_to_rownames('sample')

# plot the heatmap
heat.gen <- t(gen.meths[,3:ncol(gen.meths)]) # just the data, adjust the group as needed
heat.gen <- t(scale(t(heat.gen))) # scale before pheatmap because pheatmap scaling doesn't separate color enough
heat.gen.gs <- select(gen.meths, generation, inoculation) # groups

heat.gen.gs$generation <- factor(heat.gen.gs$generation, levels = c('adult', 'sperm'))
heat.gen.gs$inoculation <- factor(heat.gen.gs$inoculation, levels = c('placebo', 'bmc'))

pheatmap(heatr,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, 'RdYlBu'))(256)),
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2', 
         cutree_rows = 1, cutree_cols = 2,
         show_rownames = FALSE,
         angle_col = 45,
         border_color= NA,
         annotation_col = heatr.gs,
         filename = "heatmap_gen_glm.pdf",
         width = 6,
         height = 6
         )

#--

# run the heatmap for the glm_bmc
bmc.meths <- weighted.avgs.f %>% 
  filter(gene %in% glm.bmc) %>%
  pivot_wider(id_cols = c('sample', 'generation', 'inoculation'), names_from = 'gene', values_from = 'weighted_average') %>%
  column_to_rownames('sample')

# plot the heatmap
heat.bmc <- t(bmc.meths[,3:ncol(bmc.meths)]) # just the data, adjust the group as needed
heat.bmc <- t(scale(t(heat.bmc))) # scale before pheatmap because pheatmap scaling doesn't separate color enough
heat.bmc.gs <- select(bmc.meths, generation, inoculation) # groups

heat.bmc.gs$generation <- factor(heat.bmc.gs$generation, levels = c('adult', 'sperm'))
heat.bmc.gs$inoculation <- factor(heat.bmc.gs$inoculation, levels = c('placebo', 'bmc'))

pheatmap(heatr,
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(10, 'RdYlBu'))(256)),
         clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2', 
         cutree_rows = 1, cutree_cols = 2,
         show_rownames = FALSE,
         angle_col = 45,
         border_color= NA,
         annotation_col = heatr.gs,
         filename = "heatmap_bmc_glm.pdf",
         width = 6,
         height = 3
         )
#

## topGO of gen GLMs ----

library(topGO) # v2.62.0

# remember to change the folder name to point to the folder containing your genes of interest lists
folder_of_interest = "heatmap_genes/"

# exclude files with "universe" in it
mult_files = grep(list.files(folder_of_interest), pattern="*genes.txt", value = T)
universe_file = grep(list.files(folder_of_interest), pattern="*universe.txt", value =T)

for (m in mult_files) {
  annot_filename = 'acropora_ref_goIDs.tsv'
  gene_id_to_go = readMappings(file=annot_filename)
  
  # shrink list of all GO terms down to the correct universe
  #universe_file = gsub('up', 'universe', m)
  #universe_file = gsub('down', 'universe', universe_file)
  #universe_file = gsub('diff', 'universe', universe_file)
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
    output_filename = paste0("heatmap_genes/", go_category, "_", m)
    write.table(results_table, file=output_filename, quote=FALSE, sep='\t')
  }
}
#

## DSS differential methylation analysis ----

library(DSS) # v2.58.0
library(dplyr) # v1.2.1
require(bsseq) # v1.46.0

# read in the design file
design <- read.table("design.txt", sep = '\t', header = TRUE)
design
design$generation <- factor(design$generation, levels = c("adult", "sperm"))
levels(design$generation)
design$inoculation <- factor(design$inoculation, levels = c("placebo", "bmc"))
levels(design$inoculation)

# read in all the methylation data from the files
dss.dir <- "dss_input"

meth.list <- lapply(design$sample, function(s) {
  read.table(file.path(dss.dir, paste0(s, "_dss_input.txt")),
             sep = "\t", header = TRUE)
})
names(meth.list) <- design$sample

# make the bs object with the files
BSobj <- makeBSseqData(meth.list, sampleNames = design$sample)

# observe the model matrix
X = model.matrix(~generation+inoculation, design)

# fit the model (similar to multiple regression); model only needs to be fit once
DMLfit = DMLfit.multiFactor(BSobj, design=design, 
                            formula=~generation+inoculation)

# check to see which column refers to which comparison
colnames(DMLfit$X)

# look to see the positions affected by each factor 
DMLtest.gen <- DMLtest.multiFactor(DMLfit, coef = 2)
DMLtest.inoc <- DMLtest.multiFactor(DMLfit, coef = 3)

# identify DMCs
dmc.gen <- DMLtest.gen %>%
  as.data.frame() %>%
  filter(fdrs < 0.05) %>%
  arrange(fdrs)
# 61667 DMCs
dmc.inoc <- DMLtest.inoc %>%
  as.data.frame() %>%
  filter(fdrs < 0.05) %>%
  arrange(fdrs)
# 4276 DMCs

write.table(dmc.gen, "dss_dmc_generation.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dmc.inoc, "dss_dmc_inoculation.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# now we can call different methylated regions
gen.dmrs = callDMR(DMLtest.gen, p.threshold=0.001) # 3362 DMRs
bmc.dmrs = callDMR(DMLtest.inoc, p.threshold=0.001) # 400 DMRs

# visualize DMR
showOneDMR(bmc.dmrs[1,], BSobj)

# save to files
write.table(gen.dmrs, "dss_dmrs_generation.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(bmc.dmrs, "dss_dmrs_inoculation.txt", sep="\t", row.names=FALSE, quote = FALSE)

# create a Granges object from results to get the sequences associated with the DMRs
library(GenomicRanges) # v1.62.1

dmr.gr <- GRanges(seqnames = bmc.dmrs$chr,
                  ranges = IRanges(start = bmc.dmrs$start, end = bmc.dmrs$end), strand = "*")

mcols(dmr.gr)$length <- bmc.dmrs$length
mcols(dmr.gr)$nCG <- bmc.dmrs$nCG
mcols(dmr.gr)$areaStat <- bmc.dmrs$areaStat
dmr.gr

library(Rsamtools) # v2.26.0
library(Biostrings) # v2.78.0

genome <- FaFile("acropora_ref_assembly.fa")
open(genome)

dmr.seqs <- getSeq(genome, dmr.gr)
dmr.seqs

writeXStringSet(dmr.seqs, "inoculation_dmr_sequences.fa")
#

## Find overlapping DMRs between two tests  ----

library(GenomicRanges) # v1.62.1
library(dplyr) # v1.2.1
library(rtracklayer) #1.70.1

# dmr table
bmc.dmrs <- read.table("dss_dmrs_inoculation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gen.dmrs <- read.table("dss_dmrs_generation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

dmr.gen <- GRanges(seqnames = gen.dmrs$chr,
                  ranges = IRanges(start = gen.dmrs$start, end = gen.dmrs$end), strand = "*")

mcols(dmr.gen)$length <- gen.dmrs$length
mcols(dmr.gen)$nCG <- gen.dmrs$nCG
mcols(dmr.gen)$areaStat <- gen.dmrs$areaStat

dmr.bmc <- GRanges(seqnames = bmc.dmrs$chr,
                      ranges = IRanges(start = bmc.dmrs$start, end = bmc.dmrs$end),
                      strand = "*")
mcols(dmr.bmc)$length   <- bmc.dmrs$length
mcols(dmr.bmc)$nCG      <- bmc.dmrs$nCG
mcols(dmr.bmc)$areaStat <- bmc.dmrs$areaStat

overlaps.bmc.gen <- findOverlaps(dmr.bmc, dmr.gen)

n.bmc.overlap <- n_distinct(queryHits(overlaps.bmc.gen)) # 58 overlaps
n.gen.overlap  <- n_distinct(subjectHits(overlaps.bmc.gen)) # 57 overlaps

dmr.overlap <- data.frame(
  inoc_chr       = as.character(seqnames(dmr.bmc))[queryHits(overlaps.bmc.gen)],
  inoc_start     = start(dmr.bmc)[queryHits(overlaps.bmc.gen)],
  inoc_end       = end(dmr.bmc)[queryHits(overlaps.bmc.gen)],
  bmc_areaStat  = mcols(dmr.bmc)$areaStat[queryHits(overlaps.bmc.gen)],
  gen_chr        = as.character(seqnames(dmr.gen))[subjectHits(overlaps.bmc.gen)],
  gen_start      = start(dmr.gen)[subjectHits(overlaps.bmc.gen)],
  gen_end        = end(dmr.gen)[subjectHits(overlaps.bmc.gen)],
  gen_areaStat   = mcols(dmr.gen)$areaStat[subjectHits(overlaps.bmc.gen)]
)

write.table(dmr.overlap, "dmr_inoculation_generation_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#

## Determine whether DMRs fall within gene bodies/promoters ----

library(GenomicRanges) # v1.62.1
library(dplyr) # v1.2.1
library(rtracklayer) #1.70.1

# dmr table
bmc.dmrs <- read.table("dss_dmrs_inoculation.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
glm.bmc <- read.delim("glm_bmc_converged_bh.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1] %>% trimws()

# gff file
gff <- import("acropora_ref_assembly.gff")

dmr.gr <- GRanges(seqnames = bmc.dmrs$chr,
                  ranges = IRanges(start = bmc.dmrs$start, end = bmc.dmrs$end), strand = "*")

mcols(dmr.gr)$length <- bmc.dmrs$length
mcols(dmr.gr)$nCG <- bmc.dmrs$nCG
mcols(dmr.gr)$areaStat <- bmc.dmrs$areaStat
dmr.gr

genes.gr <- gff[gff$type == "gene"]
mcols(genes.gr)$gene_id <- unlist(mcols(genes.gr)$ID)
overlaps <- findOverlaps(dmr.gr, genes.gr)

dmr.gene.overlap <- data.frame(
  dmr_chr    = as.character(seqnames(dmr.gr))[queryHits(overlaps)],
  dmr_start  = start(dmr.gr)[queryHits(overlaps)],
  dmr_end    = end(dmr.gr)[queryHits(overlaps)],
  dmr_areaStat = mcols(dmr.gr)$areaStat[queryHits(overlaps)],
  dmr_nCG    = mcols(dmr.gr)$nCG[queryHits(overlaps)],
  gene_id    = mcols(genes.gr)$gene_id[subjectHits(overlaps)]
)

dmr.gene.overlap %>%
  distinct(dmr_chr, dmr_start, dmr_end, .keep_all = TRUE) %>%
  dplyr::count() # Of 400 inoculation DMRs, 66 overlap at least one gene body
dmr.gene.overlap %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::count() # 48 unique genes

# write DMRs overlapping in gene bodies to tsv file
write.table(dmr.gene.overlap, "dmr_inoculation_gene_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Confirm strandedness in genes.gr
table(strand(genes.gr))

## Determine whether the genes that are in gene bodies are in exons/introns

mrna.gr <- gff[gff$type == "mRNA"]
mrna_ids <- sapply(mcols(mrna.gr)$ID, function(x) x[1])
mrna_parents <- sapply(mcols(mrna.gr)$Parent, function(x) x[1])
mrna_to_gene <- setNames(mrna_parents, mrna_ids)

# restricting to genes with DMRs
genes_of_interest <- unique(dmr.gene.overlap$gene_id)   # 48 genes

# exons
exons.gr <- gff[gff$type == "exon"]
exon_parent_mrna <- unlist(mcols(exons.gr)$Parent)
exons.gr$gene_id <- mrna_to_gene[exon_parent_mrna]
exons.gr <- exons.gr[exons.gr$gene_id %in% genes_of_interest]

exons_by_gene <- reduce(split(exons.gr, exons.gr$gene_id))

# using the genes with DMRs
valid_genes <- intersect(genes_of_interest, names(exons_by_gene))
exons_by_gene <- exons_by_gene[valid_genes]
genes_matched <- genes.gr[match(valid_genes, mcols(genes.gr)$gene_id)]

# introns by non-exons
introns_by_gene <- psetdiff(genes_matched, exons_by_gene)
introns.gr <- unlist(introns_by_gene)
mcols(introns.gr)$gene_id <- names(introns.gr)

# classifying above DMRs are exon or intron
dmr.gene.overlap.classified <- dmr.gene.overlap %>%
  rowwise() %>%
  mutate(
    in_exon = {
      dmr_r <- GRanges(dmr_chr, IRanges(dmr_start, dmr_end))
      if (!gene_id %in% names(exons_by_gene)) {
        NA
      } else {
        any(overlapsAny(dmr_r, exons_by_gene[[gene_id]]))
      }
    },
    in_intron = {
      dmr_r <- GRanges(dmr_chr, IRanges(dmr_start, dmr_end))
      if (!gene_id %in% names(introns_by_gene)) {
        NA
      } else {
        any(overlapsAny(dmr_r, introns_by_gene[[gene_id]]))
      }
    }
  ) %>%
  ungroup()

table(dmr.gene.overlap.classified$in_exon, dmr.gene.overlap.classified$in_intron,
      useNA = "ifany", dnn = c("in_exon", "in_intron"))
# in exons 30/66 DMRs
# in introns 22/66 DMRs
# in both 16/66 DMRs

write.table(dmr.gene.overlap.classified, "dmr_gene_overlap_exons_introns.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

## Find promoters upstream of genes using the canonical window size
up <- 500
down <- 0

promoter.gr <- promoters(genes.gr, upstream = up, downstream = down)
mcols(promoter.gr)$gene_id <- mcols(genes.gr)$gene_id

overlaps <- findOverlaps(dmr.gr, promoter.gr)

dmr.promoter.overlap <- data.frame(
  dmr_chr      = as.character(seqnames(dmr.gr))[queryHits(overlaps)],
  dmr_start    = start(dmr.gr)[queryHits(overlaps)],
  dmr_end      = end(dmr.gr)[queryHits(overlaps)],
  dmr_areaStat = mcols(dmr.gr)$areaStat[queryHits(overlaps)],
  gene_id      = mcols(promoter.gr)$gene_id[subjectHits(overlaps)]
)

n_dmrs_in_promoters <- n_distinct(dmr.promoter.overlap[, c("dmr_chr","dmr_start","dmr_end")]) # 22 DMRs in promoters
n_genes_touched <- n_distinct(dmr.promoter.overlap$gene_id)
n_glm_overlap <- sum(dmr.promoter.overlap$gene_id %in% glm.bmc)

write.table(dmr.promoter.overlap, "dmr_promoter_overlap.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#

## NanoMethViz to visualize methylation levels across a ROI ----
# adapted from https://www.bioconductor.org/packages/release/bioc/vignettes/NanoMethViz/inst/doc/UsersGuide.html
                       
library(NanoMethViz) # v3.9.0
library(dplyr) # v1.2.1
library(data.table) # v1.18.4
library(tidyr) # v1.3.1
library(Rsamtools) # v2.25.0
library(patchwork) # v1.3.0

# Import the gff data into a GRanges object that can be modified to fit the format of NanoMethViz
anno <- rtracklayer::import("acropora_ref_assembly.gff")
head(anno)

genes <- anno %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  dplyr::rename(chr = "seqnames", gene_id = "ID") %>%
  dplyr::select("gene_id", "chr", "strand", "start", "end") %>%
  filter(!grepl("^tRNA", gene_id)) %>%
  mutate("transcript_id" = gene_id, "symbol" = gene_id)

head(genes)

# Define the samples annotations
sample_anno <- read.table("sample_anno.tsv", sep = '\t', header = TRUE)
                       
# Build the NanoMethResult
nmr <- NanoMethResult(
  methy = "combined_data_clair3.sorted.tsv.gz",
  samples = sample_anno,
  exons = genes
)

# loading saved results from previous bsseq analysis
bsseq_dmr <- read.table("dss_dmrs_inoculation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## DMR regions
a <- plot_region(nmr, "scaffold00015", 1782956, 1814493, smoothing_window = 500, anno_regions = bsseq_dmr)
a[[1]] <- a[[1]] + scale_y_continuous(limits = c(0.499, 0.751))
# g1713; IPR011989 Armadillo-like helical;IPR016024 Armadillo-type fold

b <- plot_region(nmr, "scaffold00018", 1063899, 1069769, smoothing_window = 500, anno_regions = bsseq_dmr)
b[[1]] <- b[[1]] + scale_y_continuous(limits = c(0.499, 0.751))
# this includes g2112+g2113 
# IPR017877 None;IPR026753 None;IPR028002 Myb/SANT-like DNA-binding domain 5
# IPR026103 Harbinger transposase-derived nuclease, animal;IPR027806 Harbinger transposase-derived nuclease domain
# I used the genes, then did 500 bases upstream for the "promoter region"

c <- plot_region(nmr, "scaffold00074", 744480, 747533, smoothing_window = 500, anno_regions = bsseq_dmr)
c[[1]] <- c[[1]] + scale_y_continuous(limits = c(0.499, 0.751))
# includes g7678; hypothetical protein/ no InterPro data

d <- plot_region(nmr, "scaffold00301", 259677, 271354, smoothing_window = 500, anno_regions = bsseq_dmr)
d[[1]] <- d[[1]] + scale_y_continuous(limits = c(0.499, 0.751))
# includes g19100; IPR008521 Magnesium transporter NIPA
dmr_plot <- (a[[1]] / a[[3]] / b[[1]] / b[[3]] / c[[1]] / c[[3]] / d[[1]] / d[[3]]) + # keep top line plot + gene bar, drop the heatmap panel
  plot_layout(heights = rep(c(4, 1), 4))
print(dmr_plot)

ggsave("dmr_plot_nanomethviz.pdf", plot = dmr_plot, width = 8, height = 11)
#
