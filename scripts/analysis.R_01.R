BiocGenerics::install (c("maftools" , "TCGAbiolinks"))
install.packages("dplyr")
library(BiocManager)
BiocManager::install(c("maftools" , "TCGAbiolinks"))
library(maftools)
library(ggplot2)
library(dplyr)
library(TCGAbiolinks)

genes <- c("KRAS" , "TP53" , "CDKN2A" , "SMAD4")
library(maftools)
library(ggplot2)
library(dplyr)

genes <- c("KRAS", "TP53", "CDKN2A", "SMAD4")
dir.create("figures", showWarnings = FALSE)

# 1. خواندن داده
maf <- read.maf(maf = "data_mutations.txt")
while(!is.null(dev.list()))dev.off()
list.files("figures")
View("01_maf_summary.png")
# 2. خلاصه کلی جهش ها
png("figures/01_maf_summary.png", width = 2400, height = 1800, res = 300)
plotmafSummary(maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
dev.off()

# 3. OncoPrint چهار ژن
png("figures/02_oncoplot.png", width = 2400, height = 1200, res = 300)
oncoplot(maf, genes = genes)
dev.off()

# 4. Lollipop plot هر ژن
for (g in genes) {
  png(paste0("figures/03_lollipop_", g, ".png"),
      width = 2400, height = 1200, res = 300)
  lollipopPlot(maf, gene = g, AACol = "HGVSp_Short", showMutationRate = TRUE)
  dev.off()
}

# 5. همانحصاری و همرخدادی (آزمون Fisher)
png("figures/04_interactions.png", width = 1800, height = 1800, res = 300)
somaticInteractions(maf, genes = genes, pvalue = c(0.05, 0.1))
dev.off()

# 6. نمودار فراوانی جهش با ggplot2
n_samples <- nrow(getSampleSummary(maf))
freq <- as.data.frame(getGeneSummary(maf)) %>%
  filter(Hugo_Symbol %in% genes) %>%
  mutate(percent = MutatedSamples / n_samples * 100)

p <- ggplot(freq, aes(x = reorder(Hugo_Symbol, -percent), y = percent)) +
  geom_col(fill = "#2c6e91") +
  geom_text(aes(label = paste0(round(percent, 1), "%")), vjust = -0.4) +
  labs(title = "Mutation frequency in TCGA-PAAD",
       x = "Gene", y = "Mutated samples (%)") +
  theme_minimal(base_size = 14)

setwd("..")
getwd()
while (!is.null(dev.list())) dev.off()
somaticInteractions(maf, genes = genes, pvalue = c(0.05, 0.1))
ggsave("figures/05_mutation_frequency.png", p, width = 6, height = 4, dpi = 300)

# جدول جهش: هر نمونه، هر ژن (TRUE = جهش دارد)
mut <- as.data.frame(maf@data) %>%
  filter(Hugo_Symbol %in% genes) %>%
  distinct(Tumor_Sample_Barcode, Hugo_Symbol)

all_samples <- as.character(getSampleSummary(maf)$Tumor_Sample_Barcode)
mat <- sapply(genes, function(g)
  all_samples %in% as.character(mut$Tumor_Sample_Barcode[mut$Hugo_Symbol == g]))

# آزمون Fisher برای هر جفت ژن
pairs <- combn(genes, 2, simplify = FALSE)
results <- do.call(rbind, lapply(pairs, function(p) {
  tab <- table(factor(mat[, p[1]], levels = c(FALSE, TRUE)),
               factor(mat[, p[2]], levels = c(FALSE, TRUE)))
  ft <- fisher.test(tab)
  data.frame(gene1 = p[1], gene2 = p[2], both_mutated = tab[2, 2],
             odds_ratio = round(unname(ft$estimate), 2),
             p_value = signif(ft$p.value, 3))
}))
results$p_adj <- round(p.adjust(results$p_value, method = "BH"), 3)

results
write.csv(results, "figures/04_interactions_table.csv", row.names = FALSE)
