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
# 2. خلاصه کلی جهش‌ها
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

# 5. هم‌انحصاری و هم‌رخدادی (آزمون Fisher)
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

ggsave("figures/05_mutation_frequency.png", p, width = 6, height = 4, dpi = 300)
