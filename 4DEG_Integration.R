####Week4
#DESeq2
library(DESeq2)
library(tibble)

#read in DESeq2 results
deseq_results <- read.table("results/GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt.gz", header=TRUE)

#remove rows with NA values
deseq_results <- na.omit(deseq_results)
#filter data
filtered_results <- deseq_results[deseq_results$padj<0.01 & abs(deseq_results$log2FoldChange)>1, ]
#count the number of DE genes
num_de_genes <- nrow(filtered_results)
#print out number
print(num_de_genes)
#1153

#Do they match the numbers reported in the paper? __ IDK cannot find number in paper
 
#Figure 2f
library(ggplot2)
library(dplyr)
library(tidyr)

#read in annotated_peaks table
annotated_peaks <- read.table("results/annotated_peaks.txt", header=TRUE, fill=TRUE, sep="\t", quote="")

#separate up/down regulated
upregulated_genes <- filtered_results[filtered_results$log2FoldChange > 0, ]
upregulated_total <- nrow(upregulated_genes)
downregulated_genes <- filtered_results[filtered_results$log2FoldChange < 0, ]
downregulated_total <- nrow(downregulated_genes)
####counting the number total DE genes not the number of genes that are also annotated

#merge tables
upregulated_genes_TSS <- upregulated_genes %>% 
  inner_join(annotated_peaks, by = c("genename"="Gene.Name"))
downregulated_genes_TSS <- downregulated_genes %>% 
  inner_join(annotated_peaks, by = c("genename"="Gene.Name"))

#5kb and 20kb
upregulated_TSS_5kb <- subset(upregulated_genes_TSS, abs(Distance.to.TSS) <= 5000)
upregulated_TSS_5kb_count <- nrow(upregulated_TSS_5kb)
upregulated_TSS_20kb <- subset(upregulated_genes_TSS, abs(Distance.to.TSS) <= 20000)
upregulated_TSS_20kb_count <- nrow(upregulated_TSS_20kb)
downregulated_TSS_5kb <- subset(downregulated_genes_TSS, abs(Distance.to.TSS) <= 5000)
downregulated_TSS_5kb_count <- nrow(downregulated_TSS_5kb)
downregulated_TSS_20kb <- subset(downregulated_genes_TSS, abs(Distance.to.TSS) <= 20000)
downregulated_TSS_20kb_count <- nrow(downregulated_TSS_20kb)

#create dataframe of counts
figure2f <- data.frame(
  Category = c("5kb upregulated", "5kb downregulated", "20kb upregulated", "20kb downregulated"),
  bound = c(upregulated_TSS_5kb_count, downregulated_TSS_5kb_count, upregulated_TSS_20kb_count, downregulated_TSS_20kb_count),
  not_bound = c(upregulated_total-upregulated_TSS_5kb_count, downregulated_total-downregulated_TSS_5kb_count, upregulated_total-upregulated_TSS_20kb_count, downregulated_total-downregulated_TSS_20kb_count),
  total = c(upregulated_total, downregulated_total, upregulated_total, downregulated_total)
)

#reshape the long
figure2f_long <- figure2f %>%
  gather(key = "Condition", value = "Count", -Category, -total)

#calculate the percentage
figure2f_long <- figure2f_long %>%
  mutate(Percentage = (Count / total) * 100)

#reorder the levels of Category based on the desired order
desired_order <- c("5kb upregulated", "5kb downregulated", "20kb upregulated", "20kb downregulated")
figure2f_long$Category <- factor(figure2f_long$Category, levels = desired_order)

#define colors and their names
colors <- c("bound" = "red", "not_bound" = "grey")

#create the stacked bar plot
ggplot(figure2f_long, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = Count), position = position_stack(reverse = TRUE, vjust = 0.5), size = 3, color = "black") +
  labs(title = "figure2f", y = "Percentage") +
  scale_fill_manual(values = colors) +
  theme_minimal()

#Supplementary figureS2D
#100kb
upregulated_TSS_100kb <- subset(upregulated_genes_TSS, abs(Distance.to.TSS) <= 100000)
upregulated_TSS_100kb_count <- nrow(upregulated_TSS_100kb)
downregulated_TSS_100kb <- subset(downregulated_genes_TSS, abs(Distance.to.TSS) <= 100000)
downregulated_TSS_100kb_count <- nrow(downregulated_TSS_100kb)

#create dataframe of counts
figureS2D <- data.frame(
  Category = c("100kb upregulated", "100kb downregulated"),
  bound = c(upregulated_TSS_100kb_count, downregulated_TSS_100kb_count),
  not_bound = c(upregulated_total-upregulated_TSS_100kb_count, downregulated_total-downregulated_TSS_100kb_count),
  total = c(upregulated_total, downregulated_total)
)

#reshape the long
figureS2D_long <- figureS2D %>%
  gather(key = "Condition", value = "Count", -Category, -total)

#calculate the percentage
figureS2D_long <- figureS2D_long %>%
  mutate(Percentage = (Count / total) * 100)

#define colors and their names
colors <- c("bound" = "red", "not_bound" = "grey")

#create the stacked bar plot
ggplot(figureS2D_long, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = Count), position = position_stack(reverse = TRUE, vjust = 0.5), size = 3, color = "black") +
  labs(title = "figureS2D", y = "Percentage") +
  scale_fill_manual(values = colors) +
  theme_minimal()

##Combined figures
#create dataframe of counts
figureRecreate <- data.frame(
  Category = c("5kb upregulated", "5kb downregulated", "20kb upregulated", "20kb downregulated", "100kb upregulated", "100kb downregulated"),
  bound = c(upregulated_TSS_5kb_count, downregulated_TSS_5kb_count, upregulated_TSS_20kb_count, downregulated_TSS_20kb_count, upregulated_TSS_100kb_count, downregulated_TSS_100kb_count),
  not_bound = c(upregulated_total-upregulated_TSS_5kb_count, downregulated_total-downregulated_TSS_5kb_count, upregulated_total-upregulated_TSS_20kb_count, downregulated_total-downregulated_TSS_20kb_count, upregulated_total-upregulated_TSS_100kb_count, downregulated_total-downregulated_TSS_100kb_count),
  total = c(upregulated_total, downregulated_total, upregulated_total, downregulated_total, upregulated_total, downregulated_total)
)

#reshape the long
figureRecreate_long <- figureRecreate %>%
  gather(key = "Condition", value = "Count", -Category, -total)

#calculate the percentage
figureRecreate_long <- figureRecreate_long %>%
  mutate(Percentage = (Count / total) * 100)

#reorder the levels of Category based on the desired order
desired_order <- c("5kb upregulated", "5kb downregulated", "20kb upregulated", "20kb downregulated", "100kb upregulated", "100kb downregulated")
figureRecreate_long$Category <- factor(figureRecreate_long$Category, levels = desired_order)

#define colors and their names
colors <- c("bound" = "red", "not_bound" = "grey")

#create the stacked bar plot
ggplot(figureRecreate_long, aes(x = Category, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = Count), position = position_stack(reverse = TRUE, vjust = 0.5), size = 3, color = "black") +
  labs(title = "figureRecreate", y = "Percentage") +
  scale_fill_manual(values = colors) +
  theme_minimal()
