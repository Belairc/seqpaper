setwd("~/Box Sync/SeqPaperDataCB/RAnalysis")
# analysis of Mature counts of control cell 4 in three libraries TJ12, TJ12hiPCR and TJ11
library('tidyverse')
theme_set(theme_minimal(base_size=8))
# sample <- "4N01_ccell4_500ng_1_Mature"

Nodedup12 <- read.table("/Volumes/BiomarkerSeq/TJ12/TJ12_COUNTS_raw/4N01_ccell4_500ng_1_Mature_counts.txt", stringsAsFactors = FALSE)
Nodedup12[is.na(Nodedup12)] <- 0
colnames(Nodedup12) <-  c("miRNA","ccell4_Mature_TJ12")
# Nodedup12 <- Nodedup12 %>% mutate(ccell4_TJ12_rpm = ccell4_Mature_TJ12/(sum(ccell4_Mature_TJ12)/1000000)) %>% select(-ccell4_Mature_TJ12)

dedup12 <-  read.table("/Volumes/BiomarkerSeq/TJ12/TJ12_collapse12/4N01_ccell4_500ng_1_Mature_count12.txt", stringsAsFactors = FALSE)
dedup12[is.na(dedup12)] <- 0
colnames(dedup12) <- c("miRNA","ccell4_Mature_TJ12_col12")

Nodedup12hiPCR <-  read.table("/Volumes/BiomarkerSeq/TJ12hi/TJ12hi_COUNTS_raw/4N01_ccell4_500ng_1_Mature_counts.txt",stringsAsFactors = FALSE)
Nodedup12hiPCR[is.na(Nodedup12hiPCR)] <- 0
colnames(Nodedup12hiPCR) <- c("miRNA","ccell4_Mature_TJ12hiPCR")

dedup12hiPCR <-  read.table("/Volumes/BiomarkerSeq/TJ12hi/TJ12hi_collapse12/4N01_ccell4_500ng_1_Mature_count12.txt", stringsAsFactors = FALSE)
dedup12hiPCR[is.na(dedup12hiPCR)] <- 0
colnames(dedup12hiPCR) <- c("miRNA","ccell4_Mature_TJ12hiPCR_col12")

NodedupTJ11 <-  read.table("/Volumes/BiomarkerSeq/TJ11/TJ11_COUNTS_raw/4N37_ctrl_cell_4_Mature_counts.txt", stringsAsFactors =FALSE)
NodedupTJ11[is.na(NodedupTJ11)] <- 0
colnames(NodedupTJ11) <- c("miRNA","ccell4_Mature_TJ11")

dedupTJ11 <-  read.table("/Volumes/BiomarkerSeq/TJ11/TJ11_collapse12/4N37_ctrl_cell_4_Mature_count12.txt", stringsAsFactors = FALSE)
dedupTJ11[is.na(dedupTJ11)] <- 0
colnames(dedupTJ11) <- c("miRNA","ccell4_Mature_TJ11_col12")



setwd("~/Box Sync/SeqPaperDataCB/RAnalysis/PCRcycle_cell4")

# ---- inner join --------
compare_counts <- list(Nodedup12, dedup12, Nodedup12hiPCR, dedup12hiPCR, NodedupTJ11, dedupTJ11 ) %>% Reduce( function(dtf1,dtf2) inner_join(dtf1,dtf2, by = "miRNA"), .)
compare_counts[is.na(compare_counts)] <-0
# mutate(ccell4_TJ12_rpm = ccell4_Mature_TJ12/(sum(ccell4_Mature_TJ12)/1000000))
cpm <- function(x)
{
  x/sum(x)*10000
}

compare_cpm <- compare_counts %>% select(-miRNA) %>% mutate_all(funs(cpm,"cpm",cpm=(.))) %>% select(ends_with("_cpm")) %>% mutate(miRNA=compare_counts$miRNA)
cutoff <- 10

# log10 axis ---------
TJ12<- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12_col12_cpm + .5), y = log10(ccell4_Mature_TJ12_cpm + .5))) +
  theme_bw(base_size = 8) + geom_point(size=.1) +
  xlab('Unique log10(cpm + .5)') + ylab('Raw log10(cpm + .5)') + ggtitle('Raw vs Collapsed cpm 14 Cycles') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))

TJ12hi <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12hiPCR_col12_cpm + .5), y = log10(ccell4_Mature_TJ12hiPCR_cpm + .5))) + 
  theme_bw(base_size = 8) + geom_point(size = .1) +
  xlab('Unique log10(cpm + .5)') + ylab('Raw log10(cpm + .5)') + ggtitle('Raw vs Collapsed cpm 24 Cycles') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))

 TJ11 <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ11_col12_cpm + .5), y = log10(ccell4_Mature_TJ11_cpm + .5))) + geom_point(size=.1) +
  xlab('TJ11 Dedup log10(cpm + .5)') + ylab('TJ11 NO dedup log10(cpm + .5)') + ggtitle('Control cell 4 TJ11\n Collapse on Complete UMI') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))

 TJ12_TJ12hi <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12_cpm + .5), y = log10(ccell4_Mature_TJ12hiPCR_cpm + .5))) + geom_point(size=.1) +
  xlab('TJ12 log10(cpm + .5)') + ylab('TJ12hiPCR log10(cpm + .5)') + ggtitle('Control cell 4\n TJ12 vs TJ12hiPCR ') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))

 TJ12_TJ12hi_coll <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12_col12_cpm + .5), y = log10(ccell4_Mature_TJ12hiPCR_col12_cpm + .5))) + geom_point(size=.05) +
  xlab('14 Cycles Unique log10(cpm + .5)') + ylab('24 Cycles Unique log10(cpm + .5)') + ggtitle('14 vs 24 Cycles Collapsed') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))

 TJ12_TJ11_coll <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12_col12_cpm + .5), y = log10(ccell4_Mature_TJ11_col12_cpm + .5))) + geom_point(size=.1) +
  xlab('TJ12_collapsed log10(cpm + .5)') + ylab("TJ11_collapsed log10(cpm + .5)") + ggtitle('Control cell 4\n TJ12  vs TJ11 collapsed') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))

 TJ12hi_TJ11_coll <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ11_col12_cpm + .5), y = log10(ccell4_Mature_TJ12hiPCR_col12_cpm + .5))) + geom_point(size=.1) +
  xlab('TJ11_collapsed log10(cpm + .5)') + ylab('TJ12hiPCR_collapsed log10(cpm + .5)') + ggtitle('Control cell 4\n TJ11  vs TJ12hiPCR collapsed') + 
  geom_vline(xintercept = log10(cutoff + .5)) + geom_hline(yintercept = log10(cutoff + .5))


## try to plot 6 on one page
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#
print(TJ12, vp = vplayout(1, 1))  # key is to define vplayout
print(TJ12hi, vp = vplayout(1, 2))
 print(TJ11, vp = vplayout(1, 3))
 print(TJ12_TJ12hi_coll, vp = vplayout(2, 1))
 print(TJ12_TJ11_coll, vp = vplayout(2, 2))
 print(TJ12hi_TJ11_coll, vp = vplayout(2, 3))
#----------
 TJ12<- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12_col12_cpm + .5), y = log10(ccell4_Mature_TJ12_cpm + .5))) +
   theme_bw(base_size = 8) + geom_point(size=.1) + 
   xlab('14 cycles Collapsed log10(cpm + .5)') + ylab('14 cycles Total log10(cpm + .5)')  + coord_cartesian(xlim = c(0, 5.5), ylim=c(0,5.5), expand = FALSE)
   
 # postscript(file="Figure3D_14cycles.eps",paper = "special", width=2.5, height=2.5, horizontal=FALSE)
 print(TJ12) # 
 # dev.off()


TJ12hi <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12hiPCR_col12_cpm + .5), y = log10(ccell4_Mature_TJ12hiPCR_cpm + .5))) + 
  theme_bw(base_size = 8) + geom_point(size = .1) +  
  xlab('24 cycles Collapsed log10(cpm + .5)') + ylab('24 cycles Total log10(cpm + .5)') + coord_cartesian(xlim = c(0, 5.5), ylim=c(0,5.5), expand = FALSE)
  
# postscript(file="Figure3E_24cycles.eps",paper = "special", width=2.5, height=2.5, horizontal=FALSE)
print(TJ12hi) # 
# dev.off()

TJ12_TJ12hi_coll <- ggplot(compare_cpm, aes(x = log10(ccell4_Mature_TJ12_col12_cpm + .5), y = log10(ccell4_Mature_TJ12hiPCR_col12_cpm + .5))) + 
  theme_bw(base_size = 8) + geom_point(size=.05) +  
  xlab('14 Cycles Collapsed log10(cpm + .5)') + ylab('24 Cycles Collapsed log10(cpm + .5)') + coord_cartesian(xlim = c(0, 5.5), ylim=c(0,5.5), expand = FALSE)

# postscript(file="Figure3F_14v24.eps",paper = "special", width=2.5, height=2.5, horizontal=FALSE)
print(TJ12_TJ12hi_coll)
# dev.off()
#
# Control Cell 4 vs exosomes
# load exosome counts from the TJ11 summary tables
library('tidyverse')
theme_set(theme_minimal(base_size=8))

# TJ11raw <- read.table("~/Box Sync/SeqPaperData/TJ11/TJ11_stats/TJ11_Mature_counts_raw.txt", header = TRUE, stringsAsFactors = FALSE)
TJ11raw <- read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ11/TJ11_stats/TJ11_Mature_counts_raw.txt", header = TRUE, stringsAsFactors = FALSE)
TJ11raw[is.na(TJ11raw)] <- 0
colnames(TJ11raw) <- gsub('X4N[0-9][0-9]_', '', colnames(TJ11raw))
colnames(TJ11raw) <- gsub('_Mature_counts.txt', '_raw', colnames(TJ11raw))
rownames(TJ11raw) <- TJ11raw$miRNA
# TJ11dedup <- read.table("~/Box Sync/SeqPaperData/TJ11/TJ11_stats/TJ11_Mature_counts12.txt", header = TRUE, stringsAsFactors = FALSE)

TJ11dedup <- read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ11/TJ11_stats/TJ11_Mature_counts12.txt", header = TRUE, stringsAsFactors = FALSE)
TJ11dedup[is.na(TJ11dedup)] <- 0
colnames(TJ11dedup) <- gsub('X4N[0-9][0-9]_', '', colnames(TJ11dedup))
rownames(TJ11dedup) <- TJ11dedup$miRNA

TJ11compare <- inner_join(TJ11dedup, TJ11raw, by = 'miRNA')
# from Jake
library('edgeR')
library('limma')
library('stringr')
library('ggplot2')
theme_set(theme_classic())

diff <- function(x) {
  count_mat <- as.matrix(x)
  layout <- data.frame(row.names = colnames(count_mat),
                       condition = rep(0, length(colnames(count_mat))))
  layout[grep('cell', rownames(layout)), 'condition'] <- 'cell'
  layout[grep('exo', rownames(layout)), 'condition'] <- 'exo'
  print(layout)
  
  count_mat_filtered <- count_mat[rowSums(count_mat > 20) >= 3, ]
  
  f <- factor(layout$condition)
  design <- model.matrix( ~ f)
  rownames(design) <- rownames(layout)
  
  # Normalize to log2(CPM + 0.5) with voom
  d <- DGEList(count_mat_filtered, genes = rownames(count_mat_filtered))
  v <- voom(d, design, plot = T)
  
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  results <- topTable(fit, number = nrow(count_mat_filtered))
  return(results)
}

dedup_results <- diff(dplyr::select(TJ11dedup, matches('ctrl_cell|ctrl_exoF')))
sum(dedup_results$adj.P.Val < 0.05)
ggplot(dedup_results, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point() +
  geom_point(data = dplyr::filter(dedup_results, adj.P.Val < 0.05), aes(x = logFC, y = -log10(adj.P.Val)), col = 'red') +
  xlab('log fold change') + ylab('-log10 (adjusted p value)') +
  ggtitle('Deduplicated')

Nodedup_results <- diff(dplyr::select(TJ11raw, matches('ctrl_cell|ctrl_exoF')))
sum(Nodedup_results$adj.P.Val < 0.05)
ggplot(Nodedup_results, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point() +
  geom_point(data = dplyr::filter(Nodedup_results, adj.P.Val < 0.05), aes(x = logFC, y = -log10(adj.P.Val)), col = 'red') +
  xlab('log fold change') + ylab('-log10 (adjusted p value)') +
  ggtitle('Raw data')

combined <- inner_join(dedup_results, Nodedup_results, by = 'genes')
FC_CellvExo<- ggplot(combined, aes(x = logFC.x, y = logFC.y)) + geom_point(size=.1) +theme_bw(base_size = 8) +
  xlab('logFC.collapsed') + ylab('logFC.raw')


# postscript(file="Figure3Fv2.eps",paper = "special", width=2.5, height=2, horizontal=FALSE)
print(FC_CellvExo) # 
# dev.off()
#

# system information at time of this rendering:
sessionInfo()