setwd("~/Box Sync/SeqPaperDataCB/RAnalysis")
# analysis of Mature counts of control cell 4 in three libraries TJ12, TJ12hiPCR and TJ11
# analysis of coefficient of variation of serum wiht multiple repeats in Tj12 and TJ12hiPCR

# this differs from Figure4.r because all cpm is taken out.  also did full join of tables, then eliminated miRNAs with low expression
library('tidyverse')
theme_set(theme_classic(base_size=8))

TJ12raw <- read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ12/TJ12_stats/TJ12_Mature_counts_raw.txt", header=TRUE, stringsAsFactors = FALSE)
TJ12raw[is.na(TJ12raw)] <- 0
colnames(TJ12raw) <- gsub("X4N[0-9][0-9]_", "TJ12raw_", colnames(TJ12raw))
# TJ12raw_cpm <- TJ12raw %>% select(-miRNA) %>% mutate_all(., (funs(cpm=(./sum(.)*10000)))) %>% select(ends_with("_cpm")) %>% mutate(miRNA=TJ12raw$miRNA)

TJ12col <-  read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ12/TJ12_stats/TJ12_Mature_counts12.txt", header=TRUE, stringsAsFactors = FALSE)
TJ12col[is.na(TJ12col)] <- 0
colnames(TJ12col) <- gsub("X4N[0-9][0-9]_", "TJ12col_", colnames(TJ12col))

TJ12hiraw <-  read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ12hi/TJ12hi_stats/TJ12hi_Mature_counts_raw.txt", header=TRUE, stringsAsFactors = FALSE)
TJ12hiraw[is.na(TJ12hiraw)] <- 0
colnames(TJ12hiraw) <- gsub("X4N[0-9][0-9]_", "TJ12hiraw_", colnames(TJ12hiraw))

TJ12hicol <-  read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ12hi/TJ12hi_stats/TJ12hi_Mature_counts12.txt", header=TRUE, stringsAsFactors = FALSE)
TJ12hicol[is.na(TJ12hicol)] <- 0
colnames(TJ12hicol) <- gsub("X4N[0-9][0-9]_", "TJ12hicol_", colnames(TJ12hicol))


TJ11raw <-  read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ11/TJ11_stats/TJ11_Mature_counts_raw.txt", header=TRUE, stringsAsFactors = FALSE)
TJ11raw[is.na(TJ11raw)] <- 0
colnames(TJ11raw) <- gsub("X4N[0-9][0-9]_", "TJ11raw_", colnames(TJ11raw))

TJ11col <-  read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ11/TJ11_stats/TJ11_Mature_counts12.txt", header=TRUE, stringsAsFactors = FALSE)
TJ11col[is.na(TJ11col)] <- 0
colnames(TJ11col) <- gsub("X4N[0-9][0-9]_", "TJ11col_", colnames(TJ11col))

setwd("~/Box Sync/SeqPaperDataCB/RAnalysis/PCRcycle_cell4")

# ---- inner join ccell4 --------
compare_counts <- list(TJ12raw, TJ12col, TJ12hiraw, TJ12hicol, TJ11raw, TJ11col)  %>% 
  Reduce( function(dtf1,dtf2) full_join(dtf1,dtf2, by = "miRNA"), .)
compare_counts[is.na(compare_counts)] <-0

# keep rows with >5 and exclude the miRNA column
# counts_filtered <- counts[rowSums(counts>5) >= 3 , ]
compare_ccell4 <- compare_counts %>% select(miRNA, contains("ccell4_500ng_1"), contains("ctrl_cell_4")) %>% 
   arrange(desc(TJ12col_ccell4_500ng_1))
compare_ccell4_filtered <- compare_ccell4[rowSums(compare_ccell4 > 10) >= 3, ]
cutoff <- 10
library(forcats)   # to convert miRNAs into factors
Unique_ccell4 <- compare_ccell4_filtered %>% mutate(PercentUnique_TJ12 =TJ12col_ccell4_500ng_1/TJ12raw_ccell4_500ng_1*100, 
         PercentUnique_TJ12hi = TJ12hicol_ccell4_500ng_1/TJ12hiraw_ccell4_500ng_1*100, 
         PercentUnique_TJ11 = TJ11col_ctrl_cell_4/TJ11raw_ctrl_cell_4*100) %>% 
  select(miRNA, starts_with("PercentUnique_")) %>% arrange(desc(PercentUnique_TJ12))
Unique_ccell4$miRNA<-parse_factor(Unique_ccell4$miRNA, levels=fct_inorder(Unique_ccell4$miRNA))


TJ12_TJ12hi_pc3 <- ggplot(Unique_ccell4) + geom_point(aes(x = miRNA, y = PercentUnique_TJ12), size=.5) +
  xlab('miRNA') + ylab("Percent Unique (collapsed/raw * 100)") + ggtitle('Percent Unique Comparison') +
  labs(subtitle="14 Cycles (black)*; 16 Cycles (red); 24 Cycles (blue)") +
  geom_point(aes(x=miRNA, y=PercentUnique_TJ12hi), color = "blue", size = .5)  +
  geom_point(aes(x=miRNA, y=PercentUnique_TJ11), color = "red", size=.5)  +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) 
print(TJ12_TJ12hi_pc3)

# for publication: # labs(subtitle="14 Cycles (black)*; 16 Cycles (red); 24 Cycles (blue)") +
TJ12_TJ12hi_pc3 <- ggplot(Unique_ccell4) + geom_point(aes(x = miRNA, y = PercentUnique_TJ12), size=.05) +
  xlab('miRNA') + ylab("Percent Unique") +
  geom_point(aes(x=miRNA, y=PercentUnique_TJ12hi), color = "blue", size = .05)  +
  geom_point(aes(x=miRNA, y=PercentUnique_TJ11), color = "red", size=.05)  +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) 
# postscript(file="Figure4c_v3.eps",paper = "special", width=2.5, height=2.5, horizontal = FALSE)
print(TJ12_TJ12hi_pc3)
 # dev.off()

# ----------------
library(forcats)   # to convert miRNAs into factors
compare_ccell4_sortcol <- compare_ccell4_filtered %>% arrange(desc(TJ12col_ccell4_500ng_1))
compare_ccell4_sortcol$miRNA<-parse_factor(compare_ccell4_sortcol$miRNA, levels=fct_inorder(compare_ccell4_sortcol$miRNA))
compare_ccell4_sortraw <- compare_ccell4_filtered %>% arrange(desc(TJ12raw_ccell4_500ng_1))
compare_ccell4_sortraw$miRNA<-parse_factor(compare_ccell4_sortraw$miRNA, levels=fct_inorder(compare_ccell4_sortraw$miRNA))
# -----

# figure 4a now figure 3G
red_blueraw2 <- ggplot(compare_ccell4_sortraw) + labs(title="Cycle Number Impact on Raw Counts") + 
  theme_classic(base_size = 8)  +
  geom_point(aes(x=miRNA, y=log10(TJ12raw_ccell4_500ng_1+.5)), color = "black", size=.5) +
  geom_point(aes( x=miRNA, y=log10(TJ12hiraw_ccell4_500ng_1+.5)), color="blue", size=.5) + 
  geom_point(aes( x=miRNA, y=log10(TJ11raw_ctrl_cell_4+.5)), color="red",size=.5) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  ylim(-0.5,5) +
  ylab("Log10 Raw Counts + 0.5") +
  labs(subtitle="14 Cycles (black)*; 16 Cycles (red); 24 Cycles (blue)")
# * ranked on 14 cycles (TJ12)    

print(red_blueraw2)

# for publication: # labs(subtitle="14 Cycles (black)*; 16 Cycles (red); 24 Cycles (blue)")
# * ranked on 14 cycles (TJ12)
red_blueraw2 <- ggplot(compare_ccell4_sortraw) + 
  theme_classic(base_size = 8)  +
  geom_point(aes(x=miRNA, y=log10(TJ12raw_ccell4_500ng_1+.5)), color = "black", size=.05) +
  geom_point(aes( x=miRNA, y=log10(TJ12hiraw_ccell4_500ng_1+.5)), color="blue", size=.05) + 
  geom_point(aes( x=miRNA, y=log10(TJ11raw_ctrl_cell_4+.5)), color="red",size=.05) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  coord_cartesian(ylim=c(0,5.5), expand = FALSE) +
  ylab("Total log10(counts + 0.5)") 
  
# postscript(file="Figure3G_v3.eps",paper = "special", width=2.5, height=2.5, horizontal = FALSE)
print(red_blueraw2)
# dev.off()

# figure 4b  now figure 3H
red_bluecollapsed2 <- ggplot(compare_ccell4_sortcol) + labs(title="Cycle Number Impact on Collapsed Counts") + 
  theme_classic(base_size = 8)  +
  geom_point(aes(x=miRNA, y=log10(TJ12col_ccell4_500ng_1+.5)), color = "black", size=.5) +
  geom_point(aes( x=miRNA, y=log10(TJ12hicol_ccell4_500ng_1+.5)), color="blue", size=.5) + 
  geom_point(aes( x=miRNA, y=log10(TJ11col_ctrl_cell_4+.5)), color="red",size=.5) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  ylab("Log10 Collapsed Counts + 0.5") +
  ylim(-0.5,5) +
  labs(subtitle="14 Cycles (black)*; 16 Cycles (red); 24 Cycles (blue)")

print(red_bluecollapsed2)

# for publication: # labs(subtitle="14 Cycles (black)*; 16 Cycles (red); 24 Cycles (blue)")
red_bluecollapsed2 <- ggplot(compare_ccell4_sortcol) + 
  theme_classic(base_size = 8)  +
  geom_point(aes(x=miRNA, y=log10(TJ12col_ccell4_500ng_1+.5)), color = "black", size=.05) +
  geom_point(aes( x=miRNA, y=log10(TJ12hicol_ccell4_500ng_1+.5)), color="blue", size=.05) + 
  geom_point(aes( x=miRNA, y=log10(TJ11col_ctrl_cell_4+.5)), color="red",size=.05) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  ylab("Collapsed log10(counts + 0.5)") +
   coord_cartesian(ylim=c(0,5.5), expand = FALSE) 
  
# * ranked on 14 cycles (TJ12)  collapsed
# postscript(file="Figure3H_v3.eps",paper = "special", width=2.5, height=2.5, horizontal = FALSE)
print(red_bluecollapsed2)
# dev.off()


# ------set up for Figure 4 E -------------
# Use TJ12 and TJ12hi libraries  
# look at all CB2 samples (adapters 23-30), CB1 (adapters 19-22), 



# ---- inner join --------
# compare_counts <- list(TJ12raw, TJ12col, TJ12hiraw, TJ12hicol) %>% Reduce( function(dtf1,dtf2) inner_join(dtf1,dtf2, by = "miRNA"), .)
# compare_counts[is.na(compare_counts)] <-0
# library(forcats)   # to convert miRNAs into factors

# compare_counts12 <- compare_counts %>% select(contains("CB2_1f_")) # %>% arrange(desc(TJ12col_CB2_1f_1))

compare_CB2 <- compare_counts %>% select(miRNA, contains("CB2_1f_")) %>% arrange(desc(TJ12col_CB2_1f_1))
compare_CB2_filtered <- compare_CB2[rowSums(compare_CB2 > 10) >= 3, ]
compare_CB2col <- compare_CB2_filtered %>% arrange(desc(TJ12col_CB2_1f_1))
compare_CB2col$miRNA<-parse_factor(compare_CB2col$miRNA, levels=fct_inorder(compare_CB2col$miRNA))
compare_CB2raw <- compare_CB2_filtered %>% arrange(desc(TJ12raw_CB2_1f_1))
compare_CB2raw$miRNA<-parse_factor(compare_CB2raw$miRNA, levels=fct_inorder(compare_CB2raw$miRNA))


 # Figure 4D
 red_blueCB2raw2 <- ggplot(compare_CB2raw) + labs(title="Cycle Number Impact on Raw Counts\nSerum Duplicates") + 
   theme_classic(base_size = 8)  +
   geom_point(aes(x=miRNA, y=log10(TJ12raw_CB2_1f_1+.5)), color = "black", size=.5) +
   geom_point(aes( x=miRNA, y=log10(TJ12raw_CB2_1f_2+.5)), color="black", shape=21, size=.5, stroke=.5) + 
   geom_point(aes( x=miRNA, y=log10(TJ12hiraw_CB2_1f_1hi+.5)), color="blue", size=.5) +
   geom_point(aes( x=miRNA, y=log10(TJ12hiraw_CB2_1f_2hi+.5)), color="blue", shape=21, size=.5, stroke=.5) +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
   ylab("Log10 (Raw Counts +.5)") +
   labs(subtitle="14 cycles 1-black*, 2-black-open; 24 cycles 1-blue; 2-blue-open")
 #  postscript(file="Figure4d.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print(red_blueCB2raw2)
 # dev.off()
 
 # Figure 4E
 red_blueCB2col2 <- ggplot(compare_CB2col) + labs(title="Cycle Number Impact on Collapsed Counts\nSerum Duplicates") + 
   theme_classic(base_size = 8)  + 
   geom_point(aes(x=miRNA, y=log10(TJ12col_CB2_1f_1+.5)), color = "black", size=.5) +
   geom_point(aes( x=miRNA, y=log10(TJ12col_CB2_1f_2+.5)), color="black", shape=21, size=.5, stroke=.5) + 
   geom_point(aes( x=miRNA, y=log10(TJ12hicol_CB2_1f_1hi+.5)), color="blue", size=.5) +
   geom_point(aes( x=miRNA, y=log10(TJ12hicol_CB2_1f_2hi+.5)), color="blue", shape=21, size=.5, stroke=.5) +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
   ylab("Log10 (Raw Counts +.5)") +
   ylim(-0.5,5) +
   labs(subtitle="14 cycles 1-black*, 2-black-open; 24 cycles 1-blue; 2-blue-open")
#  postscript(file="Figure4e.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print(red_blueCB2col2)
#  dev.off()
 


# -- Coefficient of Variation ------
 
 # coefficient of variation for variable = 100 * sd / mean
 # coefficient of variatin for a model = 100 * e(rmse) / r(mean)
 
 # CV <- function(mean, sd){
 #  (sd/mean)*100
 # }
 #
 # use cv function of raster package
 # install.packages("raster")
# library("raster")
 
 CB2meansTJ12raw <- rowMeans(compare_CB2col[,2:5])
 CB2meansTJ12coll <- rowMeans(compare_CB2col[,6:9])
 CB2meansTJ12.2raw <- rowMeans(compare_CB2col[,10:13])
 CB2meansTJ12.2coll <- rowMeans(compare_CB2col[,14:17])
 
 CB2_CVTJ12raw <- ( (apply(compare_CB2col[,2:5],1, function(x) {sd(x)})) / (rowMeans(compare_CB2col[,2:5])) * 100 )
 CB2_CVTJ12coll <- ( (apply(compare_CB2col[,6:9],1, function(x) {sd(x)})) / (rowMeans(compare_CB2col[,6:9])) * 100 )
 CB2_CVTJ12.2raw <- ( (apply(compare_CB2col[,10:13],1, function(x) {sd(x)})) / (rowMeans(compare_CB2col[,10:13])) * 100 )
 CB2_CVTJ12.2coll <- ( (apply(compare_CB2col[,14:17],1, function(x) {sd(x)})) / (rowMeans(compare_CB2col[,14:17])) * 100 )
 
 # the cv function from the raster package gives exactly the same number as teh above calculation
# CB2_CVTJ12.2col_rastercv <- apply(compare_CB2col[,14:17],1, function(x) {raster::cv(x)})
 
 df_CB2 <- data_frame(miRNA=compare_CB2col$miRNA, CB2meansTJ12raw, CB2meansTJ12coll, CB2meansTJ12.2coll, CB2_CVTJ12raw, CB2_CVTJ12coll, CB2_CVTJ12.2raw, CB2_CVTJ12.2coll)
 df_CB2[is.na(df_CB2)] <-0
 df_CB2col<- df_CB2 %>% filter(CB2meansTJ12coll>1) %>% arrange(CB2meansTJ12coll)
 df_CB2col$miRNA<-parse_factor(df_CB2col$miRNA, levels=fct_inorder(df_CB2col$miRNA))
 df_CB2col$Rank <- rank(df_CB2col$CB2meansTJ12coll, ties.method = "first")
 df_CB2col$RankHi <- rank(df_CB2col$CB2meansTJ12.2coll, ties.method = "first")
 #


 # just 14 cycle samples
 
 CV_CB2plot14 <- ggplot(df_CB2col) + labs(title="Coefficient of Variation of Serum CB2") + 
   
   geom_point(aes( x=Rank, y=CB2_CVTJ12coll), color="black", shape=21, size=.5, stroke=.5) + 
   geom_point(aes(x=Rank, y=CB2_CVTJ12raw), color = "blue", size=.5) +
   geom_smooth(aes( x=Rank, y=CB2_CVTJ12coll), colour="black", position = "identity", se=FALSE, method="loess") +
   geom_smooth(aes( x=Rank, y=CB2_CVTJ12raw), colour="blue", position = "identity", se=FALSE, method="loess") +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
   ylab("Coefficient of Variation") +
   xlab("miRNAs Ranked by Mean Expression with 14 cycles") +
   labs(subtitle="14 cycles raw-blue, collapsed*-black")
 
 # postscript(file="Figure4f.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print(CV_CB2plot14)
  # dev.off()

 CV_CB2plot24 <- ggplot(df_CB2col) + labs(title="Coefficient of Variation of Serum CB2") + 
   
   geom_point(aes( x=RankHi, y=CB2_CVTJ12.2coll), color="black", shape=21, size=.5, stroke=.5) + 
   geom_point(aes(x=RankHi, y=CB2_CVTJ12.2raw), color = "blue", size=.5) +
   geom_smooth(aes( x=RankHi, y=CB2_CVTJ12.2coll), colour="black", position = "identity", se=FALSE, method="loess") +
   geom_smooth(aes( x=RankHi, y=CB2_CVTJ12.2raw), colour="blue", position = "identity", se=FALSE, method="loess") +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
   ylab("Coefficient of Variation") +
   xlab("miRNAs Ranked by Mean Expression with 24 cycles") +
   labs(subtitle="24 cycles raw-blue, collapsed* -black")
 
 # postscript(file="Figure4g.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print(CV_CB2plot24)
  # dev.off()

 CV_CB2plot_16v24col <- ggplot(df_CB2col) + labs(title="Coefficient of Variation of Serum CB2") + 
   
   geom_point(aes( x=Rank, y=CB2_CVTJ12coll), color="black",  shape=21, size=.5, stroke=.5) + 
   geom_point(aes(x=Rank, y=CB2_CVTJ12.2coll), color = "blue", shape=21, size=.5, stroke=.5) +
   geom_smooth(aes( x=Rank, y=CB2_CVTJ12coll), colour="black", position = "identity", se=FALSE, method="loess") +
   geom_smooth(aes( x=Rank, y=CB2_CVTJ12.2coll), colour="blue", position = "identity", se=FALSE, method="loess") +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
   ylab("Coefficient of Variation") +
   xlab("miRNAs Ranked by Mean Expression with 14 cycles") +
   labs(subtitle="14 cycles Collapsed* -black, 24 cycles Collapsed -blue")
 
#  postscript(file="Figure4Hopened.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print(CV_CB2plot_16v24col)
 # dev.off()
 
 
 
 # there is essentially no difference in how the graph looks.
 #
 
 
 # -- Coefficient of Variation include more CB2 samples------
 
 compare_CB2all <- compare_counts %>% select(miRNA, contains("CB2")) %>% arrange(desc(TJ12col_CB2_1f_1))
 compare_CB2all_filtered <- compare_CB2all[rowSums(compare_CB2all > 10) >= 4, ]
 compare_CB2all_sortcol <- compare_CB2all_filtered %>% arrange(desc(TJ12col_CB2_1f_1))
 compare_CB2all_sortcol$miRNA<-parse_factor(compare_CB2all_sortcol$miRNA, levels=fct_inorder(compare_CB2all_sortcol$miRNA))

 # coefficient of variation for variable = 100 * sd / mean
 # coefficient of variatin for a model = 100 * e(rmse) / r(mean)
 
 # CV <- function(mean, sd){
 #  (sd/mean)*100
 # }
 # CB2allmeans <- compare_CB2all_sortcol %>% mutate_at(dplyr::vars(c(2:5)), dplyr::funs(round))
 CB2meansTJ12raw <- rowMeans(compare_CB2all_sortcol[,2:9])
 CB2meansTJ12coll <- rowMeans(compare_CB2all_sortcol[,10:17])
 CB2meansTJ12.2raw <- rowMeans(compare_CB2all_sortcol[,18:25])
 CB2meansTJ12.2coll <- rowMeans(compare_CB2all_sortcol[,26:33])
 
 CB2_CVTJ12raw <- ( (apply(compare_CB2all_sortcol[,2:9],1, function(x) {sd(x)})) / (rowMeans(compare_CB2all_sortcol[,2:9])) * 100 )
 CB2_CVTJ12coll <- ( (apply(compare_CB2all_sortcol[,10:17],1, function(x) {sd(x)})) / (rowMeans(compare_CB2all_sortcol[,10:17])) * 100 )
 CB2_CVTJ12.2raw <- ( (apply(compare_CB2all_sortcol[,18:25],1, function(x) {sd(x)})) / (rowMeans(compare_CB2all_sortcol[,18:25])) * 100 )
 CB2_CVTJ12.2coll <- ( (apply(compare_CB2all_sortcol[,26:33],1, function(x) {sd(x)})) / (rowMeans(compare_CB2all_sortcol[,26:33])) * 100 )
 
 df_CB2 <- data_frame(miRNA=compare_CB2all_sortcol$miRNA, CB2meansTJ12raw, CB2meansTJ12coll, CB2_CVTJ12raw, CB2_CVTJ12coll, CB2_CVTJ12.2raw, CB2_CVTJ12.2coll)
 df_CB2[is.na(df_CB2)] <-0
 df_CB2col<- df_CB2 %>% filter(CB2meansTJ12coll>1) %>% arrange(desc(CB2meansTJ12coll))
 df_CB2col$miRNA<-parse_factor(df_CB2col$miRNA, levels=fct_inorder(df_CB2col$miRNA))
 #

 
 CV_CB2plot <- ggplot(df_CB2col) + labs(title="Coefficient of Variation of Serum CB2 (8 samples per group)") + 
   geom_point(aes(x=miRNA, y=CB2_CVTJ12raw), color = "black", size=.5) +
   geom_point(aes( x=miRNA, y=CB2_CVTJ12coll), color="blue", size=.5) + 
   geom_point(aes( x=miRNA, y=CB2_CVTJ12.2raw), color="red", size=.5) +
   geom_point(aes( x=miRNA, y=CB2_CVTJ12.2coll), color="green", size=.5) +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
   ylab("Coefficient of Variation") +
   labs(subtitle="14 cycles raw-black, collapsed-blue*; 24 cycles raw-red, collapsed-green\n* miRNAs ranked on this sample")
 
 print(CV_CB2plot)
 
 #' ------- Unique comparison between cell and exo in TJ11 -----------
 #' 
 #' # keep rows with >5 and exclude the miRNA column
 # counts_filtered <- counts[rowSums(counts>1) >= 12 , ]  REMOVE ALL COUNTS=0
 
 compare_ccell_exo <- compare_counts %>% select(miRNA, contains("_ctrl_cell_"), contains("_ctrl_exoF_")) %>% 
   arrange(desc(TJ11col_ctrl_cell_4))
 compare_ccell_exo_filtered <- compare_ccell_exo[rowSums(compare_ccell_exo > 1) >= 12, ]
 cutoff <- 10
 library(forcats)   # to convert miRNAs into factors
 Unique_ccell_exo <- compare_ccell_exo_filtered %>% mutate(PercentUnique_2cell = TJ11col_ctrl_cell_2/TJ11raw_ctrl_cell_2*100, 
                                                        PercentUnique_2exo = TJ11col_ctrl_exoF_2/TJ11raw_ctrl_exoF_2*100, 
                                                        PercentUnique_3cell = TJ11col_ctrl_cell_3/TJ11raw_ctrl_cell_3*100, 
                                                        PercentUnique_3exo = TJ11col_ctrl_exoF_3/TJ11raw_ctrl_exoF_3*100, 
                                                     PercentUnique_4cell = TJ11col_ctrl_cell_4/TJ11raw_ctrl_cell_4*100, 
                                                     PercentUnique_4exo = TJ11col_ctrl_exoF_4/TJ11raw_ctrl_exoF_4*100) %>% 
   select(miRNA, starts_with("PercentUnique_"))  %>% arrange(desc(PercentUnique_4cell))
 Unique_ccell_exo$miRNA<-parse_factor(Unique_ccell_exo$miRNA, levels=fct_inorder(Unique_ccell_exo$miRNA))
 
 PCuniq_cell_exo <- ggplot(Unique_ccell_exo) + geom_point(aes(x = miRNA, y = PercentUnique_2cell), size=.5) +
   xlab('miRNA') + ylab("Percent Unique (collapsed/raw * 100)") + ggtitle('Percent Unique Comparison') +
   geom_point(aes(x=miRNA, y=PercentUnique_2exo), color = "blue", size = .5)  +
   geom_point(aes(x=miRNA, y=PercentUnique_3cell), color = "gray", size=.5)  +
   geom_point(aes(x=miRNA, y=PercentUnique_3exo), color = "purple", size=.5)  +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) 
 
 PCuniq_cellvexo <- ggplot(Unique_ccell_exo)  +  ggtitle('Percent Unique Comparison') +
   ylab("Percent Unique (collapsed/raw * 100)") +
   geom_point(aes(x = PercentUnique_2cell, y = PercentUnique_2exo), size=.5, color = "gray") +
   geom_point(aes(x = PercentUnique_3cell, y = PercentUnique_3exo), size=.5, color = "blue") +
   geom_point(aes(x = PercentUnique_4cell, y = PercentUnique_4exo), size=.5, color = "black") 
 
#  postscript(file="Figure3x3.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print( PCuniq_cellvexo)
 #  dev.off()
 
 # take average of cells and exo then calculate percent
 mean_unique_cell_exo <- compare_ccell_exo_filtered %>% 
   mutate(ColCellMean = rowMeans(select(compare_ccell_exo_filtered, starts_with("TJ11col_ctrl_cell")), na.rm = TRUE),
          RawCellMean = rowMeans(select(compare_ccell_exo_filtered, starts_with("TJ11raw_ctrl_cell")), na.rm = TRUE),
          ColExoMean = rowMeans(select(compare_ccell_exo_filtered, starts_with("TJ11col_ctrl_exoF")), na.rm = TRUE),
          RawExoMean = rowMeans(select(compare_ccell_exo_filtered, starts_with("TJ11raw_ctrl_exoF")), na.rm = TRUE)) %>%
   transmute(miRNA, ColCelluniq=ColCellMean/RawCellMean*100, ColExouniq=ColExoMean/RawExoMean*100) %>%
   select(miRNA, contains("uniq"))  %>% arrange(desc(ColCelluniq))
 mean_unique_cell_exo$miRNA<-parse_factor(mean_unique_cell_exo$miRNA, levels=fct_inorder(mean_unique_cell_exo$miRNA))

 MeanPCuniq_cell_exo <- ggplot(mean_unique_cell_exo) + geom_point(aes(x = miRNA, y = ColCelluniq), size=.5) +
   xlab('miRNA') + ylab("Percent Unique (collapsed/raw * 100)") + ggtitle('Percent Unique Comparison') +
   geom_point(aes(x=miRNA, y=ColExouniq), color = "blue", size = .5) +
   theme(axis.text.x = element_blank(),axis.ticks = element_blank()) 
 print(MeanPCuniq_cell_exo)
 
 
 fit_cellvexo <- lm(ColCelluniq ~ ColExouniq, data=mean_unique_cell_exo)
 summary(fit_cellvexo)
 
 MeanPCuniq_cellvexo <- ggplot(mean_unique_cell_exo)  +  ggtitle('Percent Unique Cells vs Exo') +
   ylab("Percent Unique Exosome (collapsed/raw * 100)") +
   xlab("Percent Unique Cell (collapsed/raw * 100)") +
   geom_point(aes(x = ColCelluniq, y = ColExouniq), size=.5)+
   coord_cartesian(xlim = c(15,80), ylim = c(15,80), expand = TRUE) + 
   geom_smooth(aes(x = ColCelluniq, y = ColExouniq),method = "lm", se = FALSE, show.legend = TRUE) +
   annotate("text", x=25, y=75, label = "R^2=0.8492") 

 # postscript(file="Figure3G.eps",paper = "special", width=4, height=3, horizontal = FALSE)
 print(MeanPCuniq_cellvexo)
#  dev.off()
 