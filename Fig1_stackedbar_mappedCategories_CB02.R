
library('tidyverse')
theme_set(theme_minimal())
basedir <- "/Volumes/BiomarkerSeq/SeqPaperData/CB02/CB02_stats"
setwd(basedir)

# read noN data files, combine contaminant specific data with summary data for noN into summary table
# then overwrite the raw data with the next dataset only keeping the combined data in memory

Summarydata <- read.table("Map_stats_CB02_noN_all.txt", header = TRUE, stringsAsFactors = FALSE)
# now break down the UMI contaminant data into markers and five prime adapter --------------
ContamSpecificdata <- read.table("CB02_Libcontam_counts_raw.txt", header=TRUE, stringsAsFactors = FALSE)
# replace NA with 0
ContamSpecificdata[is.na(ContamSpecificdata)] <- 0
names(ContamSpecificdata) <- gsub("X4", "4", names(ContamSpecificdata))
rownames(ContamSpecificdata) <- ContamSpecificdata$miRNA

MostContamData <- filter(ContamSpecificdata, miRNA=="19nt_Marker" | miRNA == "24nt_Marker" | miRNA == "5primeAdapter" ) # choose just marker and 5 prime adapter and not phiX data
rownames(MostContamData) <- MostContamData$miRNA
MostContamDataT<-data.frame(t(MostContamData[,2:51]))
contamSamples <- rownames(MostContamDataT)
MostContamDataT <- MostContamDataT %>% mutate(other_contam = Summarydata$LibraryContam-(X19nt_Marker+X24nt_Marker+X5primeAdapter)) %>% mutate(Sample = contamSamples)
Samplelist <- Summarydata$Sample
# combine the count data
completeSummary_NoUMI <- right_join(Summarydata, MostContamDataT, by = "Sample") %>% select(-LibraryContam, -Hairpins1, -DemultiCount) %>% 
  mutate(Markers=X19nt_Marker+X24nt_Marker, otherSmallRNA=Ribo+tRNA+svRNA) %>%
  select(-X19nt_Marker,-X24nt_Marker,-Ribo,-tRNA,-svRNA) %>% rename(Adapter=X5primeAdapter) %>%
  mutate(Unmapped = (TrimCount-(Markers+Adapter+other_contam+Mature1+otherSmallRNA+Genomic))) 

PercentDataComplete_NoUMI <- completeSummary_NoUMI %>% select(-(Sample)) %>% mutate_all(funs(pc=(./TrimCount*100))) %>% select(ends_with("_pc"), -TrimCount_pc) %>%  mutate(Sample = Samplelist)

 PercentDatalong_NoUMI <- gather(PercentDataComplete_NoUMI, MapCategory, PercentTotal, Markers_pc, Adapter_pc, other_contam_pc, Mature1_pc, otherSmallRNA_pc, Genomic_pc, Unmapped_pc) 
# use unite() to combine the sample name and map category into one variable/column
# PClong_NoUMI<- unite(PercentDatalong_NoUMI, Sample.Category,  Sample, MapCategory, sep=".", remove=TRUE)

# Countslong_noUMI <- completeSummary_NoUMI %>% gather(MapCategory, Counts, Markers, Adapter,other_contam, Mature1, otherSmallRNA, Genomic, Unmapped) %>% unite( Sample.Category,  Sample, MapCategory, sep=".", remove=TRUE)

rm(Summarydata)
rm(ContamSpecificdata)
rm(MostContamData)
rm(MostContamDataT)

library(RColorBrewer)
mypalette<-brewer.pal(11,"Paired")

library(forcats)   # to convert miRNAs into factors
# compare_counts12 <- compare_counts %>% arrange(desc(TJ12coll_CB2_1f_1))
# compare_counts12$Mature<-parse_factor(compare_counts12$Mature, levels=fct_inorder(compare_counts12$Mature))
# PercentDatalong_NoUMI$MapCategory <- parse_factor(PercentDatalong_NoUMI$MapCategory, levels=fct_inorder(PercentDatalong_NoUMI$MapCategory))
PercentDatalong_NoUMI$MapCategory <- parse_factor(PercentDatalong_NoUMI$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Mature1_pc","other_contam_pc","Adapter_pc","Markers_pc"))

noNbargraph <- ggplot(PercentDatalong_NoUMI) + theme_bw()+ ggtitle("Percent Category Mapped No Collapse")
noNbargraph + geom_bar(aes(y = PercentTotal, x = Sample, fill = MapCategory), stat="identity") + theme(axis.text.x  = element_text(angle=90, vjust=.5, size=8)) +
  scale_fill_manual(values=mypalette)
# modified 5p adapter 100, 3p 17 serum = 7:9, exo 16:18, cells 
# unmodified 5p adapter 100, 3p 17 serum 46:48
ModvUnMod5p<- gather(PercentDataComplete_NoUMI[c(7:9,46:48),], MapCategory, PercentTotal, Markers_pc, Adapter_pc, other_contam_pc, Mature1_pc, otherSmallRNA_pc, Genomic_pc, Unmapped_pc) 
ModvUnMod5p$MapCategory <- parse_factor(ModvUnMod5p$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Mature1_pc","other_contam_pc","Adapter_pc","Markers_pc"))

ModvUnMod5pgraph <- ggplot(ModvUnMod5p) + theme_bw()+ ggtitle("Percent Category Mapped")
ModvUnMod5pgraph + geom_bar(aes(y = PercentTotal, x = Sample, fill = MapCategory), stat="identity") + theme(axis.text.x  = element_text(angle=90, vjust=.5, size=8)) +
  scale_fill_manual(values=mypalette)
# 

Mod3ptitSerum<- gather(PercentDataComplete_NoUMI[1:9,], MapCategory, PercentTotal, Markers_pc, Adapter_pc, other_contam_pc, Mature1_pc, otherSmallRNA_pc, Genomic_pc, Unmapped_pc) 
Mod3ptitSerum$MapCategory <- parse_factor(Mod3ptitSerum$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Mature1_pc","other_contam_pc","Adapter_pc","Markers_pc"))

Mod3ptitSerumgraph <- ggplot(Mod3ptitSerum) + theme_bw()+ ggtitle("Percent Mapped per Category")
Mod3ptitSerumgraph + geom_bar(aes(y = PercentTotal, x = Sample, fill = MapCategory), stat="identity") + theme(axis.text.x  = element_text(angle=90, vjust=.5, size=8)) +
  scale_fill_manual(values=mypalette)
# ---- this one changed many times to save various serum graphs----
NoMod3ptitSerum<- gather(PercentDataComplete_NoUMI[c(19:21,30:32,40:42),], MapCategory, PercentTotal, Markers_pc, Adapter_pc, other_contam_pc, Mature1_pc, otherSmallRNA_pc, Genomic_pc, Unmapped_pc) 
NoMod3ptitSerum$MapCategory <- parse_factor(NoMod3ptitSerum$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Mature1_pc","other_contam_pc","Adapter_pc","Markers_pc"))

NoMod3ptitSerumgraph <- ggplot(NoMod3ptitSerum) + theme_bw()+ ggtitle("Percent Mapped per Category")
NoMod3ptitSerumgraph + geom_bar(aes(y = PercentTotal, x = Sample, fill = MapCategory), stat="identity") + theme(axis.text.x  = element_text(angle=90, vjust=.5, size=8)) +
  scale_fill_manual(values=mypalette)


# ---- stacked dot plot various serum graphs----

# 5p titration with lowest 3p concentration 19:21, 30:32, 40:42 
# PercentDataComplete_NoUMI[c(19:21,30:32,40:42),c(8,1)] # sample = column 8, mature miRNA = column 1
# Serum<- PercentDataComplete_NoUMI[c(19:21,30:32,40:42),c(8,1)] %>% mutate(Concentration=c(rep(10,3), rep(30,3), rep(100,3)))

# Serumgraph <- ggplot(Serum) + theme_bw(base_size=8)

# Serumgraph + geom_point(aes(y = Mature1_pc, x = Concentration)) +
#  labs(x="5' Adapter Concentration", y="Percent Mature miRNAs") +
#  coord_cartesian(ylim=c(0, 100), xlim=c(0,100)) +
#  stat_summary(aes(y = Mature1_pc, x = Concentration),fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red", size=0.1)

# Figure 1B
# 5p titration with highest 3p concentration 25:26,36:38,46:48
# Serum<- PercentDataComplete_NoUMI[c(25:26,36:38,46:48),c(8,1)] %>% mutate(Concentration=c(rep(10,2), rep(30,3), rep(100,3))) # 10, 30, 100 pmoles per 30ul rxn
Serum<- PercentDataComplete_NoUMI[c(25:26,36:38,46:48),c(8,1)] %>% mutate(Concentration=c(rep(0.33,2), rep(1.0 ,3), rep(3.3,3)))
Serum
#
Serumgraph <- ggplot(Serum) + theme_bw(base_size=8) + 
  geom_point(aes(y = Mature1_pc, x = Concentration)) +
  labs(x="5' Adapter Concentration (uM)", y="Percent Mature miRNAs") +
  coord_cartesian(ylim=c(0, 30), xlim=c(0,4)) +
  stat_summary(aes(y = Mature1_pc, x = Concentration),fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red", size=0.1)
# save file
# postscript(file="Figure1B_concentration.eps",paper = "special", width=1.5, height=1.5, horizontal=FALSE)
Serumgraph # 
# dev.off()

# 5p 10nM noMod 3p titration 19:26
# Figure 1C
# 5p 100nM noMod 3p titration 1, 4, 17nM  40:48 
# Serum3p1<- PercentDataComplete_NoUMI[c(40:48),c(8,1)] %>% mutate(Concentration=c(rep(1,3), rep(4,3), rep(17,3)))  # 1, 4, 17 pmoles per 20ul rxn
Serum3p1<- PercentDataComplete_NoUMI[c(40:48),c(8,1)] %>% mutate(Concentration=c(rep(0.05,3), rep(0.2,3), rep(0.85,3)))
Serum3p1
#
Serumgraph2 <- ggplot(Serum3p1) + theme_bw(base_size=8) + 
  geom_point(aes(y = Mature1_pc, x = Concentration)) +
  stat_summary(aes(y = Mature1_pc, x = Concentration),fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red", size=0.1)+
  labs(x="3' Adapter Concentration (uM)", y="Percent Mature miRNAs") +
  coord_cartesian(ylim=c(0, 30), xlim=c(0,1))
# save file
# postscript(file="Figure1C_concentration.eps",paper = "special", width=1.5, height=1.5, horizontal=FALSE)
Serumgraph2 # 
# dev.off()


# Figure 1D

# 7:9, 46:48 3p 17nM mod vs no mod  mutate(Concentration=c(rep("100nM 5p Modified;17nM 3p",3),rep("100nM 5p Not Modified;17nM 3p",3)))  levels=c("100nM 5p Modified;17nM 3p","100nM 5p Not Modified;17nM 3p"))
Serum5pmod<- PercentDataComplete_NoUMI[c(7:9,46:48),c(8,1)] %>%  mutate(Concentration=c(rep("Modified",3),rep("Unmodified",3)))
Serum5pmod
#
Serum5pmod$Concentration <- parse_factor(Serum5pmod$Concentration,  levels=c("Modified","Unmodified"))
Serumgraph2 <- ggplot(Serum5pmod) + theme_bw(base_size=8) + 
  geom_point(aes(y = Mature1_pc, x = Concentration)) +
  stat_summary(aes(y = Mature1_pc, x = Concentration),fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red", size=0.1)+
  labs(x="5' Adapter Modification", y="Percent Mature miRNA") +
  coord_cartesian(ylim=c(0, 30))
# save file
# postscript(file="Figure1D.eps",paper = "special", width=1.5, height=1.5, horizontal=FALSE)
Serumgraph2 # 
# dev.off()
# alternatively change above code to look at lower 3p concentration
# 1:3,40:42 3p 1NM, mod vs no mod  mutate(Concentration=c(rep("100nM 5p Modified;1nM 3p",3),rep("100nM 5p Not Modified;1nM 3p",3)))   levels=c("100nM 5p Modified;1nM 3p","100nM 5p Not Modified;1nM 3p"))

# serum old vs new protocols-------------------------------
# new protocol need to load data from TJ12 library  RB2 samples 
SummarydataTJ12 <- read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ12/TJ12_stats/Map_stats_TJ12_noN_all.txt", header = TRUE, stringsAsFactors = FALSE)
# now break down the UMI contaminant data into markers and five prime adapter --------------
ContamSpecificdata2 <- read.table("/Volumes/BiomarkerSeq/SeqPaperData/TJ12/TJ12_stats/TJ12_Libcontam_counts_raw.txt", header=TRUE, stringsAsFactors = FALSE)
# replace NA with 0
Summarydata2 <- SummarydataTJ12[23:26,]

ContamSpecificdata2[is.na(ContamSpecificdata2)] <- 0
names(ContamSpecificdata2) <- gsub("X4", "4", names(ContamSpecificdata2))
rownames(ContamSpecificdata2) <- ContamSpecificdata2$miRNA

MostContamData2 <- filter(ContamSpecificdata2, miRNA=="19nt_Marker" | miRNA == "24nt_Marker" | miRNA == "5primeAdapter" ) # choose just marker and 5 prime adapter and not phiX data
rownames(MostContamData2) <- MostContamData2$miRNA
MostContamDataT2<-data.frame(t(MostContamData2[,24:27])) # CB2_1f1_1:CB2_1f_4
contamSamples2 <- rownames(MostContamDataT2)
MostContamDataT2 <- MostContamDataT2 %>% mutate(other_contam = Summarydata2$LibraryContam-(X19nt_Marker+X24nt_Marker+X5primeAdapter)) %>% mutate(Sample = contamSamples2)
Samplelist2 <- Summarydata2$Sample
# combine the count data
completeSummary_NoUMI2 <- right_join(Summarydata2, MostContamDataT2, by = "Sample") %>% select(-DemultiCount, -LibraryContam, -Hairpins1) %>% 
  mutate(Markers=X19nt_Marker+X24nt_Marker, otherSmallRNA=Ribo+tRNA+svRNA) %>%
  select(-X19nt_Marker,-X24nt_Marker,-Ribo,-tRNA,-svRNA) %>% rename(Adapter=X5primeAdapter) %>%
  mutate(Unmapped = (TrimCount-(Markers+Adapter+other_contam+Mature1+otherSmallRNA+Genomic))) 

PercentDataComplete_NoUMI2 <- completeSummary_NoUMI2 %>% select(-(Sample)) %>% mutate_all(funs(pc=(./TrimCount*100))) %>% select(ends_with("_pc"), -TrimCount_pc) %>%  mutate(Sample = Samplelist2)

# PercentDatalong_NoUMI2 <- gather(PercentDataComplete_NoUMI2, MapCategory, PercentTotal, Markers_pc, Adapter_pc, other_contam_pc, Hairpins1_pc, otherSmallRNA_pc, Genomic_pc, Unmapped_pc)
NewData <- PercentDataComplete_NoUMI2[ ,c(8,1)] # just mature

OldData <- PercentDataComplete_NoUMI[c(46:48),c(8,1)] # just mature of 100nM no Mod, 17nM 3prime
# Concentration=c(rep("100nM 5p Not Modified;17nM 3p",3), rep("10nM 5p Modified;1nM 3p",4))
OldvNew<- rbind(OldData, NewData) %>% mutate(Concentration=c(rep(100.17,3), rep(10.01,4)))
OldvNew
# OldvNewMeans <-summarise(group_by(OldvNew, Concentration), m=mean(Hairpins1_pc), sd=sd(Hairpins1_pc))
# Serum$MapCategory <- parse_factor(Serum$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Hairpins1_pc","other_contam_pc","Adapter_pc","Markers_pc"))
# OldvNew$Concentration <- parse_factor(OldvNew$Concentration, levels=c("10nM 5p Modified;1nM 3p","100nM 5p Not Modified;17nM 3p"))
# Figure 1E
#  label x as old and optimized
OldvNew<- rbind(OldData, NewData) %>% mutate(Concentration=c(rep("Original",3), rep("Optimized",4)))
# OldvNewMeans <-summarise(group_by(OldvNew, Concentration), m=mean(Hairpins1_pc), sd=sd(Hairpins1_pc))
# Serum$MapCategory <- parse_factor(Serum$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Hairpins1_pc","other_contam_pc","Adapter_pc","Markers_pc"))
 OldvNew$Concentration <- parse_factor(OldvNew$Concentration, levels=c("Optimized","Original"))
OldvNew
 # Figure 1E
Serumgraph3b <- ggplot(OldvNew) + theme_bw(base_size=8) +
  geom_point(aes(y = Mature1_pc, x = Concentration)) +
  labs(x="Adapter Conditions", y="Percent Mature miRNA") +
  coord_cartesian(ylim=c(0, 60)) +
  stat_summary(aes(y = Mature1_pc, x = Concentration),fun.y = mean, fun.ymin = min, fun.ymax = max, colour = "red", size=0.1)

# save file
# postscript(file="Figure1E.eps",paper = "special", width=1.5, height=1.5, horizontal=FALSE)
Serumgraph3b # 
# dev.off()  

# Old vs Optimized Stacked Bar

OldvNewBar <- rbind(PercentDataComplete_NoUMI[c(46:48),],PercentDataComplete_NoUMI2) %>% 
    mutate(Concentration=(c("Original-1","Original-2","Original-3","Optimized-1","Optimized-2","Optimized-3","Optimized-4"))) %>% select(-Sample)

OldvNewBarlong <- gather(OldvNewBar, MapCategory, PercentTotal, Markers_pc, Adapter_pc, other_contam_pc, Mature1_pc, otherSmallRNA_pc, Genomic_pc, Unmapped_pc) 
OldvNewBarlong$MapCategory <- parse_factor(OldvNewBarlong$MapCategory, levels=c("Unmapped_pc","Genomic_pc","otherSmallRNA_pc","Mature1_pc","other_contam_pc","Adapter_pc","Markers_pc"))

OldvNewBgraph <- ggplot(OldvNewBarlong) + theme_bw(base_size=8) +
   geom_bar(aes(y = PercentTotal, x = Concentration, fill = MapCategory), stat="identity") + 
  theme(axis.text.x  = element_text(angle=90, vjust=.5, size=8)) +
  labs(x="Adapter Conditions", y="Percent Mapped") +
  scale_fill_manual(values=mypalette)

# postscript(file="Figure1F.eps",paper = "special", width=6, height=3, horizontal=FALSE)
OldvNewBgraph # 
# dev.off()

OldvNewBar
# system information at time of this rendering:
sessionInfo()
