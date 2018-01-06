# this script will take the uncollapsed mapped sequence data and use R to collapse 
# based on varying length of UMI
# The data will be the TJ12 library of cells and serum with 12 cycles
#


#

library(stringr)
library('tidyverse')
theme_set(theme_classic(base_size=8))

# -----repeat for Mature miRNAs -------

# basedir="~/Box Sync/SeqPaperData/TJ12" #  seq data here
basedir="/Volumes/BiomarkerSeq/SeqPaperData/TJ12"  # seq data on box moved to cloud and local version stored on attached hard drive
setwd(basedir)
outdir<-"~/Box Sync/SeqPaperDataCB/RAnalysis" # 
project="TJ12"
#

filename="4N01_ccell4_500ng_1_raw_seqs.txt" #
infoFile=read.table(filename, header=FALSE, stringsAsFactors = FALSE)
# $1=HairpinName, $2=leftSeq/insert, $3=AdapterMatch, $4=WholeSeq, $5=4N-miRNA-4Nbarcode4N
names(infoFile) <- c("miRNA","TrimmedSeq")

noCollapse <- infoFile %>% group_by(miRNA) %>% summarise(n_distinct(miRNA), uncollapse=n()) %>% select(miRNA, uncollapse)

# collapse from the 5p end 1 at a time  (I checked that these calculations are the same as the shell collapse on the same sequence)
collapsed1_01 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 1,-1)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert))%>%
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x1_m01=n()) %>% arrange(desc(x1_m01)) %>% select(miRNA, x1_m01)
# #1
collapsed1_02 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 1,-2)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x1_m02 = n()) %>% arrange(desc(x1_m02)) %>% select(miRNA, x1_m02)
# #2
collapsed1_03 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 1,-3)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x1_m03 = n()) %>% arrange(desc(x1_m03)) %>% select(miRNA, x1_m03)
# #3
collapsed1_04 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 1,-4)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x1_m04 = n()) %>% arrange(desc(x1_m04)) %>% select(miRNA, x1_m04)
# #4
# collapsed1_05 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 1,-5)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
#  group_by(miRNA) %>% summarise(n_distinct(miRNA), x1_m05 = n()) %>% arrange(desc(x1_m05)) %>% select(miRNA, x1_m05)
# #5
collapsed1_10 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 1,-10)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x1_m10 = n()) %>% arrange(desc(x1_m10)) %>% select(miRNA, x1_m10)
# #6
collapsed2_10 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 2,-10)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x2_m10 = n()) %>% arrange(desc(x2_m10)) %>% select(miRNA, x2_m10)
# #7
collapsed2_11 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 2,-11)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x2_m11 = n()) %>% arrange(desc(x2_m11)) %>% select(miRNA, x2_m11)
# #8
collapsed3_11 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 3,-11)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x3_m11 = n()) %>% arrange(desc(x3_m11)) %>% select(miRNA, x3_m11)
# #9
collapsed3_12 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 3,-12)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x3_m12 = n()) %>% arrange(desc(x3_m12)) %>% select(miRNA, x3_m12)
# #10
collapsed4_12 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 4,-12)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x4_m12 = n()) %>% arrange(desc(x4_m12)) %>% select(miRNA, x4_m12)
# #11
collapsed4_13 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 4,-13)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x4_m13 = n()) %>% arrange(desc(x4_m13)) %>% select(miRNA, x4_m13)
#
collapsed5_14 <- infoFile %>% mutate(Insert=str_sub(TrimmedSeq, 5,-14)) %>% group_by(miRNA, Insert) %>% summarise(n_distinct(Insert)) %>% 
  group_by(miRNA) %>% summarise(n_distinct(miRNA), x5_m14 = n())  %>% arrange(desc(x5_m14)) %>% select(miRNA, x5_m14)
#
AllCollapsed <- list(noCollapse,collapsed1_01,collapsed1_02,collapsed1_03,collapsed1_04,collapsed1_10,collapsed2_10,collapsed2_11,collapsed3_11,collapsed3_12,collapsed4_12,collapsed4_13,collapsed5_14) %>%
  Reduce( function(dtf1,dtf2) full_join(dtf1,dtf2, by = "miRNA"), .)
AllCollapsed_high <- filter(AllCollapsed, x1_m01>20)
library(tidyr)
AllCollapsed_highLong <- gather(AllCollapsed_high, key=Ns, value=counts, uncollapse,x1_m01,x1_m02,x1_m03,x1_m04,x1_m10,x2_m10,x2_m11,x3_m11,x3_m12,x4_m12,x4_m13,x5_m14)  

setwd(outdir)
library(ggplot2)
# library(RColorBrewer)
# mypalette<-brewer.pal(11,"Paired")
# Cell4 Mature Counts after removing indicated N from UMI\n
bargraph1 <- ggplot(AllCollapsed_highLong) +  ggtitle("All miRNAs")
# graph_file <- paste(mir,".png", sep="")
# postscript(file="Figure3A.eps",paper = "special", width=4, height=2, horizontal=FALSE)
bargraph1 + geom_bar(aes(y = counts, x = Ns), stat="identity") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) +
  theme(legend.position="none")
# dev.off()
#
# examine collapse for select mirs 
# open counts file
selectmirs <- c("hsa-miR-21-5p","hsa-let-7i-5p","hsa-miR-96-5p","hsa-miR-320a","hsa-miR-126-3p")
for(mir in selectmirs){
  print(mir)
  hsa <- filter(AllCollapsed_highLong, miRNA==mir)
  # original ggtitle(paste(mir, " Counts after removing indicated N from UMI\n Cell4-1",sep=""))
  graph<- ggplot(hsa) +  ggtitle(paste(mir))+ 
    geom_bar(aes(y = counts, x = Ns), stat="identity") + 
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) +
    theme(legend.position="none")
  #
  graph_file <- paste("Figure3_",mir,"cell4_1.eps", sep="")
#   postscript(file=graph_file, paper = "special", width=4, height=2, horizontal=FALSE)
  print(graph)
 #  dev.off()
}
# current version of figure used miR-21 and miR-320

# system information at time of this rendering:
sessionInfo()
