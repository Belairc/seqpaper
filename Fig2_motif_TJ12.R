# Note to save the Plots remove the # from the postscript() and dev.off() lines

library('tidyverse')
theme_set(theme_minimal())

Seqdir="~/Box Sync/SeqPaperData/TJ12/TJ12_collapse12" # remote directory where data is

outdir<-"~/Box Sync/SeqPaperData/TJ12/motif/cell4_mature"
# the SampleSummaryFile is a temporary file needed.
project="TJ12"
#
setwd(Seqdir)

# > head(dir(),14)
# [1] "4N01_ccell4_500ng_1_genomic_count12_seqs.txt"   "4N01_ccell4_500ng_1_genomic_count12.txt"       
# [3] "4N01_ccell4_500ng_1_Hairpin_count12_seqs.txt"   "4N01_ccell4_500ng_1_Hairpin_count12.txt"       
# [5] "4N01_ccell4_500ng_1_Libcontam_count12_seqs.txt" "4N01_ccell4_500ng_1_Libcontam_count12.txt"     
# [7] "4N01_ccell4_500ng_1_Mature_count12_seqs.txt"    "4N01_ccell4_500ng_1_Mature_count12.txt"        
# [9] "4N01_ccell4_500ng_1_ribo_count12_seqs.txt"      "4N01_ccell4_500ng_1_ribo_count12.txt"          
# [11] "4N01_ccell4_500ng_1_svrna_count12_seqs.txt"     "4N01_ccell4_500ng_1_svrna_count12.txt"         
# [13] "4N01_ccell4_500ng_1_trna_count12_seqs.txt"      "4N01_ccell4_500ng_1_trna_count12.txt"          
# > 


  # filelist=dir()

 library(stringr)
# library(dplyr)
# filename=filelist[seq(1,length(filelist), 50)] # if you don't want to analyze all, choose subset here
# use this for looping 
#
 filename="4N01_ccell4_500ng_1_Mature_count12_seqs.txt" #
infoFile=read.table(filename, header=FALSE, stringsAsFactors = FALSE)
# $1=HairpinName, $2=leftSeq/insert, $3=AdapterMatch, $4=WholeSeq, $5=4N-miRNA-4Nbarcode4N
names(infoFile) <- c("miRNA","Insert","AdapterMatch","WholeSeq","TrimmedSeq")
# library(seqLogo)
# library(Biostrings)
# 
# library("devtools")
# install_github("omarwagih/ggseqlogo")

require(ggplot2)
require(ggseqlogo)

Sample <- infoFile %>% mutate(FiveP4nt = (str_sub(infoFile$Insert, 1, 4)), ThreePUMI=(str_sub(AdapterMatch, 1, 13)))
setwd(outdir)
Sample_5pAll <- Sample %>% select(FiveP4nt)
seqsAll_5p <- as.character(Sample_5pAll[!grepl("N",Sample_5pAll$FiveP4nt),])

#  save file

 postscript(file="Figure2_5p.eps",paper = "special", width=1.25, height=1.5, horizontal=FALSE)
ggplot() + geom_logo(seqsAll_5p, method='p') + theme_logo(base_size = 8, base_family = "") # 
 dev.off()


Sample_3pAll <- Sample %>% select(ThreePUMI)
seqsAll_3p <- as.character(Sample_3pAll[!grepl("N",Sample_3pAll$ThreePUMI),])

 postscript(file="Figure2_3p.eps",paper = "special", width=3, height=1.5, horizontal = FALSE)
# par(mar=c(5,3,2,2)+0.1)
# seqLogo(consensusMatrix(seqsAll_3p, as.prob = T)[1:4,], ic.scale=F)# 
ggplot() + geom_logo(seqsAll_3p, method='p') + theme_logo(base_size = 8, base_family = "")
 dev.off()

# 
# now find motifs for top 25 expressed mature miRNAs
# open counts file
countsdir<-"~/Box Sync/SeqPaperData/TJ12/TJ12_collapse12"
setwd(countsdir)
# Countsfilelist=dir()
filename2="4N01_ccell4_500ng_1_Mature_count12.txt" #
# filename2=Countsfilelist[seq(1,length(filelist), 50)]
countsfile <- read.table(filename2, header=FALSE, stringsAsFactors = FALSE)
countsfile <- arrange(countsfile, desc(V2)) 
top25mirs <- countsfile[1:25,1]
#
# > top25mirs
# [1] "hsa-miR-21-5p"   "hsa-miR-24-3p"   "hsa-miR-29a-3p"  "hsa-let-7a-5p"   "hsa-miR-27a-3p" 
# [6] "hsa-miR-100-5p"  "hsa-miR-29b-3p"  "hsa-miR-125b-5p" "hsa-miR-16-5p"   "hsa-miR-191-5p" 
# [11] "hsa-let-7f-5p"   "hsa-miR-22-3p"   "hsa-miR-103b"    "hsa-miR-23a-3p"  "hsa-miR-96-5p"  
# [16] "hsa-let-7b-5p"   "hsa-miR-26a-5p"  "hsa-miR-221-3p"  "hsa-let-7i-5p"   "hsa-miR-320a"   
# [21] "hsa-let-7g-5p"   "hsa-miR-3074-5p" "hsa-miR-15a-5p"  "hsa-miR-103a-3p" "hsa-miR-424-5p" setwd(outdir)
setwd(outdir)
for (mir in top25mirs){
  Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
  seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
  #
  motif_file <- paste(mir,"_5p.eps", sep="")
  #  save file
  postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
  print(ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = ""))
  dev.off()
  
  Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
  seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
  motif_file2 <- paste(mir,"_3p.eps", sep="")
  
   postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
  print(ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = ""))
   dev.off()
}
# now look for the highest one in serum, which is low in the cells
mir = "hsa-miR-451a"
Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
motif_file <- paste(mir,"_5p_2.eps", sep="")
pt= ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = "")

#  save file
 postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
print(pt)
 dev.off()


Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
motif_file2 <- paste(mir,"_3p_2.eps", sep="")
# save file
pt2=ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = "")
postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
print(pt2)
dev.off()


# ccell 4 #2
setwd(Seqdir)
outdir<-"~/Box Sync/SeqPaperData/TJ12/motif/cell4_mature2/"
filename="4N02_ccell4_500ng_2_Mature_count12_seqs.txt" 

infoFile=read.table(filename, header=FALSE, stringsAsFactors = FALSE)
# $1=HairpinName, $2=leftSeq/insert, $3=AdapterMatch, $4=WholeSeq, $5=4N-miRNA-4Nbarcode4N
names(infoFile) <- c("miRNA","Insert","AdapterMatch","WholeSeq","TrimmedSeq")
Sample <- infoFile %>% mutate(FiveP4nt = (str_sub(infoFile$Insert, 1, 4)), ThreePUMI=(str_sub(AdapterMatch, 1, 13)))
# 
Sample_5pAll <- Sample %>% select(FiveP4nt)
seqsAll_5p <- as.character(Sample_5pAll[!grepl("N",Sample_5pAll$FiveP4nt),])

setwd(outdir)
#  save file
postscript(file="Figure2_5p_2.eps",paper = "special", width=1.25, height=1.5, horizontal=FALSE)
print(ggplot() + geom_logo(seqsAll_5p, method='p') + theme_logo(base_size = 8, base_family = "")) # 
dev.off()


Sample_3pAll <- Sample %>% select(ThreePUMI)
seqsAll_3p <- as.character(Sample_3pAll[!grepl("N",Sample_3pAll$ThreePUMI),])
postscript(file="Figure2_3p_2.eps",paper = "special", width=3, height=1.5, horizontal = FALSE)

print(ggplot() + geom_logo(seqsAll_3p, method='p') + theme_logo(base_size = 8, base_family = ""))
dev.off()

# 
# 
# now find motifs for top 25 expressed mature miRNAs
# open counts file
countsdir<-"~/Box Sync/SeqPaperData/TJ12/TJ12_collapse12"
setwd(countsdir)
# Countsfilelist=dir()
filename2="4N02_ccell4_500ng_2_Mature_count12.txt" #
# filename2=Countsfilelist[seq(1,length(filelist), 50)]
countsfile <- read.table(filename2, header=FALSE, stringsAsFactors = FALSE)
countsfile <- arrange(countsfile, desc(V2)) 
top25mirs <- countsfile[1:25,1]
#
# > top25mirs
# [1] "hsa-miR-21-5p"   "hsa-miR-24-3p"   "hsa-let-7a-5p"   "hsa-miR-27a-3p"  "hsa-miR-100-5p" 
# [6] "hsa-miR-191-5p"  "hsa-miR-29a-3p"  "hsa-miR-16-5p"   "hsa-miR-29b-3p"  "hsa-miR-125b-5p"
# [11] "hsa-let-7i-5p"   "hsa-miR-103b"    "hsa-miR-26a-5p"  "hsa-let-7f-5p"   "hsa-miR-22-3p"  
# [16] "hsa-miR-3074-5p" "hsa-miR-96-5p"   "hsa-miR-221-3p"  "hsa-let-7g-5p"   "hsa-miR-23a-3p" 
# [21] "hsa-miR-93-5p"   "hsa-miR-320a"    "hsa-let-7b-5p"   "hsa-miR-103a-3p" "hsa-miR-17-5p"  

setwd(outdir)
for (mir in top25mirs) {
  Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
  seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
  motif_file <- paste(mir,"_5p_2.eps", sep="")
  pt= ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = "")
  #  save file
  postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
  print(pt)
  dev.off()

  
  Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
  seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
  motif_file2 <- paste(mir,"_3p_2.eps", sep="")
  # save file
  pt2=ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = "")
  postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
  print(pt2)
  dev.off()

}
# now look for the highest one in serum, which is low in the cells
mir = "hsa-miR-451a"
Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
motif_file <- paste(mir,"_5p_2.eps", sep="")
pt= ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = "")
#  save file
postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
print(pt)
dev.off()


Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
motif_file2 <- paste(mir,"_3p_2.eps", sep="")
# save file
pt2=ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = "")
postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
print(pt2)
dev.off()

# serum CB2_1f1
setwd(Seqdir)
outdir<-"~/Box Sync/SeqPaperData/TJ12/motif/CB2_1f1_mature/"
filename="4N23_CB2_1f_1_Mature_count12_seqs.txt" 

infoFile=read.table(filename, header=FALSE, stringsAsFactors = FALSE)
# $1=HairpinName, $2=leftSeq/insert, $3=AdapterMatch, $4=WholeSeq, $5=4N-miRNA-4Nbarcode4N
names(infoFile) <- c("miRNA","Insert","AdapterMatch","WholeSeq","TrimmedSeq")

# 
Sample <- infoFile %>% mutate(FiveP4nt = (str_sub(infoFile$Insert, 1, 4)), ThreePUMI=(str_sub(AdapterMatch, 1, 13)))
Sample_5pAll <- Sample %>% select(FiveP4nt)
seqsAll_5p <- as.character(Sample_5pAll[!grepl("N",Sample_5pAll$FiveP4nt),])

setwd(outdir)
#  save file
# postscript(file="Figure2_5p_1b.eps",paper = "special", width=1.25, height=1.5, horizontal=FALSE)
print(ggplot() + geom_logo(seqsAll_5p, method='p') + theme_logo(base_size = 8, base_family = "")) # 
# dev.off()


Sample_3pAll <- Sample %>% select(ThreePUMI)
seqsAll_3p <- as.character(Sample_3pAll[!grepl("N",Sample_3pAll$ThreePUMI),])

#  save file
# postscript(file="Figure2_3p_1b.eps",paper = "special", width=3, height=1.5, horizontal = FALSE)
print(ggplot() + geom_logo(seqsAll_3p, method='p') + theme_logo(base_size = 8, base_family = ""))
# dev.off()

# 
# 
# now find motifs for top 25 expressed mature miRNAs
# open counts file
countsdir<-"~/Box Sync/SeqPaperData/TJ12/TJ12_collapse12"
setwd(countsdir)
# Countsfilelist=dir()
filename2="4N23_CB2_1f_1_Mature_count12.txt" #
# filename2=Countsfilelist[seq(1,length(filelist), 50)]
countsfile <- read.table(filename2, header=FALSE, stringsAsFactors = FALSE)
countsfile <- arrange(countsfile, desc(V2)) 
top25mirs <- countsfile[1:25,1]
#
# > top25mirs
# [1] "hsa-miR-21-5p"   "hsa-miR-24-3p"   "hsa-let-7a-5p"   "hsa-miR-27a-3p"  "hsa-miR-100-5p" 
# [6] "hsa-miR-191-5p"  "hsa-miR-29a-3p"  "hsa-miR-16-5p"   "hsa-miR-29b-3p"  "hsa-miR-125b-5p"
# [11] "hsa-let-7i-5p"   "hsa-miR-103b"    "hsa-miR-26a-5p"  "hsa-let-7f-5p"   "hsa-miR-22-3p"  
# [16] "hsa-miR-3074-5p" "hsa-miR-96-5p"   "hsa-miR-221-3p"  "hsa-let-7g-5p"   "hsa-miR-23a-3p" 
# [21] "hsa-miR-93-5p"   "hsa-miR-320a"    "hsa-let-7b-5p"   "hsa-miR-103a-3p" "hsa-miR-17-5p"  

setwd(outdir)
for (mir in top25mirs) {
  Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
  seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
  motif_file <- paste(mir,"_5p_1b.eps", sep="")
  pt= ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = "")
  #  save file
  # postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
  print(pt)
  # dev.off()
  
  
  Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
  seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
  motif_file2 <- paste(mir,"_3p_1b.eps", sep="")
  pt2=ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = "")
  # save file
  # postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
  print(pt2)
  # dev.off()
  
}
# now look for the highest one in serum, which is low in the cells
mir = "hsa-miR-96-5p"
Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
motif_file <- paste(mir,"_5p_2.eps", sep="")
pt= ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = "")
#  save file
# postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
print(pt)
# dev.off()


Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
motif_file2 <- paste(mir,"_3p_2.eps", sep="")
# save file
pt2=ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = "")

# postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
print(pt2)
# dev.off()




# serum CB2_1f2
setwd(Seqdir)
outdir<-"~/Box Sync/SeqPaperData/TJ12/motif/CB2_1f2_mature/"
filename="4N24_CB2_1f_2_Mature_count12_seqs.txt" 

infoFile=read.table(filename, header=FALSE, stringsAsFactors = FALSE)
# $1=HairpinName, $2=leftSeq/insert, $3=AdapterMatch, $4=WholeSeq, $5=4N-miRNA-4Nbarcode4N
names(infoFile) <- c("miRNA","Insert","AdapterMatch","WholeSeq","TrimmedSeq")

# 
Sample <- infoFile %>% mutate(FiveP4nt = (str_sub(infoFile$Insert, 1, 4)), ThreePUMI=(str_sub(AdapterMatch, 1, 13)))
Sample_5pAll <- Sample %>% select(FiveP4nt)
seqsAll_5p <- as.character(Sample_5pAll[!grepl("N",Sample_5pAll$FiveP4nt),])

setwd(outdir)
#  save file
# postscript(file="Figure2_5p_2b.eps",paper = "special", width=1.25, height=1.5, horizontal=FALSE)
print(ggplot() + geom_logo(seqsAll_5p, method='p') + theme_logo(base_size = 8, base_family = "")) # 
# dev.off()


Sample_3pAll <- Sample %>% select(ThreePUMI)
seqsAll_3p <- as.character(Sample_3pAll[!grepl("N",Sample_3pAll$ThreePUMI),])

# save file
# postscript(file="Figure2_3p_2b.eps",paper = "special", width=3, height=1.5, horizontal = FALSE)
print(ggplot() + geom_logo(seqsAll_3p, method='p') + theme_logo(base_size = 8, base_family = ""))
# dev.off()

# 
# 
# now find motifs for top 25 expressed mature miRNAs
# open counts file
countsdir<-"~/Box Sync/SeqPaperData/TJ12/TJ12_collapse12"
setwd(countsdir)
# Countsfilelist=dir()
filename2="4N24_CB2_1f_2_Mature_count12.txt" #
# filename2=Countsfilelist[seq(1,length(filelist), 50)]
countsfile <- read.table(filename2, header=FALSE, stringsAsFactors = FALSE)
countsfile <- arrange(countsfile, desc(V2)) 
top25mirs <- countsfile[1:25,1]
#
# > top25mirs
# [1] "hsa-miR-21-5p"   "hsa-miR-24-3p"   "hsa-let-7a-5p"   "hsa-miR-27a-3p"  "hsa-miR-100-5p" 
# [6] "hsa-miR-191-5p"  "hsa-miR-29a-3p"  "hsa-miR-16-5p"   "hsa-miR-29b-3p"  "hsa-miR-125b-5p"
# [11] "hsa-let-7i-5p"   "hsa-miR-103b"    "hsa-miR-26a-5p"  "hsa-let-7f-5p"   "hsa-miR-22-3p"  
# [16] "hsa-miR-3074-5p" "hsa-miR-96-5p"   "hsa-miR-221-3p"  "hsa-let-7g-5p"   "hsa-miR-23a-3p" 
# [21] "hsa-miR-93-5p"   "hsa-miR-320a"    "hsa-let-7b-5p"   "hsa-miR-103a-3p" "hsa-miR-17-5p"  

setwd(outdir)
for (mir in top25mirs) {
  Sample_mir5p <- Sample %>% filter(miRNA==mir) %>% select(FiveP4nt)
  seqs_mir_5p <- as.character(Sample_mir5p[!grepl("N",Sample_mir5p$FiveP4nt),])
  motif_file <- paste(mir,"_5p_2b.eps", sep="")
  pt= ggplot() + geom_logo( seqs_mir_5p, method='p') + theme_logo(base_size = 8, base_family = "")
  #  save file
  # postscript(file=motif_file,paper = "special", width=1.5, height=1.5, horizontal = FALSE)
  print(pt)
  # dev.off()
  
  
  Sample_mir3p <- Sample %>% filter(miRNA==mir) %>% select(ThreePUMI)
  seqs_mir_3p <- as.character(Sample_mir3p[!grepl("N",Sample_mir3p$ThreePUMI),])
  motif_file2 <- paste(mir,"_3p_2b.eps", sep="")
  # save file
  pt2=ggplot() + geom_logo( seqs_mir_3p, method='p') + theme_logo(base_size = 8, base_family = "")
  
  # postscript(file=motif_file2, paper = "special", width=3, height=1.5, horizontal = FALSE)
  print(pt2)
  # dev.off()
  
}


# system information at time of this rendering:
sessionInfo()