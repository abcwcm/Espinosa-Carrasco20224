# conda activate diffbind_v3 ; /athena/abc/scratch/paz2005/miniconda3/envs/diffbind_v3/bin/R

library(DiffBind) # ‘3.0.13’
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(magrittr)
library(openxlsx)
library(data.table)



DB <- dba(sampleSheet = "atac_decoder_diffBind.csv", peakCaller = "macs", peakFormat = "narrow", config=data.frame(AnalysisMethod=DBA_DESEQ2))
DB <- dba.blacklist(DB, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE) 
DB_consensus<-dba.peakset(DB,consensus=-DBA_REPLICATE, minOverlap=0.5) 
consensus_peaks <- dba.peakset(DB_consensus, bRetrieve=TRUE)
DB <- dba.count(DB,  bUseSummarizeOverlaps=T, bParallel=TRUE, peaks=consensus_peaks)
DB <- dba.contrast(DB, categories=DBA_CONDITION, minMembers=2)
saveRDS(DB, "2022_01.DB_gabriel.Rds")


# export counts
DB <- readRDS("2022_01.DB_gabriel.Rds")
db <- dba.count(DB, peaks=NULL, score=DBA_SCORE_READS)
consensus_peaks <- dba.peakset(db, bRetrieve=TRUE)
counts = as.data.frame(mcols(consensus_peaks))
row.names(counts) = paste0(as.data.frame(consensus_peaks)$seqnames, ":",as.data.frame(consensus_peaks)$start,"-",as.data.frame(consensus_peaks)$end)
decoderFile <- "atac_decoder.txt"
decoder.data <- fread(decoderFile) %>% as.data.frame()
decoder.data$group <- factor(decoder.data$group, levels=c("TCR.OTI", "TCR.OTI.SMARTA"))
decoder.data <- decoder.data[decoder.data$sample.ID %in% colnames(counts),]
counts <- counts[,c(decoder.data$sample.ID)]
if(!identical(decoder.data$sample.ID, colnames(counts))) stop()
write.csv(counts,"atac_countsPerPeakAtlas.csv")

