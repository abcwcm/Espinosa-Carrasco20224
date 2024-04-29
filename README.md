![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)

# Bioinformatic methods for Espinosa-Carrasco et al.
The quality of the sequenced reads was assessed with FastQC and QoRTs. Unless stated otherwise, plots involving high- throughput sequencing data were created using R version 4.1.0 and ggplot2.

## RNA-seq data:
## Alignment and gene counting
DNA sequencing reads were aligned with default parameters to the mouse reference genome (GRCm38.p6) using STAR v2.6.0c. Gene expression estimates were obtained with featureCounts v1.6.2 using composite gene models (union of the exons of all transcript isoforms per gene) from Gencode (version M17).

### DEGs
DEGs were determined using DESeq2 v1.34.0 with Wald tests with a q-value cutoff of 0.05 (Benjamini–Hochberg correction).

### Heatmaps
Heatmaps were created using DESeq2 normalized read counts after variance stabilizing transformation of genes identified as differentially expressed by DESeq2. Rows were centered and scaled.

### Pathway and GO term enrichment analyses
Gene set enrichment analyses were done using fgsea v1.20.0 with the fgseaMultilevel function. Genes were ranked based on the DESeq2 Wald statistic. Gene sets with an FDR < 0.05 were considered enriched. Gene ontology analysis was performed on up- and down-regulated DEGs using the clusterProfiler v4.2.2 R package. Only GO categories enriched using a 0.05 false discovery rate cutoff were considered.

## ATAC-seq data:
### Alignment and creation of peak atlas
Reads were aligned to the mouse reference genome (version GRCm38) with BWA-backtrack v0.7.17 (Li and Durbin, 2009). Post-alignment filtering was done with samtools v1.8 and Picard tools v2.18.9 (Li et al., 2009) to remove unmapped reads, improperly paired reads, nonunique reads, and duplicates. Peaks were called with MACS2 v2.1.1 (Liu, 2014),and peaks with adjusted P values smaller than 0.01 were excluded. Consensus peak sets were generated for each condition if a peak was found in at least two replicates. Reproducible peaks from each condition were merged with DiffBind v3.4.11 to create an atlas of accessible peaks, which was used for downstream analyses. The peak atlas was annotated using the ChIPseeker v1.30.3 and TxDb.Mmusculus.UCSC.mm10.knownGene. Blacklisted regions were excluded (https://sites.google.com/site/anshulkundaje/projects/blacklists).

### Differentially accessible regions
Regions where the chromatin accessibility changed between different conditions were identified with DESeq2 v1.34.0, and only Benjamini–Hochberg corrected P values < 0.05 were considered statistically significant.

### Coverage files
Genome coverage files were normalized for differences in sequencing depth (RPGC normalization) with bamCoverage from deepTools v3.1.0. Replicates were averaged together using UCSC-tools bigWigMerge. Merged coverage files were used for display in Integrated Genomics Viewer shown in Fig. 2e.

### Heatmaps
Heatmaps based on the differentially accessible peaks identified between TCROT1(-CD4) and TCROT1(+CD4) as shown in Fig. 2d were created using profileplyr v1.10.2 and ComplexHeatmap v2.15.1, by binning the region +/− 1kb around the peak summits in 20bp bins. To improve visibility, bins with read counts greater than the 75th percentile + 1.5*IQR were capped at that value.

### Motif analyses
For identifying motifs enriched in differentially accessible peaks, we utilized HOMER via marge v0.0.4. HOMER was run separately on hyper- or hypo-accessible peaks with the flags -size given and -mask. Motifs enriched in hyper- or hypo-accessible peaks were determined by comparing the rank differences (based on P value). The consensus peakset identified by DiffBind was used as the background set.

## Package versions
  
## References
