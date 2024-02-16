
The following data compiled from experts of NFDI4Microbiota members.

The goal is to gather this information centrally in NFDI4Microbiota to assess the span of tools and derive consortium-wide guidelines for the tightest integration possible.



# RNA-seq


### 🚩 Library contamination or highly homogenous library 

<ul> A high number of short, low quality or completely identical reads pollute the library. 

 Possible reasons are: contamination e.g. adapter dimers (adapter+adapter, no DNA). 

 💡 increase starting material / trim adapter /exclude by size  </ul> 


### 🚩 Batch-Wise Replicate Combination 

<ul> combine replicates from multiple batches without any replicates within each batch. 

 Possible reasons are: it is not possible to distinguish between the biological and the technical variance -> bad design. 

 💡 There is no general solution.  </ul> 


### 🚩 Wrong enrichment 

<ul> enrichment of wrong population (e.g. polyA-enrichment to compare to rRNAs). 

 Possible reasons are: There is no general reason. 

 💡 exclude rRNAs  </ul> 


### 🚩 Low 3'-End Sequencing Quality (QC)

<ul> After trimming very few reads remain or are very short and cannot be mapped unambigously.. 

 Possible reasons are: general effect esp. for older data, could be hint for 3'-adapters. 

 💡 check trim again  </ul> 


### 🚩 Low Overall Base Call Quality (QC)

<ul> Overall (>50% position) low (red area, PHRED<20) base call quality. 

 Possible reasons are: universally very good reads + wrong phred-type detection. 

 💡 check encoding detected by fastqc, could be a FP  </ul> 


### 🚩 Flow cell problems (QC)

<ul> red spots in per tile sequence quality. 

 Possible reasons are: damaged flow cell or debris. 

 💡 exclude reads from these cells, buy new flow cell  </ul> 


### 🚩 Problem with Per Base Sequence Content (QC)

<ul> high per base sequence content at the 5' end high per base sequence content. 

 Possible reasons are: adapter sequence still present. 

 💡 trim adapters  </ul> 

<ul> high per base sequence content at the 5' end high per base sequence content. 

 Possible reasons are: TSS reads -> start of read = start of gene = bias distr.. 

 💡 -  </ul> 


### 🚩 Probelm with Per Base Sequence Content (QC)

<ul> C is missing in high per base sequence content. 

 Possible reasons are: sodium bisulphite treated = C->T . 

 </ul> 

<ul> overall high high per base sequence content. 

 Possible reasons are: contamination e.g. adapter dimers. 

 </ul> 


### 🚩 Problem with GC Content (QC)

<ul> 2/more peaks in per sequence GC content. 

 Possible reasons are: multiple unrelated species present. 

 💡 remove contaminated species  </ul> 


### 🚩 Problem with Overrepresented Sequences (QC)

<ul> 2/more peaks in per sequence GC content. 

 Possible reasons are: other contamination, e.g. adapter dimers. 

 💡 adapter dimers => trim  </ul> 

<ul> high percentage of overrepresented sequences. 

 Possible reasons are: contamination present in reads. 

 💡 blast most abudant reads, can be adapter sequences too  </ul> 


### 🚩 Problem with Duplicated Sequences (QC)

<ul> high percentage of overrepresented sequences. 

 Possible reasons are: overrepresented sequence present in data e.g. rRNAs. 

 </ul> 

<ul> high percentage of duplicated sequences. 

 Possible reasons are: low complexity library. 

 💡 little starting material / heavy PCR  </ul> 

<ul> high percentage of duplicated sequences. 

 Possible reasons are: constrained library (only reads starting at TSS). 

 💡 You can use random barcoding to distigush between biol. and tech. replicates if needed  </ul> 


### 🚩 Trimmed Reads with Adapters (Trimming)

<ul> no reads are trimmed although adapters are present. 

 Possible reasons are: hard coded adapter sequences, .... 

 💡 adjust parameters  </ul> 


### 🚩 Wrongly Trimmed Reads 

<ul> trimmed reads although no adapters are present. 

 Possible reasons are: false positive result of TG. 

 💡 skip trim-galore  </ul> 


### 🚩 Extremely low read counts (Mapping)

<ul> After alignment almost all or a very large portion of the reads (>80%) are unmapped or the number of non-zero gene counts is very high.. 

 Possible reasons are: wrong target genome/transcriptome. 

 💡 -  </ul> 


### 🚩 Unmapped Reads (Mapping)

<ul> high number of unmapped reads. 

 Possible reasons are: adapter sequence still present. 

 💡 go back to trimming  </ul> 


### 🚩 Correlation between biological conditions (Post-Mapping)

<ul> different biological conditions correlate better than the replicates with each other. 

 Possible reasons are: low potential for DEGs present. 

 💡 There is no general solution.  </ul> 


### 🚩 Inconsistent Replicates (Post-Mapping)

<ul> within replicate correlation is weak (r2<0.5). 

 Possible reasons are: e.g. contamination with HVG like rRNAs. 

 💡 if only a small sub-population is confounding -> remove  </ul> 

<ul> replicates do not cluster. 

 Possible reasons are: e.g. biol. and tech. replicates are mixed up. 

 💡 END-RESTART  </ul> 

<ul> replicates do not cluster. 

 Possible reasons are: negative control study. 

 </ul> 


### 🚩 Correlation between biological conditions (Post-Mapping)

<ul> two experiment of 2 conditions (ctrl, KO) cluster in the wrong area. 

 Possible reasons are: human error, e.g. misslabeling of tubes. 

 💡 check design matrix, documentation, ...  </ul> 

<ul> no clear separation between conditions (ctrl, KO) . 

 Possible reasons are: e.g. biol. and tech. replicates are mixed up. 

 💡 END-RESTART  </ul> 


### 🚩 Dispersion Plot Anomalies (DEA)

<ul> dispersion-plot: gene estimation does not follow red fit. 

 Possible reasons are: model does not represent data well. 

 💡 Deseq is not applicable  </ul> 

<ul> dispersion-plot: high fit dispersion for high mean count. 

 Possible reasons are: low number of replication + high variability. 

 💡 careful with reported DEGs  </ul> 


### 🚩 Predominantly Differentially Expressed Genes (DEA)

<ul> almost all genes are DE (>90%). 

 Possible reasons are: e.g. biol. and tech. replicates are mixed up. 

 💡 END-RESTART  </ul> 


### 🚩 Sparse Differentially Expressed Genes (DEA)

<ul> almost no gene is DE (<5%). 

 Possible reasons are: negative control study. 

 </ul> 


### 🚩 Wrong visualization metric (Visualize)

<ul> using CPM to visualize DEGs. 

 Possible reasons are: mixed up within- and between-sample normalization. 

 💡 use correct normalization  </ul> 

<ul> use raw Ratios to visualize DEGs effect size. 

 Possible reasons are: humans are bad with ratios (0.01 = almost 0 and 100 is just large but not the largest bar ever). 

 💡 use any log transformation (e.g. log10: 0.01 => -2, 100 => +2)  </ul> 


# single cell


### 🚩 Problem with per cell histogram (quality check)

<ul> peak at left/right side in gene or reads per cell histogram or log10-cummulative-number of reads per cell id. 

 Possible reasons are: left=cell fragments, right=multiplets present. 

 💡 remove with cut-offs  </ul> 


### 🚩 inflated UMI counts (quality check)

<ul> inflated UMI counts. 

 Possible reasons are: using raw UMI counts. 

 💡 use DBEC/RSEC UMI counts  </ul> 


### 🚩 Suboptimal embedding (dimension reduction)

<ul> poor PCA/UMAP/tSNE embedding. 

 Possible reasons are: using e.g. only one assay of count data for embedding. 

 💡 for multimodal data use a WNN approach (combining both assays)  </ul> 

<ul> poor PCA/UMAP/tSNE embedding. 

 Possible reasons are: use raw data for tSNE/UMAP. 

 💡 use a significant portion of PC from the PCA as input for tSNE/UMAP  </ul> 


### 🚩 Wrong embeddings (dimension reduction)

<ul> "The first two principle components were used to perform a tSNE" (https://www.science.org/doi/full/10.1126/science.aag3009). 

 Possible reasons are: use only 2 PC from the PCA for the tSNE/UMAP projection. 

 💡 use a significant portion of PC from the PCA as input for tSNE/UMAP  </ul> 


### 🚩 Batch problems (find subpopulations)

<ul> there are batch specific cluster. 

 Possible reasons are: batch effect. 

 💡 correct for batch effect (e.g. integrate using seurat)  </ul> 


### 🚩 Abundance of DEGs (DEA)

<ul> many DEGs. 

 Possible reasons are: DEG between very small sub-populations. 

 💡 use a population size cutoff or state the number  </ul> 


### 🚩 Gene Expression Resembling Dates (DEA)

<ul> some genes look like dates (1-Mar,...). 

 Possible reasons are: some genes can be interpreted as dates when using excel for data handling (https://www.science.org/doi/10.1126/science.aah4573). 

 💡 never ever use excel or at least make sure that cell type is not "AUTO"  </ul> 


# rtqpcr


### 🚩 Non-Sigmoidal Amplification Plot 

<ul> Amplification plot is not a single sigmoidal like curve. 

 Possible reasons are: insufficient template concentration. 

 💡 increase template concentration / improve RT efficiency – try different RT or additives like Mn2+  </ul> 


### 🚩 Multiple Products in Melting Curve 

<ul> multiple products occur in the melting curve following the qPCR. 

 Possible reasons are: unspecific primer pair. 

 💡 test different primers – increase annealing temp  </ul> 


### 🚩 Inconsistent Replicates 

<ul> Great variation withing technical replicates (>10%). 

 Possible reasons are: variable pipetting – bad input RNA quality. 

 💡 check pipette seals – use electronic dispenser – use mastermix if possible – always check input RNA quality beforehand  </ul> 


### 🚩 Product Detected in -RT Control 

<ul> Low Ct value in -rt control (<40). 

 Possible reasons are: incomplete DNA digestion. 

 💡 redo DNA digestion  </ul> 


### 🚩 Variable Expression of Normalization Genes 

<ul> inconsisten expression of normalization genes within the experimental series. 

 Possible reasons are: bad gene choice for this condition. 

 💡 switch to different normalization gene – use several (3+) genes for normalization in parallel  </ul> 

