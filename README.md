# Lysine acetylation of aquaporin-3 promotes kidney water reabsorption but is not essential for concentrating urine during water deprivation
<p align="center"> Nha Van Huynh and Kelly A. Hyndman <br>
University of Alabama at Birmingham <br></c>

Kidney bulk transcriptomic profiles of AQP3 K282 acetylation and deacetylation mimetic mice under standard housing conditions. Male and female C57BL/6J wild-type (K282), acetylation (K282Q or Q) and deacetylation (K282R or R) mimetic mutant mice were group housed with ad libitum access to standard pelleted chow diet and water. Euthanasia and sample collection were done at 7 – 10 AM (Zeitgeber time 1-4). Kidney were collected, decapsulated, and cut transversely. Kidney halves were flash frozen in liquid nitrogen and stored in -80ºC until RNA extraction. Total RNA was extracted using Direct-Zol RNA miniprep kit (Zymo Research, Irvine, CA) from half of a kidney/animal and lyophilized in RNA stabilization tubes (Azenta LifeSciences, Burlington, MA) before being sent to Azenta for library preparation and RNA sequencing. cDNA libraries from Poly(A)+ RNA were prepared and next-generation sequencing with paired-end format targeting 20 million reads per sample was done. 

The fastq files were trimmed with Trim_Galore v0.6.7. The reads were mapped to the mouse genome mm39 using STAR v2.7.10a. BAI files were created using SAMtools v1.6. The raw counts were calculated using featureCounts function of the Subread v2.0.1. 

Differentially expressed genes (DEGs) were determined using R v4.2.1 and DESeq2 package v1.38.3. The whole dataset (samples were grouped as Male WT, Male Q, Male R, Female WT, Female Q, Female R) was used to generate dds object. Genes that had total counts from all samples < 100 were excluded before running DESeq2. The following pairwise comparisons were made: male Q versus male WT, male R versus male WT, male Q versus male R, female Q versus female WT, female R versus female WT, female Q versus female R. No threshold was set for Log2 fold change (Log2FC), while adjusted p-values (padj) < 0.05 were considered significant for DEG identification.

Raw files and processed files from bulk RNA sequencing can be found on Gene Expression Omnibus (GEO Accession number GSE279083). 

Publication: 
