**1. Data trimming of Illumina 150bp PE reads**

```bash
fastp -i input_1.fastq.gz -I input_2.fastq.gz -o output_clean_1.fq.gz -O output_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
```

**2. Mapping trimmed reads to the reference geneome of *C. angulata***

```bash
for i in `cat sample.list`
do
fastq-dump --split-files --gzip $i.sra && \
fastp -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ${i}_clean_1.fastq.gz -O ${i}_clean_2.fastq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 ${i}_um_1.fastq.gz —unpaired2 ${i}_um_2.fastq.gz —failed_out ${1}_failed.fastq.gz --cut_front --cut_front_window_size=1 --cut_front_mean_quality=3 --cut_tail --cut_tail_window_size=1 --cut_tail_mean_quality=3 --cut_right --cut_right_window_size=4 --cut_right_mean_quality=15 --length_required=36 --thread 2 --trim_front1 10 --trim_front2 10 --trim_poly_x && \
rm -f ${i}_{1,2}.fastq.gz && \
hisat2 -q -x ref/hisat2_index/genome_tran -p 2 -1 ${i}_clean_1.fastq.gz -2 ${i}_clean_2.fastq.gz -S hisat2.$i.sam --dta-cufflinks --summary-file hisat2.$i.log && \
samtools view -Sb hisat2.$i.sam | samtools sort -O BAM > hisat2.$i.bam && \
rm -f hisat2.$i.sam && \
samtools index -@ 4 hisat2.$i.bam && \
rm -f ${i}_clean_{1,2}.fastq.gz
done
```

**Validation of the sample information**

```bash
for i in `cut -f 1 ../sample.info`
do 
ln -s ../hisat2/$i/hisat2.$i.bam* ./
done

# snp calling
for e in `cat chr.list`
do
echo REF=/public1/ref/c.gigas/roslin/ref/ref.fa > variantcalling.$e.sh
echo bcftools mpileup --threads 6 -r $e -a AD,DP,SP -Ou -f \$REF \*.bam \| bcftools call --threads 6 -f GQ,GP -mO z -o $e.vcf.gz >> variantcalling.$e.sh
nohup time sh variantcalling.$e.sh > variantcalling.$e.log 2>&1&
done

# filter
for e in `cat chr.list`
do
echo VCF_IN=$e.vcf.gz > filter.$e.sh
echo VCF_OUT=$e.filtered.vcf.gz >> filter.$e.sh
echo MAF=0.1 >> filter.$e.sh
echo MISS=0.9 >> filter.$e.sh
echo QUAL=30 >> filter.$e.sh
echo MIN_DEPTH=10 >> filter.$e.sh
echo MAX_DEPTH=50 >> filter.$e.sh
echo # perform the filtering with vcftools >> filter.$e.sh
echo vcftools --gzvcf \$VCF_IN \\ >> filter.$e.sh
echo --remove-indels --maf \$MAF --max-missing \$MISS --minQ \$QUAL \\ >> filter.$e.sh
echo --min-meanDP \$MIN_DEPTH --max-meanDP \$MAX_DEPTH \\ >> filter.$e.sh
echo --minDP \$MIN_DEPTH --maxDP \$MAX_DEPTH --recode --stdout \| gzip -c \> \\ >> filter.$e.sh
echo \$VCF_OUT >> filter.$e.sh
nohup time sh filter.$e.sh > filter.$e.log 2>&1&
done

# merge 
bcftools concat -O z -o filtered.vcf.gz NC_04*.filtered.vcf.gz

# filter lowDP
vcftools --gzvcf filtered.vcf.gz --missing-indv > indv.miss

awk '{if($5>0.2)print $1}' out.imiss > lowDP.indv

vcftools --gzvcf filtered.vcf.gz --remove lowDP.indv --recode --stdout | gzip -c > filterlowDP.vcf.gz

# perform linkage pruning - i.e. identify prune sites
time plink --vcf filterlowDP.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out cg.vs.ca
# PCA output:
#	•	cichlids.eigenval - the eigenvalues from our analysis
#	•	cichlids.eigenvec- the eigenvectors from our analysis
# plink binary output
#	•	cichlids.bed - the cichlids bed file - this is a binary file necessary for admixture analysis. It is essentially the genotypes of the pruned dataset recoded as 1s and 0s.
#	•	cichlids.bim - a map file (i.e. information file) of the variants contained in the bed file.
#	•	cichlids.fam - a map file for the individuals contained in the bed file.

# prune and create pca
time plink --vcf filterlowDP.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract cg.vs.ca.prune.in \
--threads 10 \
--make-bed --pca --out cg.vs.ca

for i in `awk '{print $1":"$3}' ../sample.info`; do
sed -i "s/hisat2.${i%%:*}.bam/${i##*:}/g" cg.vs.ca.eigenvec
done
```
```R
# load tidyverse package
library(tidyverse)
# read in data
pca <- read.table("./cg.vs.ca.eigenvec", header = FALSE, sep=" ")
eigenval <- scan("./cg.vs.ca.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# location
loc <- rep(NA, length(pca$ind))
loc[grep("QD_AN", pca$ind)] <- "cang_n"
loc[grep("QD_QD", pca$ind)] <- "cgig_n"
loc[grep("XM_AN", pca$ind)] <- "cang_s"
loc[grep("XM_QD", pca$ind)] <- "cgig_s"
# species
spe <- rep(NA, length(pca$ind))
spe[grep("QD_AN", pca$ind)] <- "cang"
spe[grep("QD_QD", pca$ind)] <- "cgig"
spe[grep("XM_AN", pca$ind)] <- "cang"
spe[grep("XM_QD", pca$ind)] <- "cgig"
# remake data.frame
pca <- as_tibble(data.frame(pca, loc))
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

library(ggrepel)

pdf("pca.pdf",width=4,height=3)

ggplot(pca, aes(x=PC1, y=PC2, label=ind, shape=spe, color=loc)) +
    geom_point(size = 2) +
#    geom_text_repel(aes(label =ind),size = 2.5) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    xlab(paste("PC1 (",round(pve[1,2],2),"%)",sep="")) + 
    ylab(paste("PC2 (",round(pve[2,2],2),"%)",sep="")) +
    scale_color_brewer(palette="Dark2") + 
    scale_shape_manual(values = c(15, 16)) +
    theme_minimal() + 
    stat_ellipse()

dev.off()

```
<img src="https://github.com/user-attachments/assets/fe69e21a-3b5e-418b-a9cc-248c9b128762" alt="PCA" width="500"/>

<br><br>
**3. Quantification of reads aligned to genomic features and readcount to FPKM**

```bash
for i in `cat sample.list`
do 
htseq-count -s no -r pos -f bam hisat2.$i.bam ref.gtf > $i.counts
done
```


**Box plot of gene expression level**

```R
# readcount to FPKM
df <- read.table("allSamples.counts", header=T, row.names=1)
totalMappedReads <- colSums(df)

geneLengths <- read.table("gene_length.txt",header=F,row.names=1)
geneLengths <- geneLengths[rownames(df),]
geneLengthsKb <- geneLengths / 1000

fpkm <- sweep(df, 2, totalMappedReads, "/")
fpkm <- sweep(fpkm, 1, geneLengthsKb, "/")
fpkm <- fpkm * 1e6

fpkm[] <- lapply(fpkm, function(x) if(is.numeric(x)) round(x, 2) else x)
write.table(fpkm, "allSamplesFPKM.txt",quote=F,sep="\t",row.names=T,col.names=T)
```
```R
library(ggpubr)

df <- read.table("allSamplesFPKM.txt", header=T, row.names=1)
info <- read.table("sample.info", header=T, row.names=1)

# Subset the expression data for LOC128165822
gene_data <- df["LOC128165822", ]

# Convert gene_data to a data frame and melt it for ggplot2
library(reshape2)
gene_data_melted <- melt(gene_data)

# Rename columns for clarity
colnames(gene_data_melted) <- c("Sample", "Expression")
rownames(gene_data_melted) <- gene_data_melted$Sample

# Merge with sample information (info) to add species and location
merged_data <- merge(gene_data_melted, info, by = "row.names", all.x = TRUE)

# Create the boxplot
my_comparisons <- list( c("C.angulata", "C.gigas"))
pdf("LOC128165822.pdf", width=4, height=3)
ggboxplot(merged_data, x = "Species", y = "Expression",
          color = "black", palette = "jco", add = "jitter", facet.by = "Location")+ 
  stat_compare_means(comparisons = my_comparisons) + theme(legend.position = "none") + xlab("")
dev.off()

# Subset the expression data for LOC128168260
gene_data <- df["LOC128168260", ]

# Convert gene_data to a data frame and melt it for ggplot2
library(reshape2)
gene_data_melted <- melt(gene_data)

# Rename columns for clarity
colnames(gene_data_melted) <- c("Sample", "Expression")
rownames(gene_data_melted) <- gene_data_melted$Sample

# Merge with sample information (info) to add species and location
merged_data <- merge(gene_data_melted, info, by = "row.names", all.x = TRUE)

# Create the boxplot
my_comparisons <- list( c("C.angulata", "C.gigas"))
pdf("LOC128168260.pdf", width=4, height=3)
ggboxplot(merged_data, x = "Species", y = "FPKM",
          color = "black", palette = "jco", add = "jitter", facet.by = "Location")+ 
  stat_compare_means(comparisons = my_comparisons) + theme(legend.position = "none") + xlab("")
dev.off()
```
**SCD1, LOC128165822**

<img src="https://github.com/user-attachments/assets/cb877917-2b88-4a5e-94ff-2554a851e31a" alt="LOC128165822" width="500"/>


**SCD2, LOC128168260**

<img src="https://github.com/user-attachments/assets/043eed67-cfa0-4d3d-b03a-21f3f3f15d8f" alt="LOC128168260" width="500"/>
<br><br>

**4. Differential expression analysis**

```R
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
directory <- "./cang0h.vs.cang24h"
sampleFiles <- list.files(directory, pattern = ".*counts", full.names = TRUE)
sampleNames <- c("QD_AN_0_3","QD_AN_0_1","QD_AN_0_2","QD_AN_24_1","QD_AN_24_2","QD_AN_24_3","XM_AN_0_3","XM_AN_24_3","XM_AN_24_2","XM_AN_24_1","XM_AN_0_1","XM_AN_0_2")
sampleCondition <- c("heat0","heat0","heat0","heat24","heat24","heat24","heat0","heat24","heat24","heat24","heat0","heat0")
sampleTable <- data.frame(
		sampleName = sampleNames,
		fileName = sampleFiles,
		condition = sampleCondition
		)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(
		sampleTable = sampleTable, 
		directory = "./", 
		design = ~ condition
		)

ddsHTSeq$condition
treatments <- c("heat0","heat24")
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = treatments)
sumcounts <- rowSums(counts(ddsHTSeq))


# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# you can see the typically high dynamic range of RNA-Seq, with a mode in the distribution around 1000 fragments per gene, but some genes up over 1 million fragments. 

# get genes with summed counts greater than 20
keep <- sumcounts > 20

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
res <- results(dds)
write.table(res, "cang0h.vs.cang24h.result", quote=F, sep="\t")
```

**5. GO Ontology and visualization**

```bash
library(clusterProfiler)

#加载背景库文件
go_anno <- read.delim('gene2go.txt', header=FALSE, stringsAsFactors =FALSE)
names(go_anno) <- c('gene_id','ID')

#加载GO注释描述
go_class <- read.delim('go-basic.txt', header=FALSE, stringsAsFactors =FALSE)
names(go_class) <- c('ID','Description','Ontology')

#合并背景与GO文件
go_anno <-merge(go_anno, go_class, by = 'ID', all.x = TRUE)

#差异基因导入
gene_list <- read.delim('cang0h.vs.cang24h.padj05.up.list',stringsAsFactors = FALSE)
names(gene_list) <- c('gene_id')
gene_select <- gene_list$gene_id

#富集分析
go_rich <- enricher(gene = gene_select,
                 TERM2GENE = go_anno[c('ID','gene_id')],
                 TERM2NAME = go_anno[c('ID','Description')],
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'BH',
                 qvalueCutoff = 0.2,
                 maxGSSize = 200)

barplot(go_rich,showCategory = 20,drop=T)
dev.off()
write.table(go_rich, 'cang0h.vs.cang24h.go_rich.txt', sep='\t', row.names = FALSE, quote = FALSE)


library(aPEAR)
library(clusterProfiler)
library(ggplot2)
library(cols4all)

pdf("ca0hvs24h.padj05.up.significant.pdf")

df <- read.table("ca0hvs24h.padj05.up.significant.txt",header=T,sep="\t")
bio <- df[df$Ontology == "biological_process",]
enrichmentNetwork(bio, 
 colorBy="qvalue", 
 colorType=c("pval"), 
 nodeSize="Count", 
 fontSize=4, 
 drawEllipses=F, 
 pCutoff=-19, 
 verbose=T) +
 scale_color_continuous_c4a_seq("wes.zissou1", reverse = T)

dev.off()


df <- read.table("ca0hvs24h.padj05.down.significant.txt",header=T,sep="\t")
bio <- df[df$Ontology == "biological_process",]
pdf("ca0hvs24h.padj05.down.significant.pdf")
enrichmentNetwork(bio, 
 colorBy="qvalue", 
 colorType=c("pval"), 
 nodeSize="Count", 
 fontSize=4, 
 drawEllipses=F, 
 pCutoff=-16, 
 verbose=T) +
 scale_color_continuous_c4a_seq("wes.zissou1", reverse = T)
#scale_color_continuous_c4a_seq('yl_gn_bu', reverse = T)

dev.off()
```
**Upregulated gene set**

<img src="https://github.com/user-attachments/assets/265ad09c-8fae-46ff-b564-f84940e708ef" alt="PCA" width="500"/>

<br><br>
**Downregulated gene set**

<img src="https://github.com/user-attachments/assets/ac1a1abf-9ff1-4f16-bec7-997e111e45c0" alt="PCA" width="500"/>
