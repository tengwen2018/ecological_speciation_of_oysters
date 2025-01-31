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

**3. Quantification of reads aligned to genomic features and readcount to FPKM**

```bash
for i in `cat sample.list`
do 
htseq-count -s no -r pos -f bam hisat2.$i.bam ref.gtf > $i.counts
done
```

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
