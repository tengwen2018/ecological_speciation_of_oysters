**1. Data trimming of Illumina 150bp PE reads**

```bash
fastp -i input_1.fastq.gz -I input_2.fastq.gz -o output_clean_1.fq.gz -O output_clean_2.fq.gz —adapter_sequence auto —detect_adapter_for_pe —unpaired1 output_um_1.fastq.gz —unpaired2 output_um_2.fastq.gz —failed_out output_failed.fastq.gz —cut_front —cut_front_window_size=1 —cut_front_mean_quality=20 —cut_tail —cut_tail_window_size=1 —cut_tail_mean_quality=20 —cut_right —cut_right_window_size=4 —cut_right_mean_quality=20 —length_required=36 —thread 1 --trim_front1 5 --trim_front2 5
```

**2. Mapping trimmed reads to the reference geneome of *C. angulata***

```bash
bwa index ref.fa
for i in `cat sample.list`
do
bwa mem -M -t 4 ref.fa ${i}_clean_1.fq.gz ${i}_clean_2.fq.gz | samtools view -bS > $i.bam
# mark duplications with samtools v1.7
samtools sort -@ 8 -n -o namesort.bam $i.bam && \
rm -f $i.bam && \
samtools fixmate -@ 8 -m namesort.bam fixmate.bam && \
rm -f namesort.bam && \
samtools sort -@ 8 -o positionsort.bam fixmate.bam && \
rm -f fixmate.bam && \
samtools markdup -@ 8 -r positionsort.bam $i.dup.bam && \
rm -f positionsort.bam && \
samtools index -@ 8 $i.dup.bam
# remove dup and low quality mapped reads
samtools view -@ 8 -h -F 0x100 -F 0x400 $i.dup.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 8 -q 30 -bS > $i.bam && \
samtools index $i.bam && \
rm -f $i.dup.bam
done
```

**3. SNP calling and filtering**

```bash
bcftools mpileup --threads 10 -a AD,DP,SP -Ou -f ref.fa *.bam | bcftools call --threads 10 -f GQ,GP -mO z -o oyster.vcf.gz
VCF_IN=oyster.vcf.gz
VCF_OUT=oyster.filtered.vcf.gz
MAF=0.1
MISS=0.8
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT
```

**4. Estimation of per-site FST and sliding-window *F*<sup>ST</sup>**

```bash
# filter individuals with missingness above 20%
vcftools --gzvcf oyster.filtered.vcf.gz --missing-indv
sed '1d' out.imiss | awk '{if($5>0.2)print $1}' > lowDP.indv
vcftools --gzvcf oyster.filtered.vcf.gz --remove lowDP.indv --recode --stdout | gzip -c > filterlowDP.vcf.gz
# per site fst
time vcftools --gzvcf filterlowDP.vcf.gz \
--weir-fst-pop cang \
--weir-fst-pop cgig \
--out ./cang.vs.cgig
# sliding windows of 2kb
parseVCF.py -i filterlowDP.vcf.gz | bgzip > filterlowDP.geno.gz
popgenWindows.py -g filterlowDP.geno.gz -o filterlowDP.2k.Fst.Dxy.pi.csv.gz \
   -f phased -w 2000 -m 10000 -s 2000 \
   -p cang -p cgig \
   --popsFile pop.info
```
