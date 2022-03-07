# Accumulation of endosymbiont genomes in an insect autosome followed by endosymbiont replacement

Eric S. Tvedte

Last updated 2022-03-07



## Table of Contents
1. [Download D. ananassae and wAna datasets](#dl)
2. [Nuwt identification and quantification](#nuwt_id)
3. [Numt analysis](#numt)
4. [RNAseq analysis](#lgt.tx) 


### 1. Download D. ananassae and wAna datasets <a name="dl"></a>
**Reference genomes**
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_genomic.fna.gz
gunzip *.gz  
D. ananassae Hawaii PacBio Sequel II assembly available on Figshare https://doi.org/10.25387/g3.14096897
```
**Read datasets**
```
D. ananassae Hawaii datasets available on NCBI SRA under BioProject PRJNA602597
#geographic D. ana strains
fastq-dump --split-files SRR2126857
fastq-dump --split-files SRR2126916
fastq-dump --split-files SRR2127151
fastq-dump --split-files SRR2127152
fastq-dump --split-files SRR2127153
fastq-dump --split-files SRR2127154
fastq-dump --split-files SRR2127155
fastq-dump --split-files SRR2127156
fastq-dump --split-files SRR2127161
fastq-dump --split-files SRR2127162
fastq-dump --split-files SRR2127163
fastq-dump --split-files SRR2127164
fastq-dump --split-files SRR2127219
fastq-dump --split-files SRR2135551
fastq-dump --split-files SRR2135600
#wAna
fastq-dump --split-files SRR8278850
```

### 2. Nuwt identification and quantification <a name="nuwt_id"></a>

**Identify nuwt contigs**
```
mummer/nucmer -l 1000 --prefix nuwt.firstpass GCA_008033215.1_ASM803321v1_genomic.fna Dana.UMIGS.fasta
mummer/show-coords -rT nuwt.firstpass.delta > nuwt.firstpass.coords
tail -n +5 nuwt.firstpass.coords | awk '{print $9}' | sort -n | uniq > nuwt.contigs.list
seqkit grep Dana.UMIGS.fasta -f nuwt.contigs.list > Dana.UMIGS.nuwt.contigs.fasta
```

**Aligment of wAna to nuwt contigs**
```
~jdhotopp/bin/residues.pl Dana.UMIGS.nuwt.contigs.fasta > Dana.UMIGS.nuwt.contigs.residues
nucmer -l 1000 --prefix nuwt.finalpass GCA_008033215.1_ASM803321v1_genomic.fna Dana.UMIGS.nuwt.contigs.fasta
mummer/delta-filter -q nuwt.finalpass.delta > nuwt.finalpass.filter
data_viz_scripts/Dana.LGT.Rmd
```

**Estimate total nuwt content in D. ananassae**
```
mummer/show-coords -rT nuwt.finalpass.filter > nuwt.finalpass.coords
tail -n +5 nuwt.finalpass.coords | awk '{print $9"\t"$3"\t"$4}' > nuwt.finalpass.bed
Rscript fixbed.R nuwt.finalpass.bed nuwt.finalpass.fixed.bed
awk '{print $1"\t"$2-1"\t"$3}' nuwt.finalpass.fixed.bed > nuwt.finalpass.final.bed #NUCmer start coodinates are 1-based, BEDtools 0-based 
bedtools coverage -a nuwt.finalpass.fixed.bed -b nuwt.finalpass.fixed.bed -hist | grep 'all' > nuwt.finalpass.fixed.hist.out 
#estimated nuwt length is 1*1 depth + (1/2)*2 depth in hist.out file, this corrects for small overlapping segments generated using NUCmer
```

**Estimate nuwt copy number and HiFi sequencing depth**
```
tail -n +5 nuwt.finalpass.coords | awk '{print $8"\t"$1-1"\t"$2}' > wAna.finalpass.bed #NUCmer start coodinates are 1-based, BEDtools 0-based 
touch wAna.genome.bed #enter BED coordinates for whole genome - CP042904.1 0 1401460
bedtools coverage -a wAna.genome.bed -b wAna.finalpass.bed -d > wAna.genome.NUCmer.cn.bed

minimap2 -ax map-pb -t 16 GCA_008033215.1_ASM803321v1_genomic.fna PB.HiFi.fastq.gz | samtools sort -o GCA_008033215.1_ASM803321v1_mapped_HiFi_sorted.bam  
samtools depth -a GCA_008033215.1_ASM803321v1_mapped_HiFi_sorted.bam > GCA_008033215.1_ASM803321v1_mapped_HiFi_sorted.depth.txt
data_viz_scripts/Dana.LGT.Rmd
```

**Purge_haplotigs: purge/clip potential duplicated assembly sequence**
```
minimap2 -t 4 -ax map-pb Dana.UMIGS.fasta PB.HiFi.fastq.gz --secondary=no | samtools sort -m 5G -o aln.bam
purge_haplotigs hist -b aln.bam -g pilon.fasta -t 4 -d 400
purge_haplotigs cov -l 5 -m 195 -h 300 -i aln.bam.gencov -j 80 -s 80
purge_haplotigs purge -g Dana.UMIGS.fasta -b aln.bam -c coverage_stats.csv -d -t 8
purge_haplotigs clip -p contigs.fasta -h haplotigs.fasta -l 20000 -g 10000 -t 4
```

**Generate depth histograms for wAna CDS across D. anananassae strains**
```
#generate CDS BED file
gffread -C GCA_008033215.1_ASM803321v1_genomic.gff --bed > wAna.CDS.bed
#Dana cured Hawaii PacBio HiFi 
minimap2 -ax map-pb -t 16 GCA_008033215.1_ASM803321v1_genomic.fna PB.HiFi.fastq.gz | samtools sort -o Dana.UMIGS_mapped_HiFi_primary.bam
samtools depth -a -b wAna.CDS.bed -m 100000000 > samt.dp.txt 
awk '{print $3}' sam.dp.txt | sort -n | uniq -c | awk '{print $1"\t"$2}' > sam.dp.hist.txt
#Dana cured/infected Illumina 
java -Xmx20g -jar trimmomatic-0.38.jar PE -phred33 R1.fastq R2.fastq paired_R1.fastq.gz unpaired_R1.fastq.gz paired_R2.fastq.gz unpaired_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
bwa mem -k 23 -t 12 GCA_008033215.1_ASM803321v1_genomic.fna paired_R1.fastq.gz paired_R2.fastq.gz | samtools view -bho output.bam
samtools sort -o sorted.bam output.bam 
java -Xmx20g -jar picard.jar MarkDuplicates I=sorted.bam O=dedup.bam M=dups.metrics.txt REMOVE_DUPLICATES=TRUE ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
samtools depth -a -b wAna.CDS.bed -m 100000000 dedup.bam > samt.dp.txt
awk '{print $3}' sam.dp.txt | sort -n | uniq -c | awk '{print $1"\t"$2}' > sam.dp.hist.txt
```

**Generate depth histograms for D. ananassae regions across D. anananassae strains**
```
#generate CDS BED file
gffread -C GCA_008033215.1_ASM803321v1_genomic.gff --bed > wAna.CDS.bed
#Dana cured Hawaii PacBio HiFi
minimap2 -ax map-pb -t 16 Dana.UMIGS.FREEZE.fasta PB.HiFi.fastq.gz | samtools sort -o Dana.UMIGS_mapped_HiFi_primary.bam
samtools depth -a -b chr2L.bed -m 100000000 > samt.chr2L.dp.txt #chr2L
samtools depth -a -b chr4.bed -m 100000000 > samt.chr4.dp.txt #chr4
awk '{print $3}' sam.dp.txt | sort -n | uniq -c | awk '{print $1"\t"$2}' > sam.dp.hist.txt
#Dana cured/infected Illumina 
bwa mem -k 23 -t 12 Dana.UMIGS.FREEZE.fasta paired_R1.fastq.gz paired_R2.fastq.gz | samtools view -bho output.bam
samtools sort -o sorted.bam output.bam 
java -Xmx20g -jar picard.jar MarkDuplicates I=sorted.bam O=dedup.bam M=dups.metrics.txt REMOVE_DUPLICATES=TRUE ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
samtools depth -a -b chr2L.bed -m 100000000 dedup.bam > samt.chr2L.dp.txt
samtools depth -a -b chr4.bed -m 100000000 dedup.bam > samt.chr4.dp.txt
awk '{print $3}' sam.dp.txt | sort -n | uniq -c | awk '{print $1"\t"$2}' > sam.dp.hist.txt
```

### 3. NUMT analysis <a name="numt"></a>
**Assemble mitochondrial genome using NOVOPlasty**
```
perl NOVOPlasty3.7.pl -c dana.mito.config.txt
smrttools/smrtcmds/bin/pbmm2 align Circularized_assembly_mito.fasta PB.HiFi.ccs.bam Dana.UMIGS.dana.mito.mapped.HiFi_sorted.bam --sort -j 16 -J 8
smrttools/smrtcmds/bin/arrow -j 16 Dana.UMIGS.dana.mito.mapped.HiFi_sorted.bam -r Circularized_assembly_mito.fasta -o Dana.UMIGS.mito.hifi_variants.gff -o Dana.UMIGS.mito.polished.fasta
```

**Search for numts by alignments between mito genome and nuclear genome**
```
nucmer -l 100 --maxmatch --prefix dana.numt.firstpass Dana.UMIGS.mito.rotate.FREEZE.fasta Dana.UMIGS.fasta 
show-coords -rT dana.numt.firstpass.delta > dana.numt.firstpass.coords  
tail -n +5 numt.firstpass.coords | awk '{print $9}' | sort -n | uniq > numt.contigs.list
seqkit grep Dana.UMIGS.fasta -f numt.contigs.list > Dana.UMIGS.numt.contigs.fasta
samtools faidx Dana.UMIGS.fasta tig00000054:4380000-4440000 > tig00000054.numt.fasta
```

**Align mitochondrial genome to large numt in chromosome 4 contig tig00000054**
```
~jdhotopp/bin/residues.pl tig00000054.numt.fasta > tig00000054.numt.residues
nucmer -l 100 --prefix numt.finalpass Dana.UMIGS.mito.rotate.FREEZE.fasta tig00000054.numt.region.fasta
mummer/delta-filter -q numt.finalpass.delta > numt.finalpass.filter
data_viz_scripts/Dana.LGT.Rmd
```
**Determine PacBio CLR and ONT reads overlapping numt region**
```
minimap2 -ax map-pb -t 16 Dana.UMIGS.contigs.fasta PB.CLR.fastq.gz | samtools sort -o Dana.UMIGS_mapped.PB.CLR_sorted.bam 
bedtools intersect -a Dana.UMIGS_mapped_PB_CLR_sorted.bam -b numt.final.region.bed -F 1 -sorted > ONT.numt.region.bam
samtools view -F 256 -c ONT.numt.region.bam
minimap2 -ax map-ont -t 16 Dana.UMIGS.contigs.fasta ONT.fastq.gz | samtools sort -o Dana.UMIG.mapped.ONT_sorted.bam  
bedtools intersect -a Dana.UMIGS_mapped_ONT_sorted.bam -b numt.final.region.bed -F 1 -sorted > PB.CLR.numt.region.bam
samtools view -F 256 -c PB.CLR.numt.region.bam
bedtools bamtobed -i PB.CLR.numt.region.bam > PB.CLR.numt.region.bed
```

### 4. RNAseq analysis <a name="lgt.tx"></a>
**Download datasets from SRA**
```
#SRA IDs in Table SX
fastq-dump SRA_ID #single end
fastq-dump --split-files SRA_ID #paired end
```
**De-duplicate newly generated RNAseq**
```
bbtools/bbmerge.sh verystrict=t in=illumina_R1.fastq.gz in2=illumina_R2.fastq.gz out=illumina.vstrict.merged.fastq.gz outu=illumina.vstrict.unmerged.fastq.gz outa=illumina.vstrict.adapters.fa 
bbtools/clumpify.sh dedupe=t optical=f in=illumina.vstrict.merged.fastq.gz out=illumina.vstrict.merged.dedup.fastq.gz 
```
**Construct salmon index**
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/285/975/GCF_003285975.2_DanaRS2.1/GCF_003285975.2_DanaRS2.1_genomic.gff.gz
gunzip *gz
gffread -C GCA_008033215.1_ASM803321v1_genomic.gff -g GCA_008033215.1_ASM803321v1_genomic.fna -x wAna.CDS.fasta
sed -i 's/gene-//g' wAna.CDS.fasta
gffread -C GCF_003285975.2_DanaRS2.1_genomic.gff -g GCF_003285975.2_DanaRS2.1_genomic.fna -x Dana.RS2.CDS.fasta
sed -i 's/rna-//g' Dana.RS2.CDS.fasta
cat wAna.CDS.fasta Dana.RS2.CDS.fasta GCA_008033215.1_ASM803321v1_genomic.fna Dana.UMIGS.fasta > combined.CDS.genomes.fasta
grep "^>" GCA_008033215.1_ASM803321v1_genomic.fna | cut -d " " -f 1 > decoys.txt
grep "^>" Dana.UMIGS.fasta | cut -d " " -f 1 >> decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
salmon index -t combined.CDS.genomes.fasta -d decoys.txt -i salmon_index -k 25
```
**Run salmon quantification in mapping-mode with Illumina datasets**
```
#single end
salmon quant -i salmon_index -l A -r R1.fastq.gz --validateMappings --minAssignedFrags 1 -o salmon_output -p 8 --writeMappings | samtools view -bhSo mapped.transcripts.output.bam
#paired end
salmon quant -i salmon_index -l A -1 R1.fastq.gz -2 R2.fastq.gz --validateMappings --minAssignedFrags 1 -o salmon_output -p 8 --writeMappings | samtools view -bhSo mapped.transcripts.output.bam

```
**Run salmon quantification in alignment-mode with ONT datasets**
```
cat wAna.CDS.fasta Dana.RS2.CDS.fasta combined.CDS.fasta
minimap2 -ax map-ont combined.CDS.fasta cHI_ONT.direct.rna.fastq | samtools view -bho cHI_directRNA_mapped.combined.CDS_output.bam #cured Drosophila
minimap2 -ax map-ont combined.CDS.fasta WT_ONT.direct.rna.fastq | samtools view -bho WT_directRNA_mapped.combined.CDS_output.bam #wild-type Drosophila
salmon quant -l SF -t combined.CDS.fasta -a cHI_directRNA_mapped.combined.CDS_output.bam --noErrorModel --minAssignedFrags 1 -o salmon_output 
salmon quant -l SF -t combined.CDS.fasta -a WT_directRNA_mapped.combined.CDS_output.bam --noErrorModel --minAssignedFrags 1 -o salmon_output 
```

**Map RNAseq reads to D. ananassae and wAna for manual validation** 
```
#D. ananassae
hisat2-build Dana.UMIGS.fasta Dana.UMIGS.hisat2 
#single end
hisat2 -p 8 --max-intronlen 300000 -x Dana.UMIGS.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  
#paired end
hisat2 -p 8 --max-intronlen 300000 -x Dana.UMIGS.hisat2 -1 reads.R1.fastq.gz -2 reads.R2.fastq.gz | samtools view -bho output.bam -
#ont long reads
minimap2 -ax splice -uf -k14 -G 300000 Dana.UMIGS.fasta ont.fastq.gz | samtools view -bho output.bam -
#wAna
hisat2-build GCA_008033215.1_ASM803321v1_genomic.fna GCA_008033215.1_ASM803321v1_genomic.hisat2 
#single end
hisat2 -p 8 --max-intronlen 5000 -x GCA_008033215.1_ASM803321v1_genomic.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  
#paired end
hisat2 -p 8 --max-intronlen 5000 -x GCA_008033215.1_ASM803321v1_genomic.hisat2 -1 reads.R1.fastq.gz -2 reads.R2.fastq.gz | samtools view -bho output.bam -
#ont long reads
minimap2 -ax splice -uf -k14 -G 5000 GCA_008033215.1_ASM803321v1_genomic.fna ont.fastq.gz | samtools view -bho output.bam -
```

**Circos**
```
#gather coordinates of CDS on both strands
gffread -C /local/projects-t3/LGT/Dananassae_2020/wAna/GCA_008033215.1_ASM803321v1_genomic.gff --bed > wAna.CDS.bed
awk '{print $1"\t"$4"\t"$2"\t"$3"\t"$6}' wAna.CDS.bed | sed 's/gene-//g' - > wAna.CDS.slim.bed

awk '$6=="+"' wAna.CDS.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' > wAna.CDS.plus.coords.txt
awk '$6=="-"' wAna.CDS.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' > wAna.CDS.minus.coords.txt
#gather coordinates of transcribed genes
awk '$5>0' /local/projects-t3/LGT/Dananassae_2020/dana.tx/Dana-based/salmon_Dana.UMIGS_output/Illumina_cHI/wana.quant.sf | tail -n +2 - | awk '{print $1}' > Illumina.cHI.transcribed.IDs.list
awk '{print $1"\t"$2"\t"$3"\t"$4}' wAna.Illumina.cHI.transcribed.CDS.bed > wAna.Illumina.cHI.transcribed.CDS.coords.txt

/usr/local/packages/circos/bin/circos -conf circos_wAna.nuwt.conf
```
```
#GC skew
awk -v OFS='\t' {'print $1, $2'} GCA_008033215.1_ASM803321v1_genomic.fna.fai > ../wAna.circos/gc.skew/GCA_008033215.1_ASM803321v1_genomic.genomebed.bed
bedtools makewindows -g GCA_008033215.1_ASM803321v1_genomic.genomebed.bed -w 1000 > GCA_008033215.1_ASM803321v1_genomic.1kbp.windows.bed
bedtools nuc -fi ../../wAna/GCA_008033215.1_ASM803321v1_genomic.fna -bed GCA_008033215.1_ASM803321v1_genomic.1kbp.windows.bed > GCA_008033215.1_ASM803321v1_genomic.1kbp.nuc.stats.txt
tail -n +2 GCA_008033215.1_ASM803321v1_genomic.1kbp.nuc.stats.txt | awk '{print $1"\t"$2"\t"$3"\t"$8-$7"\t"$8+$7}' > GCA_008033215.1_ASM803321v1_genomic.1kbp.GC.stats.txt
sed -i 's/CP042904.1/wa1/g' GCA_008033215.1_ASM803321v1_genomic.1kbp.GC.stats.txt
awk '{print $1"\t"$2"\t"$3"\t"$4/$5}' GCA_008033215.1_ASM803321v1_genomic.1kbp.GC.stats.txt > GCA_008033215.1_ASM803321v1_genomic.1kbp.GC.skew.txt
```
