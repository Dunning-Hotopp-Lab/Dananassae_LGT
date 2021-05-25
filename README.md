# Extreme lateral gene transfer into a fly autosome

Eric S. Tvedte

2020-10-15



## Table of Contents
1. [Download D. ananassae and wAna datasets](#dl)
2. [Nuwt analysis](#nuwt)
3. [Numt analysis](#numt)
4. [LTR retrotransposon analysis](#ltr)  
5. [Transcription of LGT regions](#lgt.tx) 
11. [Data visualization] (#viz) data_viz_scripts 

### 1. Download D. ananassae and wAna datasets <a name="dl"></a>
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_genomic.fna.gz
gunzip *.gz  
D. ananassae PacBio Sequel II assembly available on Figshare https://doi.org/10.25387/g3.14096897
```

### 2. Nuwt analysis <a name="nuwt"></a>

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
bedtools coverage -a nuwt.finalpass.fixed.bed -b nuwt.finalpass.fixed.bed -hist | grep 'all' > nuwt.finalpass.fixed.hist.out 
#estimated nuwt length is 1*1 depth + (1/2)*2 depth in hist.out file, this corrects for small overlapping segments generated using NUCmer
```

**Estimate nuwt copy number and HiFi sequencing depth**
```
tail -n +5 nuwt.finalpass.coords | awk '{print $8"\t"$1-1"\t"$2}' > wAna.finalpass.bed
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

### 3. NUMT analysis <a name="numt"></a>
**Search for numts by alignments between mito genome and nuclear genome
```
nucmer -l 100 --maxmatch --prefix dana.numt.firstpass dana.mito.rotate.FREEZE.fasta /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/purge_haplotigs/dana.assembly.FREEZE.plusMITO.08.03.20.fasta

show-coords -r dana.numt.firstpass.delta > dana.numt.firstpass.coords  

tail -n +6 dana.numt.firstpass.coords | awk '{print $13}' | sort -n | uniq > dana.numt.contigs.list

samtools faidx /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/purge_haplotigs/dana.assembly.FREEZE.plusMITO.08.03.20.fasta -i contig_239:5130000-5170000 > dana.numt.contig_239.region.fasta

nucmer -l 100 --maxmatch --prefix dana.numt.finalpass dana.mito.rotate.FREEZE.fasta dana.numt.contig_239.region.fasta

```

### 4. LTR retrotransposon analysis <a name="ltr"></a>
**Create sequence database**  
```
genometools-1.5.9/bin/gt suffixerator -db Dana.UMIGS.fasta -indexname Dana.UMIGS.fasta -tis -suf -lcp -des -ssp -sds -dna
```

**Retrieve full-length LTR retrotransposons**
```
genometools-1.5.9/bin/gt ltrharvest -index Dana.UMIGS.fasta -seqids yes -tabout no -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -xdrop 7 -overlaps best > ltrharvest.bestovl.out
cp ltrharvest.bestovl.out ltrharvest.bestovl.gff
```

**Use family-based HMMs to retrieve BEL/PAO and Gypsy LTR retrotransposons**
```
genometools-1.5.9/bin/gt gff3 -sort ltrharvest.bestovl.gff > ltrharvest.bestovl.sorted.gff
#download BEL/PAO and Ty3/Gypsy HMMs from Gypsy database https://gydb.org/index.php?title=Main_Page
genometools-1.5.9/bin/gt ltrdigest -hmms *PAO.hmm -pdomevalcutoff 1E-20 -outfileprefix BELPAO.mapped.1E-20 -seqfile Dana.UMIGS.fasta -matchdescstart < ltrharvest.bestovl.sorted.gff > BELPAO.mapped.1E-20.gff
genometools-1.5.9/bin/gt ltrdigest -hmms *Gypsy.hmm -pdomevalcutoff 1E-20 -outfileprefix Ty3Gypsy.mapped.1E-20 -seqfile Dana.UMIGS.fasta -matchdescstart < ltrharvest.bestovl.sorted.gff > Ty3Gypsy.mapped.1E-20.gff
genometools-1.5.9/bin/gt select -rule_files filter_protein_match.rule < BELPAO.mapped.1E-20.gff > BELPAO.matches.gff #rule file in example_data folder
genometools-1.5.9/bin/gt select -rule_files filter_protein_match.rule < Ty3Gypsy.mapped.1E-20.gff > Ty3Gypsy.matches.gff #rule file in example_data folder
```

**Extract full-length retrotransposons and 5' + 3' LTR regions**
```
genometools-1.5.9/bin/gt extractfeat -type LTR_retrotransposon -seqfile Dana.UMIGS.fasta -matchdesc -coords -seqid BELPAO.matches.gff > BELPAO.retrotransposons.fasta
genometools-1.5.9/bin/gt extractfeat -type long_terminal_repeat -seqfile Dana.UMIGS.fasta -matchdesc -coords -seqid BELPAO.matches.gff > BELPAO.ltrs.fasta
genometools-1.5.9/bin/gt extractfeat -type LTR_retrotransposon -seqfile Dana.UMIGS.fasta -matchdesc -coords -seqid Ty3Gypsy.matches.gff > Ty3Gypsy.retrotransposons.fasta
genometools-1.5.9/bin/gt extractfeat -type long_terminal_repeat -seqfile Dana.UMIGS.fasta -matchdesc -coords -seqid Ty3Gypsy.matches.gff > Ty3Gypsy.ltrs.fasta
```

**Following manual inspection, remove duplicated records when retrotransposons score better than 1E-20 for both BEL/PAO and Ty3/Gypsy**
```
seqkit grep -v -f BELPAO.ltr.remove.list BELPAO.ltrs.fasta > BELPAO.ltrs.final.fasta
seqkit grep -v -f Ty3Gyspy.ltr.remove.list Ty3Gypsy.ltrs.fasta > Ty3Gypsy.ltrs.final.fasta
```
**Calculate pairwise distances:BEL/PAO**
```
bbtools/reformat.sh in=BELPAO.ltrs.final.fasta out1=BELPAO.5ltr.fasta out2=BELPAO.3ltr.fasta
awk '/^>/{print ">BELPAO_ltr_5_" ++i; next}{print}' < BELPAO.5ltr.fasta > BELPAO.5ltr.rn.fasta
awk '/^>/{print ">BELPAO_ltr_5_" ++i; next}{print}' < BELPAO.3ltr.fasta > BELPAO.3ltr.rn.fasta
seqkit split -i -O BELPAO.5.split BELPAO.5ltr.rn.fasta
seqkit split -i -O BELPAO.3.split BELPAO.3ltr.rn.fasta
for i in $(seq 1 942); do needle -asequence BELPAO.5.split/*5_$i.fasta -bsequence BELPAO.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile BELPAO.$i.align.fasta; done
for f in *fasta; do /usr/local/packages/gblocks/Gblocks $f -t=d -e=-g.fa -b1=2 -b3=1 -p=n ; done
rename fasta-g.fa gblocks.fasta *fa
for f in *gblocks.fasta; do sed -i 's/ //g' $f; done
for i in $(seq 1 942); do distmat -sequence BELPAO.$i.align.fasta -nucmethod 2 -outfile BELPAO.$i.k2p.distmat.out; done
for f in *distmat.out; do awk '{print $2, $3}' $f | tail -n 2 | head -n 1 >> UMIGS.all.BELPAO.k2p.distmat.out; done
sort -n -k2.15 UMIGS.all.BELPAO.k2p.distmat.out > UMIGS.all.BELPAO.k2p.distmat.sorted.out
paste <(grep '>' BELPAO.5ltr.fasta) <(grep '>' BELPAO.3ltr.fasta) <(grep '>' BELPAO.5ltr.rn.fasta) <(awk '{print $1}' UMIGS.all.BELPAO.k2p.distmat.sorted.out) | column -t | sed 's/>//g' > UMIGS.all.BELPAO.k2p.distmat.final.out
```

**Calculate pairwise distances:Ty3/Gypsy**
```
bbtools/reformat.sh in=Ty3Gypsy.ltrs.final.fasta out1=Ty3Gypsy.5ltr.fasta out2=Ty3Gypsy.3ltr.fasta
awk '/^>/{print ">Ty3Gypsy_ltr_5_" ++i; next}{print}' < Ty3Gypsy.5ltr.fasta > Ty3Gypsy.5ltr.rn.fasta
awk '/^>/{print ">Ty3Gypsy_ltr_5_" ++i; next}{print}' < Ty3Gypsy.3ltr.fasta > Ty3Gypsy.3ltr.rn.fasta
seqkit split -i -O Ty3Gypsy.5.split Ty3Gypsy.5ltr.rn.fasta
seqkit split -i -O Ty3Gypsy.3.split Ty3Gypsy.3ltr.rn.fasta
for i in $(seq 1 1172); do needle -asequence Ty3Gypsy.5.split/*5_$i.fasta -bsequence Ty3Gypsy.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile Ty3Gypsy.$i.align.fasta; done
for f in *fasta; do /usr/local/packages/gblocks/Gblocks $f -t=d -e=-g.fa -b1=2 -b3=1 -p=n ; done
rename fasta-g.fa gblocks.fasta *fa
for f in *gblocks.fasta; do sed -i 's/ //g' $f; done
for i in $(seq 1 1172); do distmat -sequence Ty3Gypsy.$i.align.fasta -nucmethod 2 -outfile Ty3Gypsy.$i.k2p.distmat.out; done
for f in *distmat.out; do awk '{print $2, $3}' $f | tail -n 2 | head -n 1 >> UMIGS.all.Ty3Gypsy.k2p.distmat.out; done
sort -n -k2.17 UMIGS.all.Ty3Gypsy.k2p.distmat.out > UMIGS.all.Ty3Gypsy.k2p.distmat.sorted.out
paste <(grep '>' Ty3Gypsy.5ltr.fasta) <(grep '>' Ty3Gypsy.3ltr.fasta) <(grep '>' Ty3Gypsy.5ltr.rn.fasta) <(awk '{print $1}' UMIGS.all.Ty3Gypsy.k2p.distmat.sorted.out) | column -t | sed 's/>//g' > UMIGS.all.Ty3Gypsy.k2p.distmat.final.out
```

### NUMT <a name="dana.numt"></a>
**De novo assembly of mitochondrial genome using Illumina data** 
```perl
perl NOVOPlasty3.7.pl -c dana.mito.config.txt
smrttools/smrtcmds/bin/pbmm2 align Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.fasta hifi.ccs.bam dana.mito.hifi.arrow1_sorted.bam --sort -j 16 -J 8
smrttools/smrtcmds/bin/arrow -j 16 dana.mito.hifi.arrow1_sorted.bam -r Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.fasta -o dana.mito.hifi_variants.gff -o dana.mito.hifi.polished.fasta
```
**Identify numt contigs**
```
mummer/nucmer -l 200 --prefix numt.firstpass dana.mito.hifi.polished.fasta Dana.UMIGS.fasta
mummer/show-coords -rT numt.firstpass.delta > numt.firstpass.coords
tail -n +5 numt.firstpass.coords | awk '{print $9}' | sort -n | uniq > numt.contigs.list
seqkit grep Dana.UMIGS.fasta -f numt.contigs.list > Dana.UMIGS.numt.contigs.fasta
```
**Estimate total numt content
```
tail -n +5 numt.firstpass.coords | awk '{print $9"\t"$3"\t"$4}' > numt.coords.bed
Rscript fixbed.R numt.coords.bed numt.coords.fixed.bed
bedtools coverage -a numt.coords.fixed.bed -b numt.coords.fixed.bed -hist | grep 'all' > numt.coords.fixed.hist.out 
#estimated numt length is 1*1 depth + (1/2)*2 depth in hist.out file, this corrects for small overlapping segments generated using NUCmer
```
**Aligment of mitochondria to numt contigs**
```
~jdhotopp/bin/residues.pl Dana.UMIGS.numt.contigs.fasta > Dana.UMIGS.numt.contigs.residues
nucmer -l 200 --prefix nuwt.finalpass GCA_008033215.1_ASM803321v1_genomic.fna Dana.UMIGS.numt.contigs.fasta
mummer/delta-filter -q numt.finalpass.delta > numt.finalpass.filter
data_viz_scripts/Dana.LGT.Rmd
```


### 5. RNAseq analysis <a name="lgt.tx"></a>
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
salmon quant -l A -t combined.CDS.fasta -a cHI_directRNA_mapped.combined.CDS_output.bam --noErrorModel --gcBias --minAssignedFrags 1 -o salmon_output 
salmon quant -l A -t combined.CDS.fasta -a WT_directRNA_mapped.combined.CDS_output.bam --noErrorModel --gcBias --minAssignedFrags 1 -o salmon_output 
```

**Map short RNA reads** 
```
hisat2-build polished.contigs.fasta polished.contigs.hisat2  
hisat2 -p 8 --max-intronlen 300000 -x polished.contigs.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  

hisat2-build /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2

for f in /local/projects-t3/RDBKO/sequencing/Dana_illumina_RNA_SRA/*.fastq; do echo "hisat2 -p 8 --max-intronlen 300000 -x /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2 -U $f | samtools view -bho ${f%_1*}_output.bam -" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N hisat2 -cwd; done

paired
for f in /local/projects-t3/LGT/Dananassae_2020/sequencing/Dana_Hawaii_RNASeq/*1.fastq; do echo "hisat2 -p 8 --max-intronlen 300000 -x /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2 -1 $f -2 ${f%_1*}_2.fastq | samtools view -bho ${f%_1*}_output.bam -" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N hisat2 -cwd; done

```

**Sort BAM**
```
java -jar picard.jar SortSam I=output.bam O=sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

for f in *output.bam; do echo "java -Xmx2g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$f O=${f%_o*}_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/etvedte/" | qsub -P jdhotopp-lab -l mem_free=2G -N SortSam -cwd; done  

```

**Nucmer alignment**
```
nucmer --prefix wAna.LGT.max.1000 -l 1000 --maxmatch wAna_v2.complete.pilon.fasta wAna.LGT.contigs.fasta

```

**Custom plot**
```
/local/projects-t3/RDBKO/scripts/Mchung.LGT.mummerplot.Rmd using nucmer delta output file and LGT contigs residues file
```

**Generate coordinates file for nucmer alignments**
```
show-coords -rl wAna.LGT.max.1000.delta > wAna.LGT.max.1000.r.coords

NUCMER parsing: grep tig00000335 wAna.LGT.max.1000.r.coords
REPEATMASKER parsing: grep tig00000335 dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta.out
```


**Retrieve putative transposable element sequences**
```
samtools faidx /local/projects-t3/LGT/Dananassae_2020/dana.postassemblyLGT/FREEZE tig00000335:416500-418000 > tig00000335_416500_418000.fasta
```

**Conduct BLAST and Dfam searches of putative transposable elements** 
```
https://blast.ncbi.nlm.nih.gov/Blast.cgi
https://dfam.org/home
```

**Salmon**
*Download files*
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_genomic.fna.gz 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_genomic.gff.gz
gunzip *gz
```
*Construct salmon index with wAna genome decoy*
```
grep "^>" GCA_008033215.1_ASM803321v1_genomic.fna | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat GCA_008033215.1_ASM803321v1_cds_from_genomic.fna GCA_008033215.1_ASM803321v1_genomic.fna > wAna.cds+genome.fasta 
salmon index -t wAna.cds+genome.fasta -d decoys.txt -i salmon_index -k 25 
```
*Map single end reads*
```
salmon quant -i salmon_index -l A -r read_1.fq --validateMappings --gcBias -o output_dir -p 8
```
*Map paired end reads*
```
salmon quant -i salmon_index -l A -1 read_1.fq -2 read_2.fq --validateMappings --gcBias -o output_dir -p 8
```


**Validating TE content in different D. ananassae strains**
fastq-dump --split-files SRA.accession

for f in *_1.fastq; do echo "java -Xmx10g -jar /usr/local/packages/trimmomatic/trimmomatic-0.38.jar PE -threads 8 $f ${f%_1*}_2.fastq ${f%_1*}_paired_1.fastq ${f%_1*}_unpaired_1.fastq ${f%_1*}_paired_2.fastq ${f%_1*}_unpaired_2.fastq ILLUMINACLIP:2:30:10:5 LEADING:3 TRAILING:3 MINLEN:70 SLIDINGWINDOW:4:15" | qsub -P jdhotopp-lab -l mem_free=10G -N trimmomatic -q threaded.q -pe thread 8 -cwd; done

for f in *paired_1.fastq; do /usr/local/packages/bbtools/reformat.sh in1=$f in2=${f%_1*}_2.fastq out=${f%_1*}_interleaved.fastq; done

use deviate

bwa index /local/projects-t3/LGT/Dananassae_2020/dana.repeats/DAn.repbase.fasta
single copy genes must be present in TE library

echo "deviaTE --threads 12 --input_fq SRR2127155_interleaved.fastq --read_type phred+33 --library /local/projects-t3/LGT/Dananassae_2020/dana.repeats/DAn.repbase.fasta --families ALL --single_copy_genes SPO11_1,SPO11_2" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 12 -V -N deviaTE -cwd
