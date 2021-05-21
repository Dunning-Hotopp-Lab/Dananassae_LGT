# Extreme lateral gene transfer into a fly autosome

Eric S. Tvedte

2020-10-15



## Table of Contents
1. [Download D. ananassae and wAna datasets] (#dl)
2. [Nuwt analysis](#nuwt)  
8. [Numt analysis](#numt)
9. [LTR retrotransposon analysis](#ltr)  
10. [Transcription of LGT regions](#lgt.tx) 
11. [Data visualization] (#viz) data_viz_scripts 

### 1. Download D. ananassae and wAna datasets <a name="dl"></a>
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/033/215/GCA_008033215.1_ASM803321v1/GCA_008033215.1_ASM803321v1_gff.gz
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
bedtools coverage -a nuwt.finalpass.fixed.bed -b nuwt.finalpass.fixed.bed -hist | grep 'all' > nuwt.finalpass.fixed.hist #estimated nuwt is 1*1 depth + (1/2)*2 depth, this corrects for small overlapping segments generated using NUCmer
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

**Annotating and dating LTR retrotransposons**





**Map long RNA reads**
```
echo -e "/usr/local/packages/minimap2-2.10/bin/minimap2 -ax splice -uf -k14 -t 8 -G 300000 /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.MinION.directRNA.fastq | samtools view -bho Dana_directRNA_output.bam -" | qsub -q threaded.q -pe thread 8 -P jdhotopp-lab -l mem_free=10G -N minimap2 -cwd 
```

**Sort BAM**
```
java -Xmx2g -jar picard.jar SortSam I=mapped.bam O=sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
```

**Build repeat database**
```
echo "/usr/local/packages/repeatmodeler-1.0.11/BuildDatabase -name Dana.repeats.database -engine wublast /local/projects-t3/RDBKO/dana.postassembly/arrow/dana.hybrid.80X.contigs.arrow.polished.fasta" | qsub -P jdhotopp-lab -l mem_free=2G -N RepBuild -cwd
```

**Classify repeat families**
```
echo "/usr/local/packages/repeatmodeler-1.0.11/RepeatModeler -database Dana.repeats.database" | qsub -P jdhotopp-lab -l mem_free=2G -N RepeatModeler -cwd
```

kseek parse TRF
```
perl /local/projects-t3/LGT/Dananassae_2020/scripts/k-seek/TRF_collapse.pl dana.assembly.FREEZE.plusMITO.6.22.20.fasta.2.7.7.80.10.50.500.dat > dana.assembly.FREEZE.plusMITO.6.22.20.TRF.collapse.fasta
/usr/local/packages/bbtools/readlength.sh in=dana.assembly.FREEZE.plusMITO.6.22.20.TRF.collapse.fasta bin=1 max=1000 out=dana.assembly.FREEZE.plusMITO.6.22.20.TRF.collapse.lengthstats.txt
grep -v '#' dana.assembly.FREEZE.plusMITO.6.22.20.TRF.collapse.lengthstats.txt > dana.assembly.FREEZE.plusMITO.6.22.20.TRF.collapse.lengthstats.final.txt
Rscript
```
fuzznuc: search for locations of enriched kmers, putative centromere repeats
```
fuzznuc -pattern @sample2.pa --complement -pmismatch 2  -sequence /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/dana.assembly.FREEZE.plusMITO.6.1.20.fasta -rformat2 gff -outfile sample2.gff
```

cent/tel RT annotation
blastn -query CentTel.RTs.fasta -db /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/purge_haplotigs/dana.assembly.FREEZE.plusMITO.6.22.20.fasta -perc_identity 75 -evalue 1e-20 -outfmt '6 sseqid sstart send qseqid' -out dana.assembly.FREEZE.plusMITO.6.22.20.CentTel.RTs.blast.out


**Soft-masking genome using repeat families**
```
/usr/local/packages/repeatmasker-4.0.7/RepeatMasker -lib rugiaPahangiNuclearGenome.fa-families.fa /local/aberdeen2rw/julie/JM_dir/PahangiPilonFASTA/Repeats/BrugiaPahangiNuclearGenome.fa

echo "/usr/local/packages/repeatmasker-4.0.7/RepeatMasker -xsmall -pa 24 -engine wublast -a -lib Dana.repeats.database-families.fa dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 24 -N RepeatMasker -cwd

echo "/usr/local/packages/repeatmasker-4.0.7/RepeatMasker -xsmall -nolow -s -no_is -norna -pa 24 -engine wublast -a -lib /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/Dmel.Dfam3.0.families.new.fa dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 24 -N RepeatMasker -cwd
```

**Genome annotation using BRAKER**
```perl 
braker.pl --species=DrosophilaAnanassae --softmasking --genome=contigs.softmasked.fasta --bam=shortRNAreads_sorted.bam, longRNAreads_sorted.bam --prot_seq=dmelanogaster.uniprot.proteins.fasta --etpmode --prg=gth --gth2traingenes --cores=8
```


### 6. Summary statistics for D. ananassae genome assembly <a name="sumstats"></a>
**QUAST**
```
quast.py dana.UMIGS.FREEZE.fasta Dana.pass.minimap2.racon.x3.pilon.x3.fasta -r GCA_000005115.1_dana_caf1_genomic.fasta --features gene:GCA_000005115.1_dana_caf1_genomic.gff -o quast_UMIGS.Miller2018.caf1 -t 24 --large -m 0 --fragmented --split-scaffolds --conserved-genes-finding
```

**BUSCO**
```python
python run_BUSCO.py -f -c 8 -t /tmp/folder -i dana.UMIGS.FREEZE.fasta -o dana.UMIGS.FREEZE -l arthropoda_odb9 -m geno
python run_BUSCO.py -f -c 8 -t /tmp/folder -i Dana.pass.minimap2.racon.x3.pilon.x3.fasta -o Miller2018 -l arthropoda_odb9 -m geno

```

**TRF**
```
trf dana.UMIGS.FREEZE.fasta 2 7 7 80 10 50 500 -h -m -f -d  
samtools faidx dana.UMIGS.FREEZE.fasta  
awk -v OFS='\t' {'print $1, $2'} dana.UMIGS.FREEZE.fasta.fai > dana.UMIGS.FREEZE.genomebed.bed  
bedtools makewindows -w 100000 -g dana.UMIGS.FREEZE.genomebed.bed > dana.UMIGS.FREEZE.100kbpwindows.bed
bedtools nuc -fi dana.UMIGS.FREEZE.trf.masked.fasta -bed dana.UMIGS.FREEZE.100kbpwindows.bed > dana.UMIGS.FREEZE.bednuc.out
grep -v '#' dana.assembly.FREEZE.plusMITO.08.03.20.bednuc.out |  awk '{print $1"\t"$2"\t"$3"\t"$10/$12}' > dana.assembly.FREEZE.plusMITO.08.03.20.bednuc.final.out
```

**BLAST for centromere/telomere retrotransposons**
```
makeblastdb dana.UMIGS.FREEZE.fasta -out dana.UMIGS.FREEZE.fasta -dbtype nucl -parse_seqids
blastn -query dana.cent.tel.RTs.fasta -db dana.UMIGS.FREEZE.fasta -outfmt ’6 qseqid staxids bitscore std’ -perc_identity 75 -evalue 1e-20 > dana.cent.tel.RTs.blast.out
```

### 7. Nuwt analysis <a name="nuwt"></a>

**Identifying LTR retrotransposons with LTRharvest and LTRdigest**  
*Create sequence database* 
```
/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt suffixerator -db /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -indexname /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -tis -suf -lcp -des -ssp -sds -dna
```

*Retrieve full-length LTR retrotransposons with LTRharvest*
```
/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt suffixerator -db chr4.contigs.fasta -indexname chr4.contigs.fasta -tis -suf -lcp -des -ssp -sds -dna     /local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt ltrharvest -index chr4.contigs.fasta -seqids yes -tabout no -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -xdrop 7 -overlaps best > chr4.contigs.ltrharvest.bestovl.out
```

*Use family-based HMMs to retrieve BEL/PAO and Gypsy LTR retrotransposons*
```
cp dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltrharvest.bestovl.out dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltrharvest.bestovl.gff

/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt gff3 -sort dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltrharvest.bestovl.gff > dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltrharvest.bestovl_sorted.gff

/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt ltrdigest -hmms ../*PAO.hmm -pdomevalcutoff 1E-20 -outfileprefix dana.FREEZE.bestovl.PAO.mapped.1E-20 -seqfile /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -matchdescstart < dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltrharvest.bestovl_sorted.gff > dana.FREEZE.bestovl.PAO.mapped.1E-20.gff

/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt select -rule_files ../filter_protein_match.rule < dana.FREEZE.bestovl.PAO.mapped.gff > dana.FREEZE.bestovl.PAO.matches.gff
```

*Extract full-length retrotransposons and 5' + 3' LTR regions*
```
/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt extractfeat -type LTR_retrotransposon -seqfile /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -matchdesc -coords -seqid dana.FREEZE.bestovl.PAO.matches.gff > dana.FREEZE.bestovl.PAO.retrotransposons.fasta

/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt extractfeat -type long_terminal_repeat -seqfile /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -matchdesc -coords -seqid dana.contigs.FREEZE.ltrdigest.mappedids.PAOhits.gff > dana.contigs.FREEZE.ltrdigest.mappedids.PAO.ltrs.fasta
```

*Following manual inspection, remove duplicated records when retrotransposons score better than 1E-20 for both BEL/PAO and Ty3/Gypsy*
```
seqkit grep -v -f BELPAO.ltr.remove.list UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.fasta > UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.fasta
seqkit grep -v -f Ty3Gyspy.ltr.remove.list UMIGS.FREEZE.2021.ltrharvest.bestovl.Ty3Gypsy.ltrs.fasta > UMIGS.FREEZE.2021.ltrharvest.bestovl.Ty3Gypsy.ltrs.final.fasta
```
*Calculate pairwise distances:Ty3/Gypsy*
```
/usr/local/packages/bbtools/reformat.sh in=dana.contigs.FREEZE.ltrdigest.mappedids.PAO.ltrs.fasta out1=dana.contigs.FREEZE.ltrdigest.mappedids.PAO.5ltr.fasta out2=dana.contigs.FREEZE.ltrdigest.mappedids.PAO.3ltr.fasta

rename sequentially 
awk '/^>/{print ">Gypsy_ltr_5_" ++i; next}{print}' < dana.FREEZE.final.Gypsy.5ltr.fasta > dana.FREEZE.final.Gypsy.5ltr.rn.fasta
awk '/^>/{print ">Gypsy_ltr_3_" ++i; next}{print}' < dana.FREEZE.final.Gypsy.3ltr.fasta > dana.FREEZE.final.Gypsy.3ltr.rn.fasta

trim sequences?
seqkit subseq -r 4:-4 dana.all.contigs.bestovl.BELPAO.5ltr.rn.fasta > dana.all.contigs.bestovl.BELPAO.5ltr.rn.trim.fasta

split into individual files
seqkit split -i -O Gypsy.5.split dana.FREEZE.final.Gypsy.5ltr.rn.fasta
seqkit split -i -O Gypsy.3.split dana.FREEZE.final.Gypsy.3ltr.rn.fasta
```

*Pairwise alignment of LTRs using Needle*
```
for i in $(seq 1 1548); do needle -asequence Gypsy.5.split/*5_$i.fasta -bsequence Gypsy.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile Gypsy.needle+distmat/dana.Gypsy.$i.align.fasta; done


ls 

for f in *fasta; do /usr/local/packages/gblocks/Gblocks $f -t=d -e=-g.fa -b1=2 -b3=1 -p=n ; done
rename fasta-g.fa gblocks.fasta *fa
for f in *gblocks.fasta; do sed -i 's/ //g' $f; done
```

*Measure pairwise distance values using emboss dismat*
```
for i in $(seq 1 1548); do distmat -sequence Gypsy.needle+distmat/dana.Gypsy.$i.align.fasta -nucmethod 0 -outfile Gypsy.needle+distmat/dana.Gypsy.$i.uncorr.distmat.out; done

K2P
for i in $(seq 1 1548); do distmat -sequence Gypsy.needle+distmat/dana.Gypsy.$i.align.fasta -nucmethod 2 -outfile Gypsy.needle+distmat/dana.Gypsy.$i.k2p.distmat.out; done

#cat PAO.needle+distmat/*uncorr.out > PAO.needle+distmat/PAO.distmat.uncorr.all.out

for f in PAO.needle+distmat/*uncorr.distmat.out; do awk '{print $2, $3}' $f | tail -n 2 | head -n 1 >> PAO.all.uncorr.distmat.out; done



sort -n -k2.12 PAO.all.uncorr.distmat.out > PAO.all.uncorr.distmat.sorted.out

sort -n -k2.17 UMIGS.Ty3Gypsy.all.uncorr.distmat.out > UMIGS.Ty3Gypsy.all.uncorr.distmat.sorted.out

paste <(grep '>' dana.FREEZE.final.PAO.5ltr.fasta ) <(grep '>' dana.FREEZE.final.PAO.3ltr.fasta ) <(grep '>' dana.FREEZE.final.PAO.5ltr.rn.fasta) <(awk '{print $1}' PAO.needle+distmat/PAO.all.uncorr.distmat.sorted.out) <(awk '{print $1}' PAO.needle+distmat/PAO.all.k2p.distmat.sorted.out) | column -t -o $'\t' | sed 's|>||g' | sed 's| ||g' > PAO.needle+distmat/PAO.all.final.distmat.out

paste <(grep '>' dana.FREEZE.final.Gypsy.5ltr.fasta ) <(grep '>' dana.FREEZE.final.Gypsy.3ltr.fasta ) <(grep '>' dana.FREEZE.final.Gypsy.5ltr.rn.fasta) <(awk '{print $1}' Gypsy.needle+distmat/Gypsy.all.uncorr.distmat.sorted.out) <(awk '{print $1}' Gypsy.needle+distmat/Gypsy.all.k2p.distmat.sorted.out) | column -t | sed 's/>//g' > Gypsy.needle+distmat/Gypsy.all.final.distmat.out

```
*Calculate pairwise distances:BEL/PAO*
```
#split interleaved LTR file
/usr/local/packages/bbtools/reformat.sh in=UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.fasta out1=UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.5ltr.fasta out2=UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.3ltr.fasta
#rename sequentially 
awk '/^>/{print ">BELPAO_ltr_5_" ++i; next}{print}' < UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.5ltr.fasta > UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.5ltr.rn.fasta
awk '/^>/{print ">BELPAO_ltr_3_" ++i; next}{print}' < UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.3ltr.fasta > UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.3ltr.rn.fasta
#split into individual files
seqkit split -i -O BELPAO.5.split UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.5ltr.rn.fasta
seqkit split -i -O BELPAO.3.split UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.3ltr.rn.fasta
#needle
for i in $(seq 1 942); do needle -asequence BELPAO.5.split/*5_$i.fasta -bsequence BELPAO.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile BELPAO.needle+distmat/UMIGS.BELPAO.$i.align.fasta; done
#gblocks
for f in *fasta; do /usr/local/packages/gblocks/Gblocks $f -t=d -e=-g.fa -b1=2 -b3=1 -p=n ; done
rename fasta-g.fa gblocks.fasta *fa
for f in *gblocks.fasta; do sed -i 's/ //g' $f; done
#distmat
for i in $(seq 1 942); do distmat -sequence BELPAO.needle+distmat/UMIGS.BELPAO.$i.align.fasta -nucmethod 2 -outfile BELPAO.needle+distmat/UMIGS.BELPAO.$i.k2p.distmat.out; done
#concatenate and sort output
for f in BELPAO.needle+distmat/*distmat.out; do awk '{print $2, $3}' $f | tail -n 2 | head -n 1 >> UMIGS.BELPAO.all.k2p.distmat.out; done
sort -n -k2.15 UMIGS.BELPAO.all.k2p.distmat.out > UMIGS.BELPAO.all.k2p.distmat.sorted.out
#build output table
paste <(grep '>' UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.5ltr.fasta ) <(grep '>' UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.3ltr.fasta  ) <(grep '>' UMIGS.FREEZE.2021.ltrharvest.bestovl.BELPAO.ltrs.final.5ltr.rn.fasta) <(awk '{print $1}' BELPAO.needle+distmat/UMIGS.BELPAO.all.k2p.distmat.sorted.out) | column -t | sed 's/>//g' > BELPAO.needle+distmat/UMIGS.BELPAO.all.k2p.distmat.final.out
```

### Transcription of LGT regions <a name="dana.lgt.tx"></a>
**Map Wolbachia-cured datasets to wAna using salmon**
```
grep "^>" wAna.genome.fasta | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat augustus.hints.codingseq ../dana.qm.merged.FREEZE.fasta > transcripts+genome.fasta 
/local/projects-t3/LGT/Dananassae_2020/scripts/salmon-latest_linux_x86_64/bin/salmon index -t transcripts+genome.fasta -d decoys.txt -i salmon_index -k 25
for f in /local/projects-t3/LGT/Dananassae_2020/sequencing/Dana_UMIGS_Hawaii_RNASeq/*dedup_R1.fastq.gz; do echo "/home/etvedte/scripts/salmon-latest_linux_x86_64/bin/salmon quant -i salmon_index -l ISR -1 $f -2 ${f%_R*}_R2.fastq.gz --validateMappings --gcBias --minAssignedFrags 1 -o /local/projects-t3/LGT/Dananassae_2020/dana.tx/wAna-based/salmon_output/cHI_UMIGS_stranded -p 8 --writeMappings | samtools view -bhSo /local/projects-t3/LGT/Dananassae_2020/dana.tx/wAna-based/salmon_output/cHI_UMIGS_dedup_stranded_output.bam" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 8 -cwd -N salmon_quant_cHI; done

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
### NUMT <a name="dana.numt"></a>
**De novo assembly of mitochondrial genome using Illumina data** 
```perl
echo "perl /local/projects-t3/RDBKO/scripts/NOVOPlasty/NOVOPlasty3.7.pl -c /local/projects-t3/RDBKO/scripts/NOVOPlasty/dana.mito.config.txt" | qsub -P jdhotopp-lab -l mem_free=20G -N NovoPlasty -cwd

```

**Polish with PacBio HiFi**
```
echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/pbmm2 align Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.fasta /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA_HiFi/cHI_Dana_2_15_19/PACBIO_DATA/RANDD_20191011_S64018_PL100122512-3_A01.ccs.bam Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.hifi.arrow1_sorted.bam --sort -j 16 -J 8" | qsub -P jdhotopp-lab -l mem_free=50G -N pbmm2.align -q threaded.q -pe thread 16 -cwd -V

echo /usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/arrow -j 16 Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.hifi.arrow1_sorted.bam -r Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.fasta -o Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.hifi.arrow1_variants.gff -o Circularized_assembly_1_dana.illumina.mito.arrow1.polished.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -N smrttools.arrow -q threaded.q -pe thread 16 -cwd

--noEvidenceConsensusCall {nocall,reference,lowercasereference}
                        The consensus base that will be output for sites with
                        no effective coverage. (default: lowercasereference)


```

##NUMT analysis
**Search for numts by alignments between mito genome and nuclear genome
```
nucmer -l 100 --maxmatch --prefix dana.numt.firstpass dana.mito.rotate.FREEZE.fasta /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/purge_haplotigs/dana.assembly.FREEZE.plusMITO.08.03.20.fasta

show-coords -r dana.numt.firstpass.delta > dana.numt.firstpass.coords  

tail -n +6 dana.numt.firstpass.coords | awk '{print $13}' | sort -n | uniq > dana.numt.contigs.list

samtools faidx /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/purge_haplotigs/dana.assembly.FREEZE.plusMITO.08.03.20.fasta -i contig_239:5130000-5170000 > dana.numt.contig_239.region.fasta

nucmer -l 100 --maxmatch --prefix dana.numt.finalpass dana.mito.rotate.FREEZE.fasta dana.numt.contig_239.region.fasta

```

**Validating TE content in different D. ananassae strains**
fastq-dump --split-files SRA.accession

for f in *_1.fastq; do echo "java -Xmx10g -jar /usr/local/packages/trimmomatic/trimmomatic-0.38.jar PE -threads 8 $f ${f%_1*}_2.fastq ${f%_1*}_paired_1.fastq ${f%_1*}_unpaired_1.fastq ${f%_1*}_paired_2.fastq ${f%_1*}_unpaired_2.fastq ILLUMINACLIP:2:30:10:5 LEADING:3 TRAILING:3 MINLEN:70 SLIDINGWINDOW:4:15" | qsub -P jdhotopp-lab -l mem_free=10G -N trimmomatic -q threaded.q -pe thread 8 -cwd; done

for f in *paired_1.fastq; do /usr/local/packages/bbtools/reformat.sh in1=$f in2=${f%_1*}_2.fastq out=${f%_1*}_interleaved.fastq; done

use deviate

bwa index /local/projects-t3/LGT/Dananassae_2020/dana.repeats/DAn.repbase.fasta
single copy genes must be present in TE library

echo "deviaTE --threads 12 --input_fq SRR2127155_interleaved.fastq --read_type phred+33 --library /local/projects-t3/LGT/Dananassae_2020/dana.repeats/DAn.repbase.fasta --families ALL --single_copy_genes SPO11_1,SPO11_2" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 12 -V -N deviaTE -cwd
