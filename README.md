# Extreme lateral gene transfer into a fly autosome

Eric S. Tvedte

2020-1-9

The repository contains Supplementary Data for the manuscript, including Tables, Figures, and Files.

## Table of Contents
1. [Pre-processing MinION sequencing data](#seq.prep)
2. [Genome assembly](#assemble)
3. [Post-assembly processing](#post)
4. [Genome annotation](#annotate)
5. [BUSCO analysis](#busco)
6. [Characterization of euchromatic regions of D. ananassae](#dana.chrom.map)  
7. [Characterization of Y contigs in D. ananassae](#dana.y)  
8. [Characterization of chromosome 4 contigs in D. ananassae](#dana.chr4)
9. [Characterization of LGT contigs in D. ananassae](#dana.lgt)  
10. [Transcription of LGT regions](#dana.lgt.tx)  
11. [NUMT](#dana.numt)


### Pre-processing MinION sequencing data <a name="seq.prep"></a>  
**Basecalling with guppy**
```
/usr/local/packages/guppy-3.1.5/bin/guppy_basecaller --input_path fast5_dir --save_path output_dir --config dna_r9.4.1_450bps_fast.cfg --fast5_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 --num_callers 8 --cpu_threads_per_caller 4  
```

**Filter lambda from sequencing data**
```
zcat minion.LIG.raw.fastq.gz | NanoLyse --reference lambda.DCS.fasta | gzip > minion.LIG.filter.lambda.fastq.gz
```

### Genome assembly <a name="assemble"></a>  
**Canu**  
```
canu -p output.prefix -d output.dir genomeSize=240m corOutCoverage=80 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq.gz -nanopore-raw minion.LIG.filter.lambda.fastq.gz
```

**Flye** 
```
echo -e "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/gcc/lib64\nexport PYTHONPATH=$PYTHONPATH:/usr/local/packages/flye-2.4.2/lib/python2.7/site-packages\n/usr/local/packages/flye-2.4.2/bin/flye -g 240m -t 24 -o /local/projects-t3/RDBKO/dana.flye/ --asm-coverage 60 --pacbio-raw /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.pbSequelII.raw.fastq.gz
```

### Post-assembly processing <a name="post"></a>
**Map PacBio Sequel II data**  
```
use python-3.5
echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/pbmm2 align dana.hybrid.contigs.fasta /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA/RANDD_20190301_S64018_PL100122512-1_C01.subreads.bam dana.hybrid.80X.contigs.mapped.pb.sqII_sorted.bam --sort -j 16 -J 8" | qsub -P jdhotopp-lab -l mem_free=50G -N pbmm2.align -q threaded.q -pe thread 16 -cwd -V
```

**polish with arrow using Sequel II data**
```
echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/pbmm2 align Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.fasta /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA/RANDD_20190301_S64018_PL100122512-1_C01.subreads.bam Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.sqII.arrow1_sorted.bam --sort -j 16 -J 8" | qsub -P jdhotopp-lab -l mem_free=50G -N pbmm2.align -q threaded.q -pe thread 16 -cwd -V
```

**Map PacBio HiFi data**  


**polish with FreeBayes**

**purge haplotigs from assembly**
```
echo "minimap2 -xmap-pb /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.fasta /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.pbSequelII.raw.fastq.gz | gzip -c - > dana.hybrid.80X.arrow.rd2.mappedsqII.paf.gz" | qsub -P jdhotopp-lab -l mem_free=10G -N minimap2 -cwd

/local/projects-t3/RDBKO/scripts/purge_dups/bin/split_fa /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.fasta > /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.split

/home/etvedte/scripts/purge_dups/bin/pbcstat dana.hybrid.80X.contigs.arrow.polished.mappedhifi.paf.gz
/home/etvedte/scripts/purge_dups/bin/calcuts PB.stat > cutoffs 2> calcuts.log
/home/etvedte/scripts/purge_dups/scripts/hist_plot.py PB.stat hist.out.pdf
```

**Convert bases to upper case**
*By default arrow outputs regions with no consensus as lower case, i.e. 'acgt'. In order to properly annotate repetitive regions as lower case, all bases must be converted to upper case. Note that this could have been accomplished in arrow using the parameter --noEvidenceConsensusCall reference* 
```
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta > dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fmt.fasta
```

**heterochromatin-sensitive assembly**
*Map long reads to genome*
```
echo "minimap2 -ax map-ont -t 16 /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta /local/projects-t3/RDBKO/sequencing/RANDD_LIG_Dana_20190405_merged_pass_filtlambda.fastq | samtools view -bho dana.hybrid.80X.arrow.rd2.contigs.FREEZE.mapped.LIG_output.bam" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -cwd -N minimap2
```

*sort Sam*
```
echo -e "java -Xmx20g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=dana.hybrid.80X.arrow.rd2.contigs.FREEZE.mapped.LIG_output.bam O=dana.hybrid.80X.arrow.rd2.contigs.FREEZE.mapped.LIG_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/etvedte" | qsub -P jdhotopp-lab -l mem_free=20G -N samsort -cwd
```

*Generate BED file for non-chromosomal contigs*
```
awk -v OFS='\t' {'print $1, $2'} ../dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta.fai > dana.hybrid.80X.arrow.rd2.contigs.FREEZE.genomebed.bed

bedtools makewindows -g dana.hybrid.80X.arrow.rd2.contigs.FREEZE.genomebed.bed -w 10000000000 > dana.hybrid.80X.arrow.rd2.contigs.FREEZE.bed
```
*Retrieve reads mapping to unlocalized contigs and unmapped*
```
echo "samtools view -@16 -L dana.hybrid.80X.arrow.rd2.contigs.FREEZE.bed -b /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.mapped.sqII_sorted.bam | samtools fasta - > dana.pb.heterochromatin.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -cwd -N sam.fasta
echo "samtools view -@16 -f 4 /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.mapped.sqII_sorted.bam | samtools fasta - > dana.pb.sqII.unmapped.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -cwd -N sam.fasta

```
**or BEDtools**
```
echo "samtools view -@16 -L nonchr.contig.list -b /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.mapped.sqII_sorted.bam | bedtools bamtofastq -i - -fq dana.pb.sqII.heterochromatin.fastq" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -cwd -N bam.to.fastq
```

**assembly**
```
/usr/local/packages/canu-1.8/bin/canu -p dana.het -d /local/projects-t3/LGT/Dananassae_2020/dana.het.assembly/pb+minion.fasta genomeSize=85m corOutCoverage=100 corMinCoverage=0 stopOnReadQuality=false MhapSensitivity=normal gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw /local/projects-t3/LGT/Dananassae_2020/dana.het.assembly/dana.pb.sqII.het.unmap.fasta -nanopore-raw /local/projects-t3/LGT/Dananassae_2020/dana.het.assembly/dana.minion.LIG.het.unmapped.fasta
```

**NUCMER**
```
nucmer --maxmatch --prefix FREEZE.chrs -l 200 dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta
```

**QUAST**
```
use python-3.5
echo "/home/etvedte/scripts/quast-5.0.2/quast.py /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/purge_haplotigs/dana.assembly.FREEZE.plusMITO.6.22.20.fasta /local/projects-t3/RDBKO/nonIGS_dana/Miller2018/Dana.pass.minimap2.racon.x3.pilon.x3.fasta /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.redux.fasta -r /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.redux.fasta --features gene:/local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.gff -o quast_FREEZE.6.22.20.Miller.caf1 -t 24 --large -m 0 --fragmented --split-scaffolds --conserved-genes-finding" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 24 -N quast.LG -cwd -V

```

**KAT**
```
use kat-2.4.0
echo "kat comp -t 16 -o FREEZE_vs_illuminaR1 /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R1.fastq /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -N kat.comp.illumina -cwd -V

```

### Genome annotation <a name="annotate"></a>  
**Map short RNA reads** 
```
hisat2-build polished.contigs.fasta polished.contigs.hisat2  
hisat2 -p 8 --max-intronlen 300000 -x polished.contigs.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  


for f in /local/projects-t3/RDBKO/sequencing/Dana_illumina_RNA_SRA/*.fastq; do echo "hisat2 -p 8 --max-intronlen 300000 -x /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2 -U $f | samtools view -bho /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/$(basename ${f%_1*})_output.bam -" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N hisat2 -cwd; done

```

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

**TRF**
```
echo "trf dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta 2 7 7 80 10 50 500 -h -m -f -d" | qsub -P jdhotopp-lab -N trf -l mem_free=20G -cwd

samtools faidx dana.qm.merged.FREEZE.fasta
awk -v OFS='\t' {'print $1, $2'} dana.qm.merged.FREEZE.fasta.fai > dana.qm.merged.FREEZE.genomebed.bed
bedtools makewindows -w 100000 -g dana.qm.merged.FREEZE.genomebed.bed > dana.qm.merged.FREEZE.100kbpwindows.bed
bedtools nuc -fi trf.masked.fasta -bed dana.qm.merged.FREEZE.100kbpwindows.bed > dana.qm.merged.bednuc.out

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

echo "/usr/local/packages/repeatmasker-4.0.7/RepeatMasker -xsmall -nolow -s -no_is -norna -gccalc -pa 24 -engine wublast -a -lib /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/Dmel.Dfam3.0.families.new.fa dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 24 -N RepeatMasker -cwd
```

**Genome annotation using BRAKER**
```perl 
braker.pl --species=DrosophilaAnanassae --softmasking --genome=contigs.softmasked.fasta --bam=shortRNAreads_sorted.bam, longRNAreads_sorted.bam --prot_seq=dmelanogaster.uniprot.proteins.fasta --etpmode --prg=gth --gth2traingenes --cores=8
```

### BUSCO analysis <a name="busco"></a>
```python
python run_BUSCO.py -f -c 8 -t /local/scratch/etvedte/tmp -i polished.contigs.fasta -o outdir_prefix -l metazoa_odb9 -m geno
```

### Characterization of euchromatic regions of D. ananassae <a name="dana.chrom.map"></a>
**Initial BLAST search to determine contig orientation**
```
makeblastdb -in polished.contigs.fasta -out polished.contigs.fasta -dbtype nucl -parse_seqids
blastn -query mapping_loci.fasta -db polished.contigs.fasta -max_target_seqs 10 -max_hsps 10 -outfmt 6 > initial.blast.out
```

**Extract matching contigs, reverse complement if necessary**
```
samtools faidx polished.contigs.fasta fwd.contig.name >> polished.contigs.correct.fasta
samtools faidx -i polished.contigs.fasta rev.contig.name >> polished.contigs.correct.fasta
```

**Final BLAST search, using matched contigs as database and custom BLAST output**
```
makeblastdb -in polished.contigs.correct.fasta -out polished.contigs.correct.fasta -dbtype nucl -parse_seqids
blastn -query mapping_loci.fasta -db polished.contigs.correct.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > final.blast.out
```

**Visualization of contigs, including locations of identified euchromatic loci**
```
contigs.Rmd using data.txt
```

### Characterization of Y contigs in D. ananassae <a name="dana.y"></a>  
**Index reference genome**
```
bwa index /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta
bwa mem -k 23 -t 8 female.R1.fq.gz female.R2.fq.gz
bwa mem -k 23 -t 8 male.R1.fq.gz male.R2.fq.gz

for f in /local/projects/JULIE/Dana*/ILLUMINA_DATA/*R1_trimmed.fastq.gz; do echo "bwa mem -k 23 -t 8 /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta $f ${f%_R*}_R2_trimmed.fastq.gz | samtools view -bho /local/projects-t3/RDBKO/dana.chrY/FREEZE/$(basename ${f%_R*})_output.bam" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N bwamem -cwd ; done
```

**Map female short reads**
``` 
echo -e "bwa mem -t 8 -k 23 /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa /local/projects/JULIE/Dana_cHi_female_2/ILLUMINA_DATA/JULIE_20190729_K00134_IL100134730_MX83_L001_R1_trimmed.fastq.gz /local/projects/JULIE/Dana_cHi_female_2/ILLUMINA_DATA/JULIE_20190729_K00134_IL100134730_MX83_L001_R2_trimmed.fastq.gz | samtools view -bho /local/projects-t3/RDBKO/dana.chrY/purged.contigs.mapped.cHi_female2_output.bam -" | qsub -P jdhotopp-lab -l mem_free=50G -N BWAMEM -q threaded.q -pe thread 8 -cwd
```

**Map male short reads**
``` 
echo -e "bwa mem -t 8 -k 23 /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa /local/projects/JULIE/Dana_cHi_male_2/ILLUMINA_DATA/JULIE_20190702_K00134_IL100134728_MX81_L008_R1_trimmed.fastq.gz /local/projects/JULIE/Dana_cHi_male_2/ILLUMINA_DATA/JULIE_20190702_K00134_IL100134728_MX81_L008_R2_trimmed.fastq.gz | samtools view -bho /local/projects-t3/RDBKO/dana.chrY/purged.contigs.mapped.cHi_male2_output.bam -" | qsub -P jdhotopp-lab -l mem_free=50G -N BWAMEM -q threaded.q -pe thread 8 -cwd
```

**Sort BAM and remove duplicates**
```
for f in *output.bam; do echo "java -Xmx2g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$f O=${f%_o*}_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/etvedte/" | qsub -P jdhotopp-lab -l mem_free=2G -N SortSam -cwd; done  
for f in *sorted.bam; do echo "java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$f O=${f%_s*}_dedup.bam M=${f%_s*}_dedup.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true" | qsub -P jdhotopp-lab -l mem_free=10G -N MarkDups -cwd; done
```

**Calculate sequencing depth, filter reads with MAPQ < 10**
```
samtools depth -Q 10 mel_f_new.bam mel_m_new.bam > mel_new.out
```

**Determine average and median female/male depth in 10kb windows (see Chang and Larracuente, 2018)** 
```perl
female
JULIE_20190729_K00134_IL100134730_MX83_L001
JULIE_20190729_K00134_IL100134731_MX84_L001
male
JULIE_20190702_K00134_IL100134728_MX81_L008
JULIE_20190702_K00134_IL100134729_MX82_L008

perl /local/projects-t3/RDBKO/scripts/Chang2019_frame_depth_new.pl mel_new.out
```

### Characterization of chromosome 4 contigs in D. ananassae <a name="dana.chr4"></a>
**Nucmer alignment to identify LGT contigs
```
nucmer --maxmatch --prefix chr4 -l 1000 ../Leung2017_chr4_scaffolds.fasta /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta
show-coords -r chr4.delta > chr4.coords
tail -n +6 chr4.coords | awk '{print $13}' | sort -n | uniq > chr4.contigs.list
xargs samtools faidx /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta < chr4.contigs.list >> chr4.contigs.fasta
```

**Nucmer aligment to chromosome 4 contigs**
```
nucmer --maxmatch --prefix chr4.match -l 1000 ../Leung2017_chr4_scaffolds.fasta chr4.contigs.fasta
mummerplot --prefix chr4.match --color --postscript chr4.match.delta
```

### Characterization of LGT contigs in D. ananassae <a name="dana.lgt"></a>
**Nucmer alignment to identify LGT contigs
```
nucmer --maxmatch --prefix wAna.LGT /local/aberdeen2rw/julie/Matt_dir/EWANA/references/wAna_v2.complete.pilon.fasta /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa
show-coords -r wAna.LGT.delta > wAna.LGT.r.coords
cat wAna.LGT.r.coords | tail -n +6 | awk '{print $13}' | sort -n | uniq > wAna.LGT.contigs.list
xargs samtools faidx /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa < wAna.LGT.contigs.list >> wAna.LGT.contigs.fasta
```

**Nucmer aligment to LGT contigs**
```
~jdhotopp/bin/residues.pl wAna.LGT.contigs.fasta > wAna.LGT.contigs.residues
nucmer --maxmatch --prefix wAna.LGT.only /local/aberdeen2rw/julie/Matt_dir/EWANA/references/wAna_v2.complete.pilon.fasta wAna.LGT.contigs.fasta
/local/projects-t3/RDBKO/scripts/Mchung.LGT.mummerplot.Rmd with residues.txt and LGT.match.delta
```

**calculate LGT segment length**
```
nucmer -l 5000 --prefix LGT.5kbp.segments wAna_v2.complete.pilon.fasta wAna.LGT.contigs.FREEZE.fasta
show-coords -rl LGT.5kbp.segments.delta > LGT.5kbp.segments.rl.coords
cat LGT.5kbp.segments.rl.coords | tail -n +6 | awk '{print $8}' | awk '{total = total + $1}END{print "Total LGT segment length = "total}' -
```

**Estimate wAna LGT segment copy number**
```
echo "bwa mem -k 23 -t 8 wAna_v2.complete.pilon.fasta /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R1.fastq /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R2.fastq | samtools view -bho wAna_v2.complete.pilon.mapped.illumina_output.bam" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N bwamem -cwd
for f in *output.bam; do echo "java -Xmx2g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$f O=${f%_o*}_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/etvedte/" | qsub -P jdhotopp-lab -l mem_free=2G -N SortSam -cwd; done  
for f in *sorted.bam; do echo "java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$f O=${f%_s*}_dedup.bam M=${f%_s*}_dedup.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true" | qsub -P jdhotopp-lab -l mem_free=10G -N MarkDups -cwd; done
samtools depth -aa -m 100000000 dedup.bam > dedup.depth.txt
```

**Read support for individual LGT regions (sliding window bed)**
```
for f in dana.assembly.FREEZE.plusMITO.6.22.20.mapped*sorted.bam; do echo "bamCoverage -p 4 -b $f -of bedgraph -o ${f%_s*}_sorted.depth.bed" | qsub -P jdhotopp-lab -N bamCov -l mem_free=10G -q threaded.q -pe thread 4 -cwd -V; done
```
**Read support for individual LGT regions (read visualization)**
```
echo "bedtools intersect -a dana.assembly.FREEZE.plusMITO.6.22.20.mapped.SQII_sorted.bam -b test.intervals.bed -F 1 > test.result.bam" | qsub -P jdhotopp-lab -l mem_free=10G -N bed.INT -cwd
```
**Mugsy**
```
source /local/projects/angiuoli/mugsy/mugsyenv.sh
/local/projects/angiuoli/mugsy/mugsy --directory /local/projects-t3/RDBKO/dana.LGT/80X.polished.rd1/mugsy --prefix ecoli 1.ecoli.fasta 2.ecoli.fasta 3.ecoli.fasta
source /home/jdhotopp/bin/jsahl_mugsy_to_tree_dir/pythonenv.sh
/home/jdhotopp/bin/jsahl_mugsy_to_tree_dir/process_maf.sh 432481_433504_CDS.maf
```

**Identifying LTR retrotransposons with LTRharvest and LTRdigest**  
*Create sequence database* 
```
/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt suffixerator -db /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -indexname /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -tis -suf -lcp -des -ssp -sds -dna
```

*Retrieve full-length LTR retrotransposons with LTRharvest*
```
echo "/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt ltrharvest -index /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta -seqids yes -tabout no -mindistltr 2000 -maxlenltr 2000 -overlaps best > dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltrharvest.bestovl.out" | qsub -P jdhotopp-lab -l mem_free=10G -cwd -N ltrharvest
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

*parse data into usable format*
```
/usr/local/packages/bbtools/reformat.sh in=dana.contigs.FREEZE.ltrdigest.mappedids.PAO.ltrs.fasta out1=dana.contigs.FREEZE.ltrdigest.mappedids.PAO.5ltr.fasta out2=dana.contigs.FREEZE.ltrdigest.mappedids.PAO.3ltr.fasta

rename sequentially 
awk '/^>/{print ">Gypsy_ltr_5_" ++i; next}{print}' < dana.FREEZE.final.Gypsy.5ltr.fasta > dana.FREEZE.final.Gypsy.5ltr.rn.fasta
awk '/^>/{print ">Gypsy_ltr_3_" ++i; next}{print}' < dana.FREEZE.final.Gypsy.3ltr.fasta > dana.FREEZE.final.Gypsy.3ltr.rn.fasta

split into individual files
seqkit split -i -O Gypsy.5.split dana.FREEZE.final.Gypsy.5ltr.rn.fasta
seqkit split -i -O Gypsy.3.split dana.FREEZE.final.Gypsy.3ltr.rn.fasta
```

*Pairwise alignment of LTRs using Needle*
```
for i in $(seq 1 1548); do needle -asequence Gypsy.5.split/*5_$i.fasta -bsequence Gypsy.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile Gypsy.needle+distmat/dana.Gypsy.$i.align.fasta; done

for i in $(seq 1 765); do needle -asequence PAO.5.split/*5_$i.fasta -bsequence PAO.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile PAO.needle+distmat/dana.PAO.$i.align.fasta; done
ls 
```

*Measure pairwise distance values using emboss dismat*
```
for i in $(seq 1 1548); do distmat -sequence Gypsy.needle+distmat/dana.Gypsy.$i.align.fasta -nucmethod 0 -outfile Gypsy.needle+distmat/dana.Gypsy.$i.uncorr.distmat.out; done

K2P
for i in $(seq 1 1548); do distmat -sequence Gypsy.needle+distmat/dana.Gypsy.$i.align.fasta -nucmethod 2 -outfile Gypsy.needle+distmat/dana.Gypsy.$i.k2p.distmat.out; done

cat PAO.needle+distmat/*uncorr.out > PAO.needle+distmat/PAO.distmat.uncorr.all.out

for f in PAO.needle+distmat/*uncorr.distmat.out; do awk '{print $2, $3}' $f | tail -n 2 | head -n 1 >> PAO.all.uncorr.distmat.out; done


sort -n -k2.12 PAO.all.uncorr.distmat.out > PAO.all.uncorr.distmat.sorted.out

sort -n -k2.14 Gypsy.all.uncorr.distmat.out > Gypsy.all.uncorr.distmat.sorted.out

paste <(grep '>' dana.FREEZE.final.PAO.5ltr.fasta ) <(grep '>' dana.FREEZE.final.PAO.3ltr.fasta ) <(grep '>' dana.FREEZE.final.PAO.5ltr.rn.fasta) <(awk '{print $1}' PAO.needle+distmat/PAO.all.uncorr.distmat.sorted.out) <(awk '{print $1}' PAO.needle+distmat/PAO.all.k2p.distmat.sorted.out) | column -t -o $'\t' | sed 's|>||g' | sed 's| ||g' > PAO.needle+distmat/PAO.all.final.distmat.out

paste <(grep '>' dana.FREEZE.final.Gypsy.5ltr.fasta ) <(grep '>' dana.FREEZE.final.Gypsy.3ltr.fasta ) <(grep '>' dana.FREEZE.final.Gypsy.5ltr.rn.fasta) <(awk '{print $1}' Gypsy.needle+distmat/Gypsy.all.uncorr.distmat.sorted.out) <(awk '{print $1}' Gypsy.needle+distmat/Gypsy.all.k2p.distmat.sorted.out) | column -t | sed 's/>//g' > Gypsy.needle+distmat/Gypsy.all.final.distmat.out

```


pipeline 2: LTRharvest + custom HMMer searches

echo "/local/projects-t3/LGT/Dananassae_2020/scripts/genometools-1.5.9/bin/gt ltrharvest -index /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.ltr.db -overlaps all -seed 30 -minlenltr 100 -maxlenltr 2000 -mindistltr 2000 -maxdistltr 15000 -v -gff3 dana.contigs.FREEZE.ltrharvest.allovl.gff -out dana.contigs.FREEZE.ltrharvest.allovl.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -cwd -N ltrharvest.alloverlaps


### Transcription of LGT regions <a name="dana.lgt.tx"></a>

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
```
grep "^>" ../dana.qm.merged.FREEZE.fasta | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat augustus.hints.codingseq ../dana.qm.merged.FREEZE.fasta > transcripts+genome.fasta 
/local/projects-t3/LGT/Dananassae_2020/scripts/salmon-latest_linux_x86_64/bin/salmon index -t transcripts+genome.fasta -d decoys.txt -i salmon_index -k 25 

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

**Search for numts by alignments between mito genome and nuclear genome
```
nucmer -l 50 --prefix dana.numt.firstpass dana.mito.complete.FREEZE.fasta /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta

tail -n +6 dana.numt.firstpass.r.coords | awk '{print $13}' | sort -n | uniq > numt.contigs.list

xargs samtools faidx /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta < numt.contigs.list > numt.contigs.fasta

nucmer -l 50 --prefix dana.numt.finalpass dana.mito.complete.FREEZE.fasta numt.contigs.fasta

```
