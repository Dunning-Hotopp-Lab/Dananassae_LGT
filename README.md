# INSERT TITLE HERE

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
canu -p output.prefix -d output.dir genomeSize=240m corOutCoverage=60 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq.gz -nanopore-raw minion.LIG.filter.lambda.fastq.gz
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
echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/arrow -j 16 dana.hybrid.80X.contigs.mapped.pb.sqII_sorted.bam -r dana.hybrid.80X.contigs.fasta -o dana.hybrid.80X.contigs.mapped.pb.sqII_variants.gff -o dana.hybrid.80X.contigs.arrow.polished.fasta" | qsub -P jdhotopp-lab -l mem_free=50G -N smrttools.arrow -q threaded.q -pe thread 16 -cwd -V
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

**NUCMER**
nucmer --maxmatch --prefix FREEZE.chrs -l 100 dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta


**QUAST**
```
use python-3.5
echo "/home/etvedte/scripts/quast-5.0.2/quast.py --large --k-mer-stats --fragmented --threads 12 /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd1/dana.hybrid.80X.contigs.arrow.polished.rn.fasta /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa -r /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic_scaffolds.fna -g /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.gff -o /local/projects-t3/RDBKO/dana.correctness/quast" | qsub -P jdhotopp-lab -l mem_free=10G -N quast-LG -q threaded.q -pe thread 12 -cwd -V
```

**KAT**
```
use kat-2.4.0
echo "kat comp -t 16 -o FREEZE_vs_illuminaR1 /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R1.fastq /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -N kat.comp.illumina -cwd -V

```

**Rename FASTA**  
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done



### Genome annotation <a name="annotate"></a>  
**Map short RNA reads** 
```
hisat2-build polished.contigs.fasta polished.contigs.hisat2  
hisat2 -p 8 --max-intronlen 300000 -x polished.contigs.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  

for f in /local/projects-t3/RDBKO/sequencing/Dana_illumina_RNA_SRA/*.fastq; do echo "hisat2 -p 8 --max-intronlen 300000 -x dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2 -U $f | samtools view -bho ${f%_1*}_output.bam -" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N hisat2 -cwd; done
```

**Map long RNA reads**
```
echo -e "/usr/local/packages/minimap2-2.10/bin/minimap2 -ax splice -uf -k14 -t 8 -G 300000 /local/projects-t3/RDBKO/dana.postassembly/arrow/dana.hybrid.80X.contigs.arrow.polished.fasta /local/projects-t3/RDBKO/sequencing/Dana_directRNA.fastq | samtools view -bho output_bam -" | qsub -q threaded.q -pe thread 8 -P jdhotopp-lab -l mem_free=10G -N minimap2 -cwd 
```

**Sort BAM**
```
java -Xmx2g -jar picard.jar SortSam I=mapped.bam O=sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
```

**Find and Mask repeats**
Build repeat database
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
```

**Soft-masking genome using repeat families
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
cat chr4.coords | tail -n +6 | awk '{print $13}' | sort -n | uniq > chr4.contigs.list
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

**Mugsy**
```
source /local/projects/angiuoli/mugsy/mugsyenv.sh
/local/projects/angiuoli/mugsy/mugsy --directory /local/projects-t3/RDBKO/dana.LGT/80X.polished.rd1/mugsy --prefix ecoli 1.ecoli.fasta 2.ecoli.fasta 3.ecoli.fasta
source /home/jdhotopp/bin/jsahl_mugsy_to_tree_dir/pythonenv.sh
/home/jdhotopp/bin/jsahl_mugsy_to_tree_dir/process_maf.sh 432481_433504_CDS.maf
```
### Transcription of LGT regions <a name="dana.lgt.tx"></a>

**Map short RNA reads** 
```
hisat2-build polished.contigs.fasta polished.contigs.hisat2  
hisat2 -p 8 --max-intronlen 300000 -x polished.contigs.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  

hisat2-build /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2

for f in /local/projects-t3/RDBKO/sequencing/Dana_illumina_RNA_SRA/*.fastq; do echo "hisat2 -p 8 --max-intronlen 300000 -x /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/braker/FREEZE/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.hisat2 -U $f | samtools view -bho ${f%_1*}_output.bam -" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N hisat2 -cwd; done
```

**Sort BAM**
```
java -jar picard.jar SortSam I=output.bam O=sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

for f in *output.bam; do echo "java -Xmx2g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$f O=${f%_o*}_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/etvedte/" | qsub -P jdhotopp-lab -l mem_free=2G -N SortSam -cwd; done  

```
**Count reads in each LGT region**
#count how many reads map to LGT region, include conditional statement to ensure coordinates are in correct order

BAM=/local/aberdeen2rw/julie/ben/Dana_transcriptome/hisat2/SRR921454_sorted.bam
LIST=/local/aberdeen2rw/julie/ben/Dana_transcriptome/hisat2/SRR921454.list
OLD_CSV=/local/aberdeen2rw/julie/ben/Dana_transcriptome/hisat2/Dana_LGT_transcription_6.csv
NEW_CSV=/local/aberdeen2rw/julie/ben/Dana_transcriptome/hisat2/Dana_LGT_transcription_7.csv

while read Line
do
contig=$(echo $Line | awk '{print $2}')

start=$(echo $Line | awk '{print $3}')

stop=$(echo $Line | awk '{print $4}')

if [ "$stop" -gt "$start" ]
then
    samtools view -c "$BAM" $contig:$start-$stop >> "$LIST"
else
    samtools view -c "$BAM" $contig:$stop-$start >> "$LIST"
fi

pr -mts "$OLD_CSV" "$LIST" > "$NEW_CSV"

done < /local/aberdeen2rw/julie/ben/Dana_transcriptome/hisat2/wAna.LGT.only.filter.samtools.coords

```
