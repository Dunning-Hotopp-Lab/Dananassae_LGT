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
echo "minimap2 -xmap-pb purged.fa /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.pbSequelII.raw.fastq.gz | gzip -c - > dana.hybrid.80X.purged.mappedsqII.paf.gz" | qsub -P jdhotopp-lab -l mem_free=10G -N minimap2 -cwd
/home/etvedte/scripts/purge_dups/bin/pbcstat dana.hybrid.80X.contigs.arrow.polished.mappedhifi.paf.gz
/home/etvedte/scripts/purge_dups/bin/calcuts PB.stat > cutoffs 2> calcuts.log
/home/etvedte/scripts/purge_dups/scripts/hist_plot.py PB.stat hist.out.pdf
```

**QUAST**
```
use python-3.5
echo "/home/etvedte/scripts/quast-5.0.2/quast.py --large --k-mer-stats --fragmented --threads 12 /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd1/dana.hybrid.80X.contigs.arrow.polished.rn.fasta /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa -r /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic_scaffolds.fna -g /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.gff -o /local/projects-t3/RDBKO/dana.correctness/quast" | qsub -P jdhotopp-lab -l mem_free=10G -N quast-LG -q threaded.q -pe thread 12 -cwd -V
```

**KAT**
```
use kat-2.4.0
echo "kat comp -t 16 -o purged_vs_illuminaR1 /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.illuminaPE_R1.fastq.gz /local/projects-t3/RDBKO/dana.postassembly/purge_dups/purged.fa" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -N kat.comp.illumina -cwd -V
```

**Rename FASTA**  
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done


### Genome annotation <a name="annotate"></a>  
**Map short RNA reads** 
```
hisat2-build polished.contigs.fasta polished.contigs.hisat2  
hisat2 -p 8 --max-intronlen 300000 -x polished.contigs.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  
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

**Soft-masking genome using repeat families
```
/usr/local/packages/repeatmasker-4.0.7/RepeatMasker -lib rugiaPahangiNuclearGenome.fa-families.fa /local/aberdeen2rw/julie/JM_dir/PahangiPilonFASTA/Repeats/BrugiaPahangiNuclearGenome.fa
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

**Sort SAM and remove duplicates**
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
