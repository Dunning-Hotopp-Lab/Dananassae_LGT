# INSERT TITLE HERE

Eric S. Tvedte

2020-1-9

The repository contains Supplementary Data for the manuscript, including Tables, Figures, and Files.

## Table of Contents
1. [Pre-processing MinION sequencing data](#seq.prep)
2. [Genome assembly](#assemble)
3. [Post-assembly processing](#ecoli.uni)
4. [Genome annotation](#annotate)
5. [BUSCO analysis](#busco)
6. [Characterization of euchromatic regions of D. ananassae](#dana.chrom.map)  
7. [Characterization of Y contigs in D. ananassae](#dana.y)  
8. [Characterization of LGT contigs in D. ananassae](#dana.lgt)  

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

### Genome polishing  
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
pilon_iter.sh canu/assembly.contigs.fasta illumina.reads.R1.fastq illumina.reads.R2.fastq canu/assembly.trimmedReads.fasta  
circlator minimus2 pilon5.fasta circularise.fasta  
circlator fixstart --genes_fa ecoli.dnaA.DNA.fasta circularise.fasta rotated.fasta  

**Rename FASTA**  
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done



**Map short read RNA**  
hisat2-build polished.contigs.fasta polished.contigs.hisat2  
hisat2 -p 8 --max-intronlen 300000 -x polished.contigs.hisat2 -U reads.fastq.gz | samtools view -bho output.bam -  

**Map direct RNA sequencing**
{bash, eval = F}
echo -e "/usr/local/packages/minimap2-2.10/bin/minimap2 -ax splice -uf -k14 -t 8 -G 300000 /local/projects-t3/RDBKO/dana.postassembly/arrow/dana.hybrid.80X.contigs.arrow.polished.fasta /local/projects-t3/RDBKO/sequencing/Dana_directRNA.fastq | samtools view -bho output_bam -" | qsub -q threaded.q -pe thread 8 -P jdhotopp-lab -l mem_free=10G -N minimap2 -cwd 

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
/usr/local/packages/repeatmodeler-1.0.11/RepeatModeler -database BrugiaPahangiNuclearGenome.fa
    /usr/local/packages/repeatmasker-4.0.7/RepeatMasker -lib rugiaPahangiNuclearGenome.fa-families.fa /local/aberdeen2rw/julie/JM_dir/PahangiPilonFASTA/Repeats/BrugiaPahangiNuclearGenome.fa
```

