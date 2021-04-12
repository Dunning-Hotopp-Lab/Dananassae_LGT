# Extreme lateral gene transfer into a fly autosome

Eric S. Tvedte

2020-10-15

## Table of Contents
1. [Pre-processing sequencing data](#seq.prep)
2. [Genome assembly](#assemble)
3. [Post-assembly processing](#post)
4. [Anchoring assembly contigs](#anchor)  
5. [Genome annotation](#annotate)
6. [Summary statistics for D. ananassae genome assembly](#sumstats)
7. [Nuwt analysis](#nuwt)  
8. [Numt analysis](#numt)
9. [LTR retrotransposon analysis](#ltr)  
10. [Transcription of LGT regions](#lgt.tx)  


4. [Genome annotation](#annotate)

### 1. Pre-processing sequencing data <a name="seq.prep"></a>  
**ONT basecalling with guppy**
```
/usr/local/packages/guppy-3.1.5/bin/guppy_basecaller --input_path fast5_dir --save_path output_dir --config dna_r9.4.1_450bps_fast.cfg --fast5_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 --num_callers 8 --cpu_threads_per_caller 4  
```

### 2. Genome assembly <a name="assemble"></a>  
**Canu**  
```
canu -p output.prefix -d output.dir genomeSize=240m corOutCoverage=80 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq.gz -nanopore-raw minion.LIG.filter.lambda.fastq.gz
```

canu -p 'PB.CLR.het' 'genomeSize=60m' 'corMinCoverage=0' 'corOutCoverage=100' 'ovlMerSize=31' 'correctedErrorRate=0.035' 'utgOvlErrorRate=0.065' 'trimReadsCoverage=2' 'trimReadsOverlap=500' 'gridEngineThreadsOption=-pe thread THREADS' 'gridEngineMemoryOption=-l mem_free=MEMORY' 'gridOptions=-P jdhotopp-lab -q threaded.q' -pacbio '/autofs/projects-t3/LGT/Dananassae_2020/dana.nuwt/het.assembly/het.plus.unmapped.PB.CLR.fasta' canuIteration=0



**Flye** 
```
echo -e "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/packages/gcc/lib64\nexport PYTHONPATH=$PYTHONPATH:/usr/local/packages/flye-2.4.2/lib/python2.7/site-packages\n/usr/local/packages/flye-2.4.2/bin/flye -g 240m -t 24 -o /local/projects-t3/RDBKO/dana.flye/ --asm-coverage 60 --pacbio-raw /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.pbSequelII.raw.fastq.gz
```

### 3. Post-assembly processing <a name="post"></a>
**Polish Canu assembly with Arrow using Sequel II CLR data**
```
use python-3.5

echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/pbmm2 align dana.hybrid.contigs.fasta /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA/RANDD_20190301_S64018_PL100122512-1_C01.subreads.bam dana.hybrid.80X.contigs.mapped.pb.sqII_sorted.bam --sort -j 16 -J 8" | qsub -P jdhotopp-lab -l mem_free=50G -N pbmm2.align -q threaded.q -pe thread 16 -cwd -V

echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/arrow dana.pb.sqII.contigs.arrow1_sorted.bam -r dana.pb.sqII.contigs.fasta -o dana.pb.sqII.contigs.variants.gff -o dana.pb.sqII.contigs.polish.rd1.fasta -j 16 --noEvidenceConsensusCall reference" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 16 -N arrow -cwd -V

sed 's/|arrow//g' dana.pb.sqII.contigs.polish.rd1.fasta > dana.pb.sqII.contigs.polish.rd1.rn.fasta
```

**Merge Canu and Flye assemblies with Quickmerge** 
```
merge_wrapper.py canu_assembly.fasta flye_assembly.fasta -l 2000000 -ml 10000
```
**Polish merged genome**
```
pbmm2 align merged.genome.fasta pb.sequelII.subreads.bam merged.genome.mapped.pb.sequelII_sorted.bam --sort -j 16 -J 8

arrow merged.genome.mapped.pb.sequelII_sorted.bam -r merged.genome.fasta -o merged.genome.pb.sequelII.variants.gff -o merged.genome.polish.rd1.fasta -j 16 --noEvidenceConsensusCall reference  

minimap2 -ax map-pb -t 16 merged.genome.fasta pb.CCS.fasta.gz | samtools sort -o merged.genome.mapped.pb.CCS_sorted.bam  

pilon-1.22.jar --threads 16 --genome merged.genome.polish.rd1.fasta --output dana.genome.final.polish --unpaired merged.genome.mapped.pb.CCS_sorted.bam --changes --vcf --fix bases --mindepth 5 --K 85 --minmq 0 --minqual 35
```
**Purge haplotigs from assembly**
```
minimap2 -ax map-pb -t 16 dana.genome.final.polish.fasta pb.sequelII.fasta.gz | samtools sort -o merged.genome.mapped.pb.sequelII_sorted.bam  

purge_haplotigs hist -b merged.genome.mapped.pb.sequelII_sorted.bam -g dana.genome.final.polish.fasta -t 16  
purge_haplotigs cov -i merged.genome.mapped.pb.sequelII_sorted.bam.gencov -l 15 -m 150 -h 600 -o pb.sequelII.coverage.stats.csv  
purge_haplotigs purge -g dana.genome.final.polish.fasta -c pb.sequelII.coverage.stats.csv

blastn -query dana.genome.final.polish.fasta -db nt -outfmt ’6 qseqid staxids bitscore std’ -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 > merged.genome.polish.blast.out  
blobtools create -i dana.genome.final.polish.fasta -b merged.genome.mapped.pb.sequelII_sorted.bam -t merged.genome.polish.blast.out -o dana.sequelII.blobplot
blobtools view -i dana.sequelII.blobplot.blobDB.json -o dana.sequelII.blobplot.view
blobtools plot -i dana.sequelII.blobplot.blobDB.json -o dana.sequelII.blobplot.plot
```

### 4. Anchoring assembly contigs <a name="anchor"></a>
**Chromosome X, 2, 3**
*Initial BLAST search to determine contig orientation*
```
makeblastdb dana.genome.final.polish+purge.fasta -out dana.genome.final.polish+purge.fasta -dbtype nucl -parse_seqids
blastn -query Schaeffer2008.map.loci.fasta -db dana.genome.final.polish+purge.fasta -max_target_seqs 10 -max_hsps 10 -outfmt 6 > initial.blast.out
```
*Extract matching contigs, reverse complement if necessary*
```
samtools faidx dana.genome.final.polish+purge.fasta fwd.contig.name >> dana.chrX23.fasta
samtools faidx -i dana.genome.final.polish+purge.fasta rev.contig.name >> dana.chrX23.fasta
```
*Final BLAST search, using matched contigs as database and custom BLAST output*
```
makeblastdb dana.chrX23.fasta -out dana.chrX23.fasta -dbtype nucl -parse_seqids
blastn -query Schaeffer2008.map.loci.fasta -db dana.chrX23.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > dana.chrX23.Schaeffer2008.map.loci.blast.out
```
*Add remaining contigs to final assembly*
```
awk '{print $2}' dana.UMIGS.FREEZE.Schaeffer2008.map.loci.blast.out | sort -n | uniq > dana.UMIGS.FREEZE.chrX23.list
seqkit grep -v -f dana.UMIGS.FREEZE.chrX23.list dana.genome.final.polish+purge.fasta > dana.UMIGS.nonchrX23.fasta
cat dana.chrX23.fasta dana.UMIGS.nonchrX23.fasta > dana.UMIGS.FREEZE.fasta
```
**Chromosome Y**  
*Map male and female short reads*
```
seqkit subsample male 1
seqkit subsample male 2
seqkit subsample female 1
seqkit subsample female 2

bwa index dana.UMIGS.FREEZE.fasta
bwa mem -k 23 -t 8 dana.UMIGS.FREEZE.fasta female.R1.fq.gz female.R2.fq.gz
bwa mem -k 23 -t 8  dana.UMIGS.FREEZE.fasta male.R1.fq.gz male.R2.fq.gz

for f in /local/projects/JULIE/Dana*/ILLUMINA_DATA/*R1_trimmed.fastq.gz; do echo "bwa mem -k 23 -t 8 /local/projects-t3/RDBKO/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta $f ${f%_R*}_R2_trimmed.fastq.gz | samtools view -bho /local/projects-t3/RDBKO/dana.chrY/FREEZE/$(basename ${f%_R*})_output.bam" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 8 -N bwamem -cwd ; done
```
*Sort BAM and remove duplicates*
```
for f in *output.bam; do echo "java -Xmx2g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam I=$f O=${f%_o*}_sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/local/scratch/etvedte/" | qsub -P jdhotopp-lab -l mem_free=2G -N SortSam -cwd; done  
for f in *sorted.bam; do echo "java -Xmx10g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates I=$f O=${f%_s*}_dedup.bam M=${f%_s*}_dedup.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true" | qsub -P jdhotopp-lab -l mem_free=10G -N MarkDups -cwd; done
```
*Calculate sequencing depth, filter reads with MAPQ < 10*
```
samtools depth -Q 10 mel_f_new.bam mel_m_new.bam > mel_new.out
```
*Determine average and median female/male depth in 10kb windows (see Chang and Larracuente, 2018)* 
```perl
female
JULIE_20190729_K00134_IL100134730_MX83_L001
JULIE_20190729_K00134_IL100134731_MX84_L001
male
JULIE_20190702_K00134_IL100134728_MX81_L008
JULIE_20190702_K00134_IL100134729_MX82_L008

perl /local/projects-t3/RDBKO/scripts/Chang2019_frame_depth_new.pl mel_new.out
```

**Chromosome 4**
```
nucmer --maxmatch -l 1000 --prefix chr4.firstpass Leung2017.chr4.scaffolds.fasta dana.UMIGS.FREEZE.fasta  
show-coords -r chr4.firstpass.delta > chr4.firstpass.coords  
tail -n +6 chr4.firstpass.coords | awk '{print $13}' | sort -n | uniq > chr4.contigs.list  
seqkit grep -f chr4.contigs.list dana.UMIGS.FREEZE.fasta > dana.chr4.fasta

nucmer --maxmatch -l 1000 --prefix chr4.finalpass Leung2017.chr4.scaffolds.fasta dana.chr4.fasta
show-coords -r chr4.finalpass.delta > chr4.finalpass.coords  
mummerplot --color --postscript --prefix chr4.finalpass chr4.finalpass.delta  
```

### 5. Genome annotation <a name="annotate"></a>  
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
*parse data into usable format*
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

for i in $(seq 1 765); do needle -asequence PAO.5.split/*5_$i.fasta -bsequence PAO.3.split/*3_$i.fasta -gapopen 10 -gapextend 0.5 -aformat3 fasta -outfile PAO.needle+distmat/dana.PAO.$i.align.fasta; done
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

cat PAO.needle+distmat/*uncorr.out > PAO.needle+distmat/PAO.distmat.uncorr.all.out

for f in PAO.needle+distmat/*uncorr.distmat.out; do awk '{print $2, $3}' $f | tail -n 2 | head -n 1 >> PAO.all.uncorr.distmat.out; done


sort -n -k2.12 PAO.all.uncorr.distmat.out > PAO.all.uncorr.distmat.sorted.out

sort -n -k2.17 UMIGS.Ty3Gypsy.all.uncorr.distmat.out > UMIGS.Ty3Gypsy.all.uncorr.distmat.sorted.out

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
