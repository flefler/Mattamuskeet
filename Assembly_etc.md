# Assemble, bin, and refine the cyanobacterial genomes

## Packages used 
[fastQC](https://github.com/s-andrews/FastQC)
[multiQC](https://multiqc.info/)
[singleM](https://github.com/wwood/singlem)
[megahit](https://github.com/voutcn/megahit)
[Quast](https://github.com/ablab/quast)
[BASALT](https://github.com/EMBL-PKU/BASALT)
[NanoPhase](https://github.com/Hydro3639/NanoPhase/blob/main/README.md)

# Short reads

## First lets run QC on the data
```
mamba activate Assemble_Bin
mkdir 01_QC
mkdir 99_logs
```

## fastQC and multiQC
```
mamba activate Assemble_Bin

SAMPLES=`cut -f 1 samples.txt | sed '1d'`

mamba activate Assemble_Bin
for SAMPLE in $SAMPLES; do
    fastqc 00_Reads/${SAMPLE}/*.fq.gz -o 01_QC
done

multiqc 01_QC -o 01_QC --interactive
```
## Now lets assemble with megahit
```
mkdir 03_ASSEMBLIES

SAMPLES=`cut -f 1 samples3.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
  megahit -1 00_Reads/${SAMPLE}/*_1.fq.gz \
          -2 00_Reads/${SAMPLE}/*_2.fq.gz \
          --out-dir 03_ASSEMBLIES/${SAMPLE} \
          --min-contig-len 1000 \
          -m 0.99
done
```
## Assess Assemblies with Quast
```
mkdir 04_ASSEMBLIES_QC

SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
  metaquast 03_ASSEMBLIES/${SAMPLE}/final.contigs.fa \
        --output-dir 04_ASSEMBLIES_QC/${SAMPLE} \
        --max-ref-number 0 \
        --threads 16
done
```
## BASALT
This bins genomes and then refines them
```
conda activate BASALT
mkdir 04_BASALT
```
You must move the files into a directory. The commands are NOT path friendly.
```
SAMPLES=`cut -f 1 samples_6.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
mkdir 04_BASALT/${SAMPLE}
cp 03_ASSEMBLIES/${SAMPLE}/final.contigs.fa 04_BASALT/${SAMPLE}/${SAMPLE}_final.contigs.fa
cp 00_Reads/${SAMPLE}/*_1.fq.gz 04_BASALT/${SAMPLE}/${SAMPLE}_1.fq.gz
cp 00_Reads/${SAMPLE}/*_2.fq.gz 04_BASALT/${SAMPLE}/${SAMPLE}_2.fq.gz
done

SAMPLES=`cut -f 1 /blue/hlaughinghouse/flefler/BLCC_Genomes/samples_2_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
	
	# Set a useful name for log files
	N="${SAMPLE}";

	# Set up command
	CMD="cd ${SAMPLE} && BASALT -a ${SAMPLE}_final.contigs.fa -s ${SAMPLE}_1.fq.gz,${SAMPLE}_2.fq.gz -t 28 -m 218 --min-cpn 80 --max-ctn 20 && cd .."
	
	# If you have SLURM engine, do something like this:
	sbatch -A hlaughinghouse -J ${N} -c 28 --mem=218G -o ${N}.o -e ${N}.e --export=ALL --mail-type=ALL --mail-user=flefler@ufl.edu -t 196:00:00  --wrap="${CMD}"
	
done
```
Yay were done! Onto the...

# Long reads

## Concat data
```
cd 4956-selected
cat 4956_1-2/fastq_pass/4956_2/*.fastq.gz 4956_1-2_Resequencing/fastq_pass/4956_2/*.fastq.gz  > finalreads/F210_nano.fastq.gz
```
## Nanophase 
### testing long reads only
This really wasnt that good :/
```
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    CMD="nanophase meta -t 28 -l finalreads/${SAMPLE}_nano.fastq.gz -o ${SAMPLE}_output"
    sbatch -A hlaughinghouse -J ${SAMPLE} -c 28 --mem=218G -o ${SAMPLE}_o -e ${SAMPLE}_e --export=ALL --mail-type=ALL --mail-user=flefler@ufl.edu -t 196:00:00  --wrap="${CMD}"
done
```
### testing hybrid
This was really good, one contig!!1!
```
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    R1="/blue/hlaughinghouse/flefler/BLCC_Genomes/00_Reads/${SAMPLE}/*_1.fq.gz"
    R2="/blue/hlaughinghouse/flefler/BLCC_Genomes/00_Reads/${SAMPLE}/*_2.fq.gz"
    CMD="nanophase meta -t 28 -l finalreads/${SAMPLE}_nano.fastq.gz --hybrid -1 ${R1} -2 ${R2} -o ${SAMPLE}_hybrid"
    sbatch -A hlaughinghouse -J ${SAMPLE} -c 28 --mem=218G -o ${SAMPLE}_h_o -e ${SAMPLE}_h_e --export=ALL --mail-type=ALL --mail-user=flefler@ufl.edu -t 196:00:00  --wrap="${CMD}"
done
```
## BASALT hybrid assembly
Not so good
```
cd ${SAMPLE} && BASALT -a ${SAMPLE}_contigs.fa -s ${SAMPLE}_1.fq,${SAMPLE}_2.fq -l ${SAMPLE}_lr.fq -t 8 -m 100 --min-cpn 80 --max-ctn 20 && cd ..
```

# Genome quality and other stastics
## Packages used 
[gtdbtk](https://github.com/Ecogenomics/GTDBTk)
[CheckM2](https://github.com/chklovski/CheckM2)
[CoverM](https://github.com/wwood/CoverM)
[SeqTK](https://github.com/lh3/seqtk)

## CheckM2
Used to assess completion and contamination
```
conda activate checkm2
mkdir 06_checkM

#Set job name
N="checkm2"
checkm2 predict --threads 28 --input 05_GENOMES --output-directory 06_checkM -x .fa.gz --quiet --remove_intermediates
```
## CoverM
Used to determine mean coverage
```
conda activate reassemble
mkdir 07_BAMFILES

for SAMPLE in $SAMPLES; do
    N=${SAMPLE}_coverm
    R1="/blue/hlaughinghouse/flefler/BLCC_Genomes/00_Reads/${SAMPLE}/*_1.fq.gz"
    R2="/blue/hlaughinghouse/flefler/BLCC_Genomes/00_Reads/${SAMPLE}/*_2.fq.gz"
    CMD="coverm genome -1 ${R1} -2 ${R2} --genome-fasta-files 05_GENOMES/${SAMPLE}_* --output-file 09_COVERM/${SAMPLE}_output.tsv \
    -x .gz --threads 8 --methods mean relative_abundance --bam-file-cache-directory process/bamcache"
    sbatch -A hlaughinghouse -J ${SAMPLE} -c 8 --mem=100G -o 99_logs/${SAMPLE}_h_o -e 99_logs/${SAMPLE}_h_e --export=ALL --mail-type=ALL --mail-user=flefler@ufl.edu -t 196:00:00  --wrap="${CMD}"
done
```

## GTDB
Used to assign taxonomy
```
conda activate gtdbtk-2.4.0
mkdir 08_gtdbtk
N="gtdb"
gtdbtk classify_wf --genome_dir 05_GENOMES --mash_db /blue/hlaughinghouse/flefler --out_dir 08_gtdbtk -x .fa.gz --cpus 28 --pplacer_cpus 28 --force
```

## SeqTK
Used to gather other stats, e.g., N50, length
```
seqkit stats -Ta *.fa.gz | csvtk tab2csv -o SeqTK_Output.csv
```
## Massage the outputs and make a single file
```
#CheckM file
csvtk tab2csv quality_report.tsv | csvtk rename -f 1 -n file -o quality_report.csv

#GTDB file
csvtk tab2csv gtdb/gtdbtk.bac120.summary.tsv | csvtk rename -f 1 -n file -o gtdb/gtdbtk.bac120.summary.csv

#coverM file
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    cat ${SAMPLE}_output.tsv | csvtk rename -t -f 1,2,3 -n file,meancoverage,relativeabundance -o ${SAMPLE}_output2.csv
done

csvtk concat -t *2.csv | csvtk tab2csv | csvtk filter2 -f '$file!="unmapped"' -o coverMoutput.csv
csvtk join -f 1 SeqTK_Output.csv coverm/coverMoutput.csv checkm/quality_report.csv gtdb/gtdbtk.bac120.summary.csv -o merged.csv
csvtk concat -k merged.csv process/selectedmags.csv -o genomeinfo_230524.csv
```
