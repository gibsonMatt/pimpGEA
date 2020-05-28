# Solanum pimpinellifolium Gene-Environment Association
Matt Gibson (gibsomat@indiana.edu)

Gibson, MJS & Moyle, LC. (2020). Regional differences in the abiotic environment contribute to genomic divergence within a wild tomato species. *Molecular Ecology*  https://doi.org/10.1111/mec.15477

Stacks pipeline. 140 *Solanum pimpinellifolium* individuals (assembled into two ddRAD libraries using PstI and EcoRI) were sequenced across two NextSeq flowcells (mid output; PE; 150bp).

## Pre-process data

### Filter adapters, fix bases in read overlap
* Seems easiest to filter adapters prior to running process radtags. So I'm doing that with `fastp` first, but not trimming for quality.


* See `./pimpGEA/scripts/adapter_trimming`

```
#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=50gb,walltime=20:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N trim_adapted_pimp_poolA

cd /N/dc2/projects/gibsonTomato/pimpGEA/data/rawdata/poolA


#ACAGTG
./N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001_noadapters.fastq -Q -L -w 6

#ATCACG
./N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-ATCACG_S1_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-ATCACG_S1_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ATCACG_S1_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ATCACG_S1_R2_001_noadapters.fastq -Q -L -w 6


#CGATGT
./N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-CGATGT_S2_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-CGATGT_S2_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-CGATGT_S2_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-CGATGT_S2_R2_001_noadapters.fastq -Q -L -w 6

#TGACCA
./N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-TGACCA_S4_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-TGACCA_S4_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TGACCA_S4_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TGACCA_S4_R2_001_noadapters.fastq -Q -L -w 6

#TTAGGC
./N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R2_001_noadapters.fastq -Q -L -w 6

#Same done for poolB
```

### Fix RE site with `recoverRE.py`

* See `./pimpGEA/scripts/recoverREsite`

```
#Read 1
python recoverRE.py --truesite TGCAG --startpos 6 --endpos 10 --out /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/pimpGEA/data/adapter_clean_rawdata/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001_noadapters.fastq

#Read 2
python recoverRE.py --truesite AATTC --startpos 1 --endpos 5 --out /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001_noadapters_refix.fastq /N/dc2/projects/gibsonTomato/pimpGEA/data/adapter_clean_rawdata/poolA/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001_noadapters.fastq

#...done for all index libraries
```


**At this point, I continue here with the stacks pipeline that I ran. See `README_pyrad_refix.md` for the new ipyrad pipeline that incorporates the RE-site recovered reads + more relaxed parameters.**

## Demultiplexing with `process_radtags`
* See `./pimpGEA/STACKS`

### Pool A
```
#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=2,vmem=20gb,walltime=20:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N process_radtags_poolA

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/

module load stacks

process_radtags -P -p /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolA/ -b /N/dc2/projects/gibsonTomato/pimpGEA/data/barcodes/poolA/poolA_stacks_barcodes.txt -o ./demultiplexed/poolA/ -r --inline_index --renz_1 pstI --renz_2 ecoRI -y gzfastq -s 5 --len_limit 30
```

### PoolB
```
#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=4,vmem=50gb,walltime=08:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N process_radtags_poolB

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/

module load stacks

process_radtags -P -p /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolB/ -b /N/dc2/projects/gibsonTomato/pimpGEA/data/barcodes/poolB/poolB_stacks_barcodes.txt -o ./demultiplexed/poolB/ -r --inline_index --renz_1 pstI --renz_2 ecoRI -y gzfastq -s 5 --len_limit 30
```

## Map with bwa

See `./pimpGEA/STACKS/map/` for scripts and mapped reads

```
#pimpGEA-stacks-map_no_rem.sh

#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=8,vmem=140gb,walltime=08:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N pimpGEA-stacks-map-bwa-no_rem

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/demultiplexed/demultiplexed_all/

module load bwa

for r1 in *.1.fq.gz
do
	echo "Mapping..."
	echo "$r1"
	toremove='.1.'
	newsubstr='.2.'
	r2="${r1/$toremove/$newsubstr}"
	
	echo "$r2"

	outname="/N/dc2/projects/gibsonTomato/pimpGEA/STACKS/map/mapped_reads_no_rem/${r1::-8}.sam"
	echo "$outname"

	bwa mem -t 8 /N/dc2/projects/gibsonTomato/pimpGEA/data/genome/S_lycopersicum_chromosomes.3.00.fa $r1 $r2 > $outname

done
wait
```

## Filter mapped reads

* See `./pimpGEA/STACKS/filter_mapped/` for scripts

```
#filter-bams-stacks-assembly_no_rem.sh

#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=30gb,walltime=05:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N filterBams-STACKS-assembly_no_rem

module load samtools

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/filter_mapped/

./filter_bams.sh /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/map/mapped_reads_no_rem/
```

```
#filter_bams.sh

#!/usr/bin/env bash
#Filters bams (sams in this case actually) based on specific criteria. In this case MQ > 3, no supplementary reads, and proper PE sequencing

cd $1

for f in *.sam
do
	echo "$f" 1>&2
	outname="/N/dc2/projects/gibsonTomato/pimpGEA/STACKS/filter_mapped/filtered_bams_no_rem/${f::-4}.filtered.bam"
	
	samtools view -f 2 -q 4 -h -b -o $outname $f &
done
wait

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/filter_mapped/filtered_bams_no_rem/

for f in *.filtered.bam
do

	outname="/N/dc2/projects/gibsonTomato/pimpGEA/STACKS/filter_mapped/filtered_bams_no_rem/${f::-4}.filtered.bam"
	samtools flagstat $f > "$outname.stats" &
done
wait
```

## Sort bams

* See `./pimpGEA/STACKS/sort_bams` for scripts

Code for sorting and indexing--I didnt index these bams, though. 

```
#sort_index_bams.sh

#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=1,vmem=60gb,walltime=08:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N sortNindexbams

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/filter_mapped/filtered_bams_no_rem

module load samtools

echo "Sorting..."
for f in *.bam
do
	echo "$f"
	outname="/N/dc2/projects/gibsonTomato/pimpGEA/STACKS/sorted_bams/${f::-4}.sorted.bam"
	
	samtools sort $f > $outname
done


#echo "Indexing..."

#cd /N/dc2/projects/gibsonTomato/pimpGEA/data/sorted_bams/

#for f in *.sorted.bam
#do
#	samtools index $f &
#done
#wait
```

## Stack Reference Mapping Pipeline

### `ref_map.pl`

* See `./pimpGEA/STACKS/gstacks/`

```
#pimpGEA-stacks-refmap.sh

#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=100gb,walltime=10:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N STACKS_refmap_no_rem

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks

#module load stacks

/N/dc2/projects/gibsonTomato/bin/ref_map.pl -T 6 --samples /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/sorted_bams/ --popmap /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/popmap.txt -o /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/outfiles/
```


### `populations`

* See `./pimpGEA/STACKS/populations/`

I pretty much requested all output files with minimal filters applied:

1. Site needs to be in at least 1 population
2. Within a population, the site much be typed in at least 20% of individuals. 

```
#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=60gb,walltime=2:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N STACKS_populations

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/populations

#module load stacks

/N/dc2/projects/gibsonTomato/bin/populations -P /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/outfiles/ -M /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/popmap.txt -t 6 -O /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/populations/outfiles1/ -p 1 -r 0.20 --hwe --fstats --fasta_loci --fasta_samples --vcf --vcf_haplotypes --structure --fasta_samples_raw
```


## Vcf filtering

**First, I changed 0772-101 to 0722-101, since it was a mislabeled plant from the get-go**

### Remove .filtered.sorted from indv names
```
sed -i '' 's/.filtered.sorted//g' populations.snps.vcf
```


### Initial stats:

Total # snps:  366,466
All biallelic


### Filter 1: Remove chromosome 0

Total # snps: 359,796

```
grep -Ev 'SL3.0ch00' populations.snps.vcf > populations.snps.filter1.vcf
```

### Filter 2: For genotype allele depth greater than 4

Total # of snps: 359,796

Doesn't remove sites, but changes genotypes to missing if DP < 4

```
vcftools --vcf populations.snps.filter1.vcf --recode --recode-INFO-all --out populations.snps.filter2 --minDP 4
```

### Split up by species in TASSEL GUI:

1. With all species
2. With just arc
3. With just chl
4. With just chm
5. With just hab
6. With just pimp

### Filter 3: For missing data

Remove sites when fewer than 20% of individuals are typed. 

Total # of snps full dataset: 109,184
Total # of snps just pimp: 112,252
```
#For all
vcftools --vcf populations.snps.filter2.alltaxa.vcf --recode --recode-INFO-all --out populations.snps.filter3.alltaxa --max-missing-count 128

#For pimp
vcftools --vcf populations.snps.filter2.pimp.vcf --recode --recode-INFO-all --out populations.snps.filter3.pimp --max-missing-count 120
```

### Filter 4: Remove invariant sites:

Total # of snps full dataset: 106,440
Total # of snps just pimp: 44,533

```
vcftools --vcf populations.snps.filter3.alltaxa.recode.vcf  --non-ref-ac-any 1 --recode --recode-INFO-all --out populations.snps.filter4.alltaxa.recode

vcftools --vcf populations.snps.filter3.pimp.recode.vcf  --non-ref-ac-any 1 --recode --recode-INFO-all --out populations.snps.filter4.pimp
```



### Filter 5: Remove highly heterozygous snps

Total # of snps: 44,243

**Done in TASSEL** Max heterozygous = 0.6



### Filter 6: Remove individuals with >60% missing data + remove invariant sites


142 individuals remaining. 
Total # of SNPs: 44,064

```
vcftools --vcf populations.snps.filter5.pimp.vcf --remove-indv 365911-101 --remove-indv 1601-101 --remove-indv 1416-101 --remove-indv 1987-101 --remove-indv 2649-101 --remove-indv 2933-101  --remove-indv 390695-101 --remove-indv 503517-101 --recode --recode-INFO-all --out populations.snps.filter6.pimp


#Remove invariant sites post indv filter
vcftools --vcf populations.snps.filter6.pimp.recode.vcf  --non-ref-ac-any 1 --recode --recode-INFO-all --out populations.snps.filter7.pimp
```

### Filter 8: Remove SNPs in full LD

Total # of SNPs: 17,358

```
bcftools +prune -l 0.9 -w 1000 populations.snps.filter7.pimp.recode.vcf -Ov -o populations.snps.filter8.pimp.vcf
```
