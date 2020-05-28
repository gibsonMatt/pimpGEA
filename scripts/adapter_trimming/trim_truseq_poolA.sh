#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=50gb,walltime=20:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N trim_adapted_pimp_poolA

#Trim adapters and correct overlap seqs

cd /N/dc2/projects/gibsonTomato/pimpGEA/data/rawdata/poolA


#ACAGTG
/N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ACAGTG_S5_R2_001_noadapters.fastq -Q -L -w 6 -c

#ATCACG
/N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-ATCACG_S1_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-ATCACG_S1_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ATCACG_S1_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-ATCACG_S1_R2_001_noadapters.fastq -Q -L -w 6 -c


#CGATGT
/N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-CGATGT_S2_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-CGATGT_S2_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-CGATGT_S2_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-CGATGT_S2_R2_001_noadapters.fastq -Q -L -w 6 -c 

#TGACCA
/N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-TGACCA_S4_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-TGACCA_S4_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TGACCA_S4_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TGACCA_S4_R2_001_noadapters.fastq -Q -L -w 6 -c

#TTAGGC
/N/dc2/projects/gibsonTomato/fastp -i GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R1_001.fastq -I GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R1_001_noadapters.fastq -O ../../adapter_clean_rawdata/GSF1870-Moyle-PoolA-Index-TTAGGC_S3_R2_001_noadapters.fastq -Q -L -w 6 -c
