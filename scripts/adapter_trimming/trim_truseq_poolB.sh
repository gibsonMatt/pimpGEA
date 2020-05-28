#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=50gb,walltime=20:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N trim_adapted_pimp_poolB

#Trim adapters and correct overlap seqs

cd /N/dc2/projects/gibsonTomato/pimpGEA/data/rawdata/poolB


#365960
/N/dc2/projects/gibsonTomato/fastp -i GSF1486-1780-Moyle-PoolB-365960-101_S5_R1_001.fastq -I GSF1486-1780-Moyle-PoolB-365960-101_S5_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1486-1780-Moyle-PoolB-365960-101_S5_R1_001_noadapters.fastq -O ../../GSF1486-1780-Moyle-PoolB-365960-101_S5_R2_001_noadapters.fastq -Q -L -w 6 -c

#1924
/N/dc2/projects/gibsonTomato/fastp -i GSF1486-1780-Moyle-PoolB-1924-101_S4_R1_001.fastq -I GSF1486-1780-Moyle-PoolB-1924-101_S4_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1486-1780-Moyle-PoolB-1924-101_S4_R1_001_noadapters.fastq -O ../../GSF1486-1780-Moyle-PoolB-1924-101_S4_R2_001_noadapters.fastq -Q -L -w 6 -c


#1576
/N/dc2/projects/gibsonTomato/fastp -i GSF1486-1780-Moyle-PoolB-1576-101_S3_R1_001.fastq -I GSF1486-1780-Moyle-PoolB-1576-101_S3_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1486-1780-Moyle-PoolB-1576-101_S3_R1_001_noadapters.fastq -O ../../GSF1486-1780-Moyle-PoolB-1576-101_S3_R2_001_noadapters.fastq -Q -L -w 6 -c

#1381
/N/dc2/projects/gibsonTomato/fastp -i GSF1486-1780-Moyle-PoolB-1381-101_S1_R1_001.fastq -I GSF1486-1780-Moyle-PoolB-1381-101_S1_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1486-1780-Moyle-PoolB-1381-101_S1_R1_001_noadapters.fastq -O ../../GSF1486-1780-Moyle-PoolB-1381-101_S1_R2_001_noadapters.fastq -Q -L -w 6 -c

#0122
/N/dc2/projects/gibsonTomato/fastp -i GSF1486-1780-Moyle-PoolB-0122-101_S2_R1_001.fastq -I GSF1486-1780-Moyle-PoolB-0122-101_S2_R2_001.fastq -o ../../adapter_clean_rawdata/GSF1486-1780-Moyle-PoolB-0122-101_S2_R1_001_noadapters.fastq -O ../../GSF1486-1780-Moyle-PoolB-0122-101_S2_R2_001_noadapters.fastq -Q -L -w 6 -c
