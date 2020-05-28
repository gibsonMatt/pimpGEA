#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=30gb,walltime=05:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N filterBams-STACKS-assembly_no_rem

module load samtools

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/filter_mapped/

./filter_bams.sh /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/map/mapped_reads_no_rem/


