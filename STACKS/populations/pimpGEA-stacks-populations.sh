#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=60gb,walltime=2:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N STACKS_populations

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/populations

#module load stacks

/N/dc2/projects/gibsonTomato/bin/populations -P /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/outfiles/ -M /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/popmap.txt -t 6 -O /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/populations/outfiles1/ -p 1 -r 0.20 --hwe --fstats --fasta_loci --fasta_samples --vcf --vcf_haplotypes --structure --fasta_samples_raw
