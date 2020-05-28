#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=6,vmem=100gb,walltime=10:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N STACKS_refmap_no_rem

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks

#module load stacks

/N/dc2/projects/gibsonTomato/bin/ref_map.pl -T 6 --samples /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/sorted_bams/ --popmap /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/popmap.txt -o /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/gstacks/outfiles/ --unpaired
