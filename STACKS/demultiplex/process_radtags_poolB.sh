#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=1,vmem=50gb,walltime=08:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N process_radtags_poolB

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/

module load stacks

process_radtags -P -p /N/dc2/projects/gibsonTomato/pimpGEA/data/recoverREsite/poolB/ -b /N/dc2/projects/gibsonTomato/pimpGEA/data/barcodes/poolB/poolB_stacks_barcodes.txt -o ./demultiplexed/poolB/ -r --inline_index --renz_1 pstI --renz_2 ecoRI -y gzfastq -s 5 --len_limit 30
