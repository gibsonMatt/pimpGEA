#PBS -k oe #Keep the output and errors from the run
#PBS -l nodes=1:ppn=8,vmem=140gb,walltime=08:00:00
#PBS -M gibsomat@iu.edu
#PBS -m abe # Get notifications for the job
#PBS -N pimpGEA-stacks-map-bwa

cd /N/dc2/projects/gibsonTomato/pimpGEA/STACKS/demultiplexed_merged

module load bwa

for r1 in *.c.1.fq.gz
do
	echo "Mapping..."
	echo "$r1"
	toremove='.1.'
	newsubstr='.2.'
	r2="${r1/$toremove/$newsubstr}"
	
	echo "$r2"

	outname="/N/dc2/projects/gibsonTomato/pimpGEA/STACKS/map/mapped_reads/${r1::-8}.sam"
	echo "$outname"

	bwa mem -t 8 /N/dc2/projects/gibsonTomato/pimpGEA/data/genome/S_lycopersicum_chromosomes.3.00.fa $r1 $r2 > $outname

done
wait
