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
