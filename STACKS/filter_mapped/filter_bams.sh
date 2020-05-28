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
