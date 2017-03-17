#!/bin/bash

runs=$(ls | grep -v .fasta)

for run in $runs; do
	#echo $run
	cd $run
	samples=$(ls | grep _NS5A.fasta)
	for sample in $samples; do
		#echo $sample
		
		#tr '^M' '' < $sample > ../$sample
		sed 's/\r/\n/g' < $sample > ../$sample
		
		sample_short=$(echo $sample | sed 's/_NS5A.fasta//')
		#echo $sample_short
		ln -sv /data/MiSeq/MiSeqOutput/${run}_*/Data/Intensities/BaseCalls/${sample_short}*.gz ../
	done
	cd ..
done
exit