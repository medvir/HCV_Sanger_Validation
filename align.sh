#!/bin/bash

n_reads=200000

### make list of files to analyse
list=$(ls ./*_NS5A_disamb.fasta)

#list=("1000282841_NS5A_disamb.fasta 1000335797_NS5A_disamb.fasta 1000336164_NS5A_disamb.fasta")


### loop over list	
for i in $list; do
	ref=$i
	name=$(basename $i | sed 's/_NS5A_disamb.fasta//')
	fastq_file=$(ls ./*gz | grep $name)

	seqtk sample $fastq_file $n_reads > ${name}_reads.fastq
	n_sample=$(wc -l ${name}_reads.fastq | cut -f 1 -d " ")
	n_sample=$(($n_sample / 4))
	
	echo
	echo sample $name
	if [ -s $ref ]; then
		echo reference $ref
		echo fastq file $fastq_file, $n_sample reads
		echo "*******************************************************************************"
	else
		echo ERROR: No reference
		exit
	fi
		
	### align with smalt to reference ${name}_reads.fastq
	bwa index $ref
	samtools faidx $ref
	bwa mem -t 12 $ref ${name}_reads.fastq | samtools view -u -@ 4 - | samtools sort -@ 4 -O bam -T tmp -o ${name}.bam -
	samtools index ${name}.bam

	### create vcf with lofreq
	rm -f ${name}_lofreq.vcf
	lofreq call-parallel --pp-threads 10 -f $ref -o ${name}_lofreq.vcf ${name}.bam
	
	### calculate depth
	samtools depth ${name}.bam > ${name}.depth
	
	### remove temporary files of this iteration
	rm -f ${name}.sam
	rm -f ${name}.vcf
	rm -f ${name}_smalt_index.*
	rm -f ${name}_*_cons.fasta.fai
	rm -f ${name}_reads.fastq
	rm -f ${name}_reads_contigs.fasta
	rm -rf ${name}
done