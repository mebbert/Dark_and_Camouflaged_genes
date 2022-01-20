#!/bin/bash


crams=$1
ref_tag=$2

#mkdir -p $out

echo ${crams}
echo ${crams[@]}
echo 'this is a test'

for cram in ${crams[@]}
do
	echo ${cram}
	samtools view ${cram} >> 'testing_out.cram'
done	

