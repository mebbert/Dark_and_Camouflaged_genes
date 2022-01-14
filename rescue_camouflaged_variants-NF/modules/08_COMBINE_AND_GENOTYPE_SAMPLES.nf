process step_08_COMBINE_AND_GENOTYPE_SAMPLES {
	
	label 'step_08'

	input:
		path(realign_bed)
		path(gatk_bed)
		path(gatk_jar)
		val(ref_prefix)
		path(filtered_list)

	output:
		path '*.dark.low_depth.bed', emit: low_depth_out
		path '*.dark.low_mapq.bed', emit: low_mapq_out

	script:
	"""
		bash combine_and_genotype.sh $ $ $ $ $ $
	"""
}

process step_08_PREP {

	label 'local'

	input:
		path(filtered_bam_list)

	output:
		path 'split_list*', emit: filtered_lists
		tuple file(), val(), file(), file(), emit: set_of_values

	script:
	"""
		for dir in \${}/* 
		do
			ploidy=\$(basename \$dir)
			repeat=\$((\${ploidy##*_} / 2))
			
			gvcf_list="\${ploidy}.gvcfs.list"
			find \$dir -name "*.g.vcf" > \$gvcf_list

			REF="\${camo_mask_ref_prefix}.\${ploidy}.fa"
			awk "\$5 == $repeat {print \$1\":\"\$2\"-\"\$3}" $GATK_bed | \
				while read region
				do
					echo \$region
					echo \$REF
					echo \$gvcf_list
					echo \$result_dir
				done

		done
	"""
}
