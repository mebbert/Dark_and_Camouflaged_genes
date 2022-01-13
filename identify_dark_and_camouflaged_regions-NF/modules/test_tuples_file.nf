process test_tuples {
	
	label 'local'

	input:
		val(gvcfs_result_dir)
		path(gatk_bed)
		val(result_dir_main)
		val(mask_ref_prefix)

	output:
		tuple path(gvcf_list), val(region), val(result_dir), path(ref) , emit: test_tuples_out

	script:
	"""
	for dir in $gvcfs_result_dir/*
	do
		ploidy=\$(basename \$dir)
		repeat=\$((\${ploidy##*_} / 2))
		result_dir="\${result_dir_main}/genotyped_by_region/\${ploidy}"
		gvcf_list="\${ploidy}.gvcfs.list"
		find \$dir -name "*.gvcf" > \$gvcf_list

		ref="\${camo_mask_ref_prefix}.\${ploidy}.fa"

		awk "\$5 == \$repeat {print \$1\":\"\$2\"-\"\$3}" $gatk_bed | \
			while read region
			do
				echo \$region
			done

	done
	"""
}
