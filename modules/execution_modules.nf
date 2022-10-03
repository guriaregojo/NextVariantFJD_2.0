




process BS_CHECK {
	input:
		val project
		val baseuser
		path samples
		val analysis

	output:
		path "projects.txt", emit: bsproyects
		path "controlsamples.txt", emit: controlsamples
		path "samples2analyce.txt", emit: samples2analyce

	script:
		def baseuser_field = baseuser ? "--config ${baseuser} " : ''
		
		if(samples && analysis.contains("C"))
			"""
			bs list projects -f csv -F Name ${baseuser_field} > projects.txt


			if grep -Fq ${project} projects.txt; then
				bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName \
				${baseuser_field} --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi


			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
					exit 1
				fi
			done
			"""


		else if(samples) 
			"""
			bs list projects -f csv -F Name ${baseuser_field} > projects.txt


			if grep -Fq ${project} projects.txt; then
				bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName \
				${baseuser_field} --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi


			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
					exit 1
				fi   
			done


			cat samples2analyce.txt > controlsamples.txt
			"""
		

		else 
			"""
			bs list projects -f csv -F Name ${baseuser_field} > projects.txt


			if grep -Fq ${project} projects.txt; then
				bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName \
				${baseuser_field} --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			cat controlsamples.txt > samples2analyce.txt
			"""
		
}












process LOCAL_CHECK {
	input:
		path input
		path samples
		val analysis

	output:
		path "controlsamples.txt", emit: controlsamples
		path "samples2analyce.txt", emit: samples2analyce

	script:
		if ( analysis.contains("M") )      { extension_local_check = "(.fq|.fq.gz|.fastq|.fastq.gz)" }
		else if ( analysis.contains("Q") ) { extension_local_check = "(.bam)" }
		else if ( analysis.contains("S") ) { extension_local_check = "(.bam)" }
		else if ( analysis.contains("G") ) { extension_local_check = "(.bam)" }
		else if ( analysis.contains("C") ) { extension_local_check = "(.bam)" }
		else if ( analysis.contains("X") ) { extension_local_check = "(.bam)" }
		else if ( analysis.contains("A") ) { extension_local_check = "(.vcf|.vcf.gz)" }
		else if ( analysis.contains("N") ) { extension_local_check = "(.tsv|.bed)" }

		if(samples && analysis.contains("C"))
			"""
			ls ${input} | grep "[${extension_local_check}]\$" -P | sed 's/_.*//' | sort | uniq > controlsamples.txt
			

			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep -o \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the input folder '${input}'\n"
					exit 1
				fi
			done
			"""


		else if(samples) 
			"""
			ls ${input} | grep "[${extension_local_check}]\$" -P | sed 's/_.*//' | sort | uniq > controlsamples.txt
			

			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep -o \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the input folder '${input}'\n"
					exit 1
				fi
			done


			cat samples2analyce.txt > controlsamples.txt
			"""
		

		else 
			"""
			ls ${input} | grep "${extension_local_check}\$" -P | sed 's/_.*//' | sed 's/\\..*//' | sort | uniq > controlsamples.txt


			cat controlsamples.txt > samples2analyce.txt
			"""
		
}




process BS_COPY {

	errorStrategy 'retry'
	maxRetries 3

	input:
		val project
		val sample2download
		val baseuser

	output:
		tuple \
			val(sample2download_config), \
			path("${sample2download_config}_R1.fastq.gz"), \
			path("${sample2download_config}_R2.fastq.gz"), emit: fastq

	script:
		sample2download_config = sample2download[0]
		def baseuser_config = baseuser ? "--config ${baseuser} " : ''
		"""
		echo ${sample2download_config} > ${sample2download_config}.sample.txt
		bs download biosample -q -n "${sample2download_config}" --exclude '*' --include '*.fastq.gz'


		cat ${sample2download_config}*/*_R1*fastq.gz > ${sample2download_config}_R1.fastq.gz
		cat ${sample2download_config}*/*_R2*fastq.gz > ${sample2download_config}_R2.fastq.gz
		"""
}



process FASTQ_CONCATENATION {
	input:
		path inputdir
		val sample2concat

	output:
		tuple \
			val(sample2concat_config), \
			path("${sample2concat_config}_R1.fastq.gz"), \
			path("${sample2concat_config}_R2.fastq.gz"), emit: fastq

	script:
		sample2concat_config = sample2concat[0]
		"""
		echo ${sample2concat_config}
		nR1="\$(ls ${inputdir}/${sample2concat_config}*R1*.f*q* | wc -l)"
		nR2="\$(ls ${inputdir}/${sample2concat_config}*R2*.f*q* | wc -l)"

		if [[ \$nR1 != \$nR2 ]]; then
			>&2 echo "Error: Number of R1 (\${nR1}) files do not match with the number of R2 (\${nR2}) files for the sample ${sample2concat_config}.\n"
			exit 1
		fi

		if [[ \$nR1 == 1 ]]; then
			ln -s ${inputdir}/${sample2concat_config}*R1*.f*q* ${sample2concat_config}_R1.fastq.gz
			ln -s ${inputdir}/${sample2concat_config}*R2*.f*q* ${sample2concat_config}_R2.fastq.gz
		elif [[ \$nR1 -gt 1 ]]; then
			cat ${inputdir}/${sample2concat_config}*R1*.f*q* > ${sample2concat_config}_R1.fastq.gz
			cat ${inputdir}/${sample2concat_config}*R2*.f*q* > ${sample2concat_config}_R2.fastq.gz 
		else
			>&2 echo "Error: wrong number of fastq. fastq files must match the pattern '*R[12]*.f*q*'\n"
			exit 1
		fi
		"""
}

		











process BWA {
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(forward), path(reverse)
		path ref
		path bwa_amb
		path bwa_ann
		path bwa_pac
		path bwa_bwt
		path bwa_sa
		val bwa_threads

	output:
		tuple \
			val(sample), \
			path("${sample}.mapped.bam"), emit: mapped_bam

	script:
		// sample  = fastq[0]
		// forward = fastq[1]
		// reverse = fastq[2]
		"""
		bwa mem -v 3 -t ${bwa_threads} -Y \\
		${ref} \\
		${forward} \\
		${reverse} |  samtools view -1 > ${sample}.mapped.bam

		"""
}





process FASTQTOSAM {
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(forward), path(reverse)
		path scratch

	output:
		tuple \
			val(sample), \
			path("${sample}.unmapped.bam"), emit: unmapped_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_FastqToSam" : ""
		"""

		gatk FastqToSam ${scratch_field} \
		--FASTQ ${forward} \
		--FASTQ2 ${reverse} \
		--OUTPUT ${sample}.unmapped.bam \
		--READ_GROUP_NAME ${sample} \
		--SAMPLE_NAME ${sample} \
		--PLATFORM illumina 
		"""
}



process MERGEBAMALIGNMENT {
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(mapped_bam), path(unmapped_bam)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch

	output:
		tuple \
			val(sample), \
			path("${sample}.mapped.merged.bam"), emit: merged_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_MergeBamAlignment" : ""
		"""
		gatk MergeBamAlignment ${scratch_field} \
		--VALIDATION_STRINGENCY SILENT \
		--EXPECTED_ORIENTATIONS FR \
		--ATTRIBUTES_TO_RETAIN X0 \
		--ALIGNED_BAM ${mapped_bam} \
		--UNMAPPED_BAM ${unmapped_bam}  \
		--OUTPUT ${sample}.mapped.merged.bam \
		--REFERENCE_SEQUENCE ${ref} \
		--PAIRED_RUN true \
		--SORT_ORDER "unsorted" \
		--IS_BISULFITE_SEQUENCE false \
		--ALIGNED_READS_ONLY false \
		--CLIP_ADAPTERS false \
		--ADD_MATE_CIGAR true \
		--MAX_INSERTIONS_OR_DELETIONS -1 \
		--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
		--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
		--ALIGNER_PROPER_PAIR_FLAGS true \
		--UNMAP_CONTAMINANT_READS true
		"""
}









process MARKDUPLICATESSPARK {
	errorStrategy 'ignore'

	publishDir "${params.output}/mapping_stats", mode: 'copy', pattern: "marked_dup_metrics*"
	

	input:
		tuple val(sample), path(merged_bam)


	output:
		tuple \
			val(sample), \
			path("${sample}.dedupped.sorted.bam"), \
			path("${sample}.dedupped.sorted.bam.bai"), emit: deduppedsorted_bam
		tuple \
			val(sample), \
			path("marked_dup_metrics_${sample}.txt"), emit: dedupped_txt

	script:
		"""
		gatk MarkDuplicatesSpark \
		-I ${merged_bam} \
		-O ${sample}.dedupped.sorted.bam \
		-M marked_dup_metrics_${sample}.txt \
		--remove-all-duplicates false \
		--optical-duplicate-pixel-distance 2500 \
		--read-validation-stringency SILENT \
		--create-output-bam-index true

		chmod 777 \$(find . -user root) 
		"""


}








process SORTSAM {	
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(merged_bam)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.sorted.bam"), \
			path("${sample}.sorted.bai"), emit: sorted_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_SortSam" : ""
		"""
		gatk SortSam ${scratch_field} \
		--INPUT ${merged_bam} \
		--OUTPUT ${sample}.sorted.bam \
		--SORT_ORDER "coordinate" \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE false 
		"""
}








process SETTAGS {	
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(sorted_bam), path(sorted_bai)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.tag.bam"), \
			path("${sample}.tag.bai"), emit: tagged_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_SortSam" : ""
		"""
		gatk  SetNmMdAndUqTags ${scratch_field} \
		--INPUT ${sorted_bam} \
		--OUTPUT ${sample}.tag.bam \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE false \
		--REFERENCE_SEQUENCE ${ref}
		"""
}









process BASERECALIBRATOR {	
	publishDir "${params.output}/mapping_stats", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(deduppedsorted), path(deduppedsorted_bai)
		path ref
		path index
		path dict
		path reference_gzi
		path g1000_knownsites
		path g1000_knownsites_idx
		path mills_knownsites
		path mills_knownsites_idx
		path dbsnp_knownsites
		path dbsnp_knownsites_idx
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("before_recalibrated_bqsr_data_${sample}.recal.table"), emit: bqsr_table

	script:
		// def scratch_field = scratch ? "--TMP-DIR ${scratch}/${sample}_SortSam" : ""
		def g1000_knownsites_field = g1000_knownsites ? "--known-sites ${g1000_knownsites}" : ""
		def mills_knownsites_field = mills_knownsites ? "--known-sites ${mills_knownsites}" : ""
		def dbsnp_knownsites_field = dbsnp_knownsites ? "--known-sites ${dbsnp_knownsites}" : ""
		"""
		gatk BaseRecalibrator \
		-R ${ref} \
		-I ${deduppedsorted} \
		--use-original-qualities \
		${g1000_knownsites_field} \
		${mills_knownsites_field} \
		${dbsnp_knownsites_field} \
		-O before_recalibrated_bqsr_data_${sample}.recal.table
		"""
}









process APPLYBQSR {	
	publishDir "${params.output}/bams", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(deduppedsorted), path(deduppedsorted_bai), path(bqsr_table)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.bam"), \
			path("${sample}.bai"), emit: bam
		tuple \
			val(sample), \
			path("${sample}.bam.md5"), emit: md5 

	script:
		// def scratch_field = scratch ? "--TMP-DIR ${scratch}/${sample}_SortSam" : ""
		"""
		gatk ApplyBQSR \
		-R ${ref} \
		-I ${deduppedsorted}  \
		--bqsr ${bqsr_table} \
		-O ${sample}.bam \
		--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
		--add-output-sam-program-record \
		--create-output-bam-md5 \
		--use-original-qualities \
		--create-output-bam-index
		"""
}







process LOCALBAM {	

	input:
		path inputdir
		val sample2analyce
		
	output:
		tuple \
			val(sample2analyce_config), \
			path("${sample2analyce_config}.bam"), \
			path("${sample2analyce_config}.bai"), emit: bam

	script:
		sample2analyce_config = sample2analyce[0]
		"""
		ln -s ${inputdir}/${sample2analyce_config}*.bam ${sample2analyce_config}.bam
		ln -s ${inputdir}/${sample2analyce_config}*.bai ${sample2analyce_config}.bai
		"""
}





























process MOSDEPTH {
	publishDir "${params.output}/qc/mosdepth/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed

	output:
		path("${sample}_qc.*"), emit: mosdepth

	script:
		def bed_field = bed ? "--by ${bed}" : ""	
		
		"""
		mosdepth -x --no-per-base ${bed_field} ${sample}_qc ${bam}
		"""
}



process MOSDEPTH_JOIN {
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path ""
		val readthreshold
		val analysis
		val run

	output:
		path("minCovFilterResults_${analysis}.txt"), emit: coverage
		path("discarded_samples_${analysis}.txt"), emit: discarded_samples optional true 


	script:
		
		"""
		for file in *_qc.mosdepth.region.dist.txt; do

			sample="\$(basename \${file} _qc.mosdepth.region.dist.txt)"

			cov=\$(awk -v reads="${readthreshold}" '{if(\$1=="total" && \$2==reads){print \$3}}' \${file})
			pass=\$(awk 'BEGIN{if ('\${cov}'>'0.9') print 0}')

			if [ "\${pass}" = "0" ]; then
				echo -e \${sample}"\\t"${readthreshold}"\\t"\${cov}"\\tPASSED\\t"${run} >> minCovFilterResults_${analysis}.txt
			else
				echo -e \${sample}"\\t"${readthreshold}"\\t"\${cov}"\\tFAILED\\t"${run} >> minCovFilterResults_${analysis}.txt
				echo -e \${sample} >> discarded_samples_${analysis}.txt
			fi
		done
		"""
}



process MOSDEPTH_PLOT {
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path ""
		path tasks

	output:
		path("mosdepth.region.dist.html"), emit: plot
		

	script:
		
		"""
		python ${tasks}/plot-dist.py *.mosdepth.region.dist.txt -o mosdepth.region.dist.html
		"""
}





process MOSDEPTH_COV {
	publishDir "${params.output}/qc/mosdepth_cov", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val padding

	output:
		tuple \
			val(sample), \
			path("${sample}*bed"), emit: mosdepth

	script:
		
		if(bed)
			"""
			mosdepth --quantize 10: -n -x ${sample}_cov ${bam}
			zcat ${sample}_cov.quantized.bed.gz > ${sample}.global.quantized.bed

			awk -v pad="${padding}" '{print \$1"\\t"\$2-pad"\\t"\$3+pad"\\t"\$4}' ${bed} | bedtools merge -i - > bedpadding.bed
			bedtools intersect -a bedpadding.bed -b ${sample}_cov.quantized.bed.gz -wb | awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$NF}' > ${sample}.padding.quantized.bed

			
			#bedtools coverage -a BED -b BAM -mean -sorted > MeanCoverageBED.bedgraph
			#bedtools map -a BED -b BAM -mean -sorted > MeanCoverageBED.bedgraph
			"""

		else
			"""
			mosdepth --quantize 10: -n -x ${sample}_cov ${bam}
			zcat ${sample}_cov.quantized.bed.gz > ${sample}.global.quantized.bed
			"""
}








process GENOMECOV {
	publishDir "${params.output}/qc/genomecov/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed

	output:
		tuple \
			val(sample), \
			path("${sample}*bed"), emit: genomecov

	script:
		
		if(bed)
			"""
			samtools view -b -F 3332 ${bam} | bedtools genomecov -bga -ibam stdin > ${sample}.genomecov.bed
			#bedtools intersect -a ${sample}.genomecov.bed -b ${bed} > ${sample}.panelcov.bed
			bedtools intersect -a ${sample}.genomecov.bed -b ${bed} -wb > ${sample}.panelcov.bed
			"""

		else
			"""
			samtools view -b -F 3332 ${bam} | bedtools genomecov -bga -ibam stdin > ${sample}.genomecov.bed
			"""
}
// 3332
// read unmapped (0x4)
// not primary alignment (0x100)
// read is PCR or optical duplicate (0x400)
// supplementary alignment (0x800)

// 3076
// read unmapped (0x4)
// read is PCR or optical duplicate (0x400)
// supplementary alignment (0x800)





process SAMTOOLS_FLAGSTAT {
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.samtools_flagstat.txt"), emit: flagstat

	script:

		"""
		samtools flagstat ${bam} > ${sample}.samtools_flagstat.txt
		"""
}






process READ_LENGTH_STATS {
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path tasks

	output:
		tuple \
			val(sample), \
			path("${sample}.read_lenghts.txt"), emit: read_length_stats

	script:

		"""
		samtools view -F 3332 ${bam} | cut -f 10 | perl -ne 'chomp;print length(\$_) . "\\n"' > ${sample}.read_lenghts.txt
		"""
}



process SEQUENCING_QUALITY_SCORES {
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.quality.txt"), emit: quality

	script:

		"""
		samtools view -F 3332 ${bam} | cut -f 11 | grep -o . | awk '{c[\$0]++}END{for(l in c){print c[l], l}}' | sort -n > ${sample}.quality.txt
		"""
}



process SEQUENCING_CG_AT_CONTENT {
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.CG_AT.txt"), emit: cg_at

	script:

		"""
		samtools view -F 3332 ${bam} | cut -f 10 | grep -o . | awk '{c[\$0]++}END{for(l in c){print c[l], l}}' | sort -n > ${sample}.CG_AT.txt
		"""
}




process NREADS_NONDUP_UNIQ {
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.nreads_nondup_uniq.txt"), emit: nreads_nondup_uniq

	script:

		"""
		samtools view -F 3332 ${bam} | wc -l > ${sample}.nreads_nondup_uniq.txt
		"""
}



process QC_SUMMARY {
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(""), path("${sample}.samtools_flagstat.txt"), path("${sample}.read_lenghts.txt"), path("${sample}.quality.txt"), path("${sample}.CG_AT.txt"), path("${sample}.nreads_nondup_uniq.txt")
		path bed
		path tasks

	output:
		tuple \
			val(sample), \
			path("${sample}*txt"), emit: quality_summary

	script:
		def bed_path = bed ? "--bed_path ${bed} " : ''
		def panelcov_path = bed ? "--panelcov_path ${sample}.panelcov.bed " : ''
		def bedoutpath = bed ? "--bedout ${sample}.library.stats.txt " : ''

		"""
		nreads_nondup_uniq="\$(cat ${sample}.nreads_nondup_uniq.txt)"

		Rscript ${tasks}/quality_summary.R \\
		--samplename ${sample} \\
		--genomecov_path ${sample}.genomecov.bed \\
		--quality_path ${sample}.quality.txt \\
		--samtools_flagstat_path ${sample}.samtools_flagstat.txt \\
		--read_lenghts_path ${sample}.read_lenghts.txt \\
		--cg_at_path ${sample}.CG_AT.txt \\
		${panelcov_path}\\
		${bed_path}\\
		${bedoutpath}\\
		--output ${sample}.quality.summary.txt \\
		--nreads_nondup_uniq \${nreads_nondup_uniq}
		"""
}




process RUN_QC_CAT {
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path("")
		val runname

	output:
		tuple \
			val(runname), \
			path("${runname}.quality.summary.txt"), emit: quality_summary

	script:
		"""
		first_file="\$(ls *quality.summary.txt | head -n 1)"
		head -n 1 \${first_file} > ${runname}.quality.summary.txt.tmp
		tail -n 1 -q *quality.summary.txt >> ${runname}.quality.summary.txt.tmp
		mv ${runname}.quality.summary.txt.tmp ${runname}.quality.summary.txt
		"""
}



















process HAPLOTYPECALLER {	
	// publishDir "${params.output}/", mode: 'copy'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val intervals
		val padding
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.vcf"), \
			path("${sample}.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_HaplotypeCaller" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_HaplotypeCaller" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}

		gatk HaplotypeCaller ${scratch_field} \
		--java-options "-Xmx${params.bigmem}g" \
		-R ${ref} \
		-I ${bam} \
		-O ${sample}.vcf \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}
		"""
}






process SELECT_SNV {	

	input:
		tuple val(sample), path(vcf), path(idx)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.snp.vcf"), \
			path("${sample}.snp.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_selectsnvs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_selectsnvs" : ""

		"""
		${scratch_mkdir}

		gatk SelectVariants ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--select-type-to-include SNP \
		-O ${sample}.snp.vcf
		"""
}





process SELECT_INDEL {	

	input:
		tuple val(sample), path(vcf), path(idx)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.indel.vcf"), \
			path("${sample}.indel.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_selectindel" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_selectindel" : ""

		"""
		${scratch_mkdir}

		gatk SelectVariants ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--select-type-to-include INDEL \
		-O ${sample}.indel.vcf
		"""
}




process SELECT_MIX {	

	input:
		tuple val(sample), path(vcf), path(idx)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.mix.vcf"), \
			path("${sample}.mix.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_selectmix" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_selectmix" : ""

		"""
		${scratch_mkdir}

		gatk SelectVariants ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--select-type-to-include MIXED \
		--select-type-to-include MNP \
		--select-type-to-include SYMBOLIC \
		-O ${sample}.mix.vcf
		"""
}





process FILTRATION_SNV {	

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.filtered.snv.vcf"), \
			path("${sample}.filtered.snv.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filtersnvs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filtersnvs" : ""

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O ${sample}.filtered.snv.vcf
		"""
}







process FILTRATION_INDEL {	

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.filtered.indel.vcf"), \
			path("${sample}.filtered.indel.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filterindel" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filterindel" : ""

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-O ${sample}.filtered.indel.vcf
		"""
}








process FILTRATION_MIX {	

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.filtered.mix.vcf"), \
			path("${sample}.filtered.mix.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filtermix" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filtermix" : ""

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-O ${sample}.filtered.mix.vcf
		"""
}






process MERGE_VCF {	
	
	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(vcf_snv), path(idx_snv), path(vcf_indel), path(idx_indel), path(vcf_mix), path(idx_mix)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.gatkLabeled.vcf"), \
			path("${sample}.gatkLabeled.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--TMP_DIR ${scratch}/${sample}_filtermix" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filtermix" : ""

		"""
		${scratch_mkdir}

		gatk MergeVcfs ${scratch_field} \
		-R ${ref} \
		-I ${vcf_snv} \
		-I ${vcf_indel} \
		-I ${vcf_mix} \
		-O ${sample}.gatkLabeled.vcf
		"""
}



process FILTER_VCF {	

	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(vcf_snv), path(idx_snv)
		
	output:
		tuple \
			val(sample), \
			path("${sample}.final.vcf"), emit: vcf

	script:

		"""
		filter_vep \
		-i ${vcf_snv} -o ${sample}.final.vcf \
		--filter "(FILTER = PASS) and (DP > 10) and \
		(CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY)" \
		--force_overwrite
		"""
}









process LOCALVCF {	

	input:
		path inputdir
		val sample2analyce
		
	output:
		tuple \
			val(sample2analyce_config), \
			path("${sample2analyce_config}.vcf"), emit: vcf

	script:
		sample2analyce_config = sample2analyce[0]
		"""
		if [ -f ${inputdir}/${sample2analyce_config}*.vcf.gz ]; then
			zcat ${inputdir}/${sample2analyce_config}*.vcf.gz > ${sample2analyce_config}.vcf
		else
			cp ${inputdir}/${sample2analyce_config}*.vcf ${sample2analyce_config}.vcf
		fi
		"""
}







process FORMAT2INFO {
	
	input:
		tuple val(sample), path(final_vcf)

	output:
		tuple \
			val(sample), \
			path("${sample}.vcf_to_annotate.vcf.gz"), \
			path("${sample}.vcf_to_annotate.vcf.gz.tbi"), \
			path("${sample}.fields.txt"), emit: sample_info

	
	script:

		"""
		bcftools view -h ${final_vcf} | grep "##" > ${sample}.vcf_to_annotate.vcf
		echo "##INFO=<ID=variant_id,Number=.,Type=String,Description=\\"variant identification\\">" >> ${sample}.vcf_to_annotate.vcf
		echo "##INFO=<ID=Original_pos,Number=.,Type=String,Description=\\"original position\\">" >> ${sample}.vcf_to_annotate.vcf
		for sample in \$(bcftools query -l ${final_vcf})
		do
			echo "##INFO=<ID=\${sample}_GT,Number=.,Type=String,Description=\\"\${sample} Genotype\\">" >> ${sample}.vcf_to_annotate.vcf
			echo "##INFO=<ID=\${sample}_AD,Number=.,Type=String,Description=\\"\${sample} Allelic depths for the ref and alt alleles in the order listed\\">" >> ${sample}.vcf_to_annotate.vcf
			echo "##INFO=<ID=\${sample}_DP,Number=.,Type=String,Description=\\"\${sample} Approximate read depth (reads with MQ=255 or with bad mates are filtered)\\">" >> ${sample}.vcf_to_annotate.vcf
			echo "##INFO=<ID=\${sample}_GQ,Number=.,Type=String,Description=\\"\${sample} Genotype Quality\\">" >> ${sample}.vcf_to_annotate.vcf
		done
		bcftools view -h ${final_vcf} | grep "#CHROM" | cut -f1-8 >> ${sample}.vcf_to_annotate.vcf


		bcftools query -f 'variant_id=%CHROM\\_%POS\\_%REF\\_%ALT;Original_pos=%POS;[;%SAMPLE\\_GT=%GT][;%SAMPLE\\_AD=%AD][;%SAMPLE\\_DP=%DP][;%SAMPLE\\_GQ=%GQ]\\n' ${final_vcf} | sed 's/;//2' | sed 's/,/_/g' > new_info.txt
		bcftools view -H ${final_vcf} | cut -f1-8 > old_info.txt
		paste -d ';' old_info.txt new_info.txt >> ${sample}.vcf_to_annotate.vcf
		
		bgzip ${sample}.vcf_to_annotate.vcf
		tabix -p vcf ${sample}.vcf_to_annotate.vcf.gz


		if bcftools view -h ${final_vcf} | grep -Fq hiConfDeNovo; then
			fields=",hiConfDeNovo,loConfDeNovo,variant_id,Original_pos"
		else
			fields=",variant_id,Original_pos"
		fi
		 
		for sample in \$(bcftools query -l ${final_vcf}); do fields="\$(echo "\${fields},\${sample}_GT,\${sample}_AD,\${sample}_DP,\${sample}_GQ")"; done
		echo \${fields} > ${sample}.fields.txt
		"""
}












process VEP {
	publishDir "${params.output}/snvs/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path(dbscSNV)
		path(dbscSNV_tbi)
		path loFtool
		path exACpLI
		path(dbNSFP)
		path(dbNSFP_tbi)
		path maxEntScan
		path(cADD_INDELS)
		path(cADD_INDELS_tbi)
		path(cADD_SNVS)
		path(cADD_SNVS_tbi)
		path(kaviar)
		path(kaviar_tbi)
		path(cCRS_DB)
		path(cCRS_DB_tbi)
		path(dENOVO_DB)
		path(dENOVO_DB_tbi)
		path(cLINVAR)
		path(cLINVAR_tbi)
		path(gNOMADg)
		path(gNOMADg_tbi)
		path(gNOMADe)
		path(gNOMADe_tbi)
		path(gNOMADg_cov)
		path(gNOMADg_cov_tbi)
		path(gNOMADe_cov)
		path(gNOMADe_cov_tbi)
		path(cSVS)
		path(cSVS_tbi)
		path(mutScore)
		path(mutScore_tbi)
		path(mAF_FJD_COHORT)
		path(mAF_FJD_COHORT_tbi)
		path(spliceAI_SNV)
		path(spliceAI_SNV_tbi)
		path(spliceAI_INDEL)
		path(spliceAI_INDEL_tbi)
		val vep_threads
		path vep_cache
		path vep_plugins
		path vep_fasta
		path vep_fai
		path vep_gzi
		val vep_assembly
		tuple val(sample), path(final_vcf), path(sample_info), path(sample_info_index), path(sample_info_fields)
		
	
	output:
		tuple \
			val(sample), \
			path("${sample}.vep.tsv"), emit: vep_tsv
	
	script:
		def dbscSNV_config = dbscSNV ? "--plugin dbscSNV,${dbscSNV} " : ''
		def loFtool_config = loFtool ? "--plugin LoFtool,${loFtool} " : ''
		def exACpLI_config = exACpLI ? "--plugin ExACpLI,${exACpLI} " : ''
		def dbNSFP_config  = dbNSFP  ? "--plugin dbNSFP,${dbNSFP},\
LRT_pred,M-CAP_pred,MetaLR_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,\
FATHMM_pred,MetaRNN_pred,PrimateAI_pred,DEOGEN2_pred,BayesDel_addAF_pred,BayesDel_noAF_pred,ClinPred_pred,\
LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,\
phyloP30way_mammalian,phastCons30way_mammalian,GERP++_RS,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue " : ''
		def maxEntScan_config     = maxEntScan    ? "--plugin MaxEntScan,${maxEntScan} " : ''
		def cADD_config           = cADD_INDELS && cADD_SNVS ? "--plugin CADD,${cADD_INDELS},${cADD_SNVS} " : ''
		def kaviar_config         = kaviar         ? "--custom ${kaviar},kaviar,vcf,exact,0,AF,AC,AN " : ''
		def cCRS_DB_config        = cCRS_DB        ? "--custom ${cCRS_DB},gnomAD_exomes_CCR,bed,overlap,0 " : ''
		def dENOVO_DB_config      = dENOVO_DB      ? "--custom ${dENOVO_DB},denovoVariants,vcf,exact,0,SAMPLE_CT " : ''
		def cLINVAR_config        = cLINVAR        ? "--custom ${cLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN " : ''
		def gNOMADg_config        = gNOMADg        ? "--custom ${gNOMADg},gnomADg,vcf,exact,0,AF,AC,AN,nhomalt,popmax,AF_popmax,AC_popmax,AF_nfe,AC_nfe,filt " : ''
		def gNOMADe_config        = gNOMADe        ? "--custom ${gNOMADe},gnomADe,vcf,exact,0,AF,AC,AN,nhomalt,popmax,AF_popmax,AC_popmax,AF_nfe,AC_nfe,filt " : ''
		def gNOMADg_cov_config    = gNOMADg_cov    ? "--custom ${gNOMADg_cov},gnomADg_cov,vcf,overlap,0,median,perc_20x " : ''
		def gNOMADe_cov_config    = gNOMADe_cov    ? "--custom ${gNOMADe_cov},gnomADe_cov,vcf,overlap,0,median,perc_20x " : ''
		def cSVS_config           = cSVS           ? "--custom ${cSVS},CSVS,vcf,exact,0,AF,AC " : ''
		def mutScore_config       = mutScore       ? "--custom ${mutScore},Mut,vcf,exact,0,Score " : ''
		def mAF_FJD_COHORT_config = mAF_FJD_COHORT ? "--custom ${mAF_FJD_COHORT},FJD_MAF,vcf,exact,0,AF,AC " : ''
		def spliceAI_SNV_config   = spliceAI_SNV   ? "--custom ${spliceAI_SNV},SpliceAI_SNV,vcf,exact,0,SpliceAI " : ''
		def spliceAI_INDEL_config = spliceAI_INDEL ? "--custom ${spliceAI_INDEL},SpliceAI_INDEL,vcf,exact,0,SpliceAI " : ''
		def sample_info_config    = sample_info    ? "--custom ${sample_info},SAMPLE,vcf,exact,0\$(cat ${sample_info_fields}) " : ''

		"""
		vep \\
		--cache --offline --dir_cache ${vep_cache} --dir_plugins ${vep_plugins} \\
		--refseq --species homo_sapiens --assembly ${vep_assembly} --force_overwrite --use_transcript_ref \\
		--verbose --fork ${vep_threads} --tab --format vcf --no_stats \\
		--fasta ${vep_fasta} \\
		--input_file ${final_vcf} \\
		--output_file ${sample}.vep.tsv \\
		--check_existing --canonical --numbers --hgvs --biotype --regulatory --symbol --protein \\
		--sift p --polyphen p --allele_number --variant_class --pubmed \\
		${dbscSNV_config}\\
		${loFtool_config}\\
		${exACpLI_config}\\
		${dbNSFP_config}\\
		${maxEntScan_config}\\
		${cADD_config}\\
		${kaviar_config}\\
		${cCRS_DB_config}\\
		${dENOVO_DB_config}\\
		${cLINVAR_config}\\
		${gNOMADg_config}\\
		${gNOMADe_config}\\
		${gNOMADg_cov_config}\\
		${gNOMADe_cov_config}\\
		${cSVS_config}\\
		${mutScore_config}\\
		${mAF_FJD_COHORT_config}\\
		${spliceAI_SNV_config}\\
		${spliceAI_INDEL_config}\\
		${sample_info_config}

		"""
}








process AUTOMAP {

	publishDir "${params.output}/automap/", mode: 'copy'
	errorStrategy 'retry'
	
	input:
		tuple val(sample), path(final_vcf)
		val automap_assembly
		

	output:
	tuple \
		val(sample), \
		path("*HomRegions*"), emit: roh_automap 
	
	script:
		if ( task.attempt == 1 )
			"""
			if [[ \$(bcftools query -l ${final_vcf} | wc -l) -gt 1 ]]; then
				for sample in \$(bcftools query -l ${final_vcf}); do

					bcftools view -s \${sample} -O v -o \${sample}.indv.vcf ${final_vcf}

					if [[ \$(bcftools view -H \${sample}.indv.vcf | wc -l) -gt 10000 ]]; then
						
						/home/docker/AutoMap/AutoMap_v1.2.sh \\
						--vcf \${sample}.indv.vcf \\
						--out . \\
						--genome ${automap_assembly}

						mv \${sample}/* .

					else
						echo "Less than 10k variants" > \${sample}_no_automap_HomRegions.txt
					fi
				done
			
			else
				if [[ \$(bcftools view -H ${final_vcf} | wc -l) -gt 10000 ]]; then
					
					/home/docker/AutoMap/AutoMap_v1.2.sh \\
					--vcf ${final_vcf} \\
					--out . \\
					--genome ${automap_assembly}

					mv */* .
				
				else
					echo "Less than 10k variants" > ${sample}_no_automap_HomRegions.txt
				fi
			fi
			"""
		else
			"""
			echo "Less than 10k variants (with quality)" > ${sample}_no_automap_HomRegions.txt
			"""
}








process PVM {

	publishDir "${params.output}/snvs/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(vep_tsv), path(roh_automap)
		path dbNSFP_gene 
		path omim 
		path regiondict 
		val maf 
		path genefilter 
		path glowgenes 
		path tasks 

	output:
		tuple \
			val(sample), \
			path("${sample}.pvm.tsv"), \
			path("${sample}.pvm.tsv.xlsx"), emit: pvm_tsv 
	
	script:
	
		def omim_field       = omim       ? "--omim ${omim} " : ''
		def genefilter_field = genefilter ? "--genefilter ${genefilter} " : ''
		def glowgenes_field  = glowgenes  ? "--glowgenes ${glowgenes} " : ''

		"""
		header_row="\$(head -n 1000 ${vep_tsv} | grep "#Uploaded_variation" -n | sed 's/:.*//')"

		Rscript ${tasks}/post-VEP_modification.R \\
		--input ${vep_tsv} \\
		--output ${sample}.pvm.tsv \\
		--numheader \${header_row} \\
		--dbNSFPgene ${dbNSFP_gene} \\
		--regiondict ${regiondict} \\
		--automap ./ \\
		--maf ${maf} \\
		${omim_field}\\
		${genefilter_field}\\
		${glowgenes_field}
		"""
}
























process BEDPROCCESING {

	publishDir "${params.output}/cnvs/", mode: 'copy'
	
	input:
		path bed
		val min_target 
		val window
		val chrx


	output:
		path("*cnv.bed"), emit: bed 
	
	script:
		if ( chrx && window )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && (\$1=="chrX" || \$1=="X") ){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.chrx.cnv.bed

			python $tasksPath/CNV_windowSize.py \${panel}.min${min_target}bp.chrx.cnv.bed \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted
			sort -V -k1,1 -k2,2 \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted | uniq > \${panel}.window125bp.min${min_target}bp.chrx.cnv.bed

			rm \${panel}.min${min_target}bp.chrx.cnv.bed \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted
			"""

		else if ( chrx )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && (\$1=="chrX" || \$1=="X") ){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.chrx.cnv.bed
			"""

		else if ( window )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && \$1!="chrX" && \$1!="chrY" && \$1!="Y" && \$1!="X"){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.cnv.bed

			python $tasksPath/CNV_windowSize.py \${panel}.min${min_target}bp.cnv.bed \${panel}.min${min_target}bp.cnv.bed_unsorted
			sort -V -k1,1 -k2,2 \${panel}.min${min_target}bp.cnv.bed_unsorted | uniq > \${panel}.window125bp.min${min_target}bp.cnv.bed

			rm \${panel}.min${min_target}bp.cnv.bed \${panel}.min${min_target}bp.cnv.bed_unsorted
			"""

		else
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && \$1!="chrX" && \$1!="chrY" && \$1!="Y" && \$1!="X"){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.cnv.bed
			"""
}







process EXOMEDEPTH {

	publishDir "${params.output}/cnvs/exomedepth", mode: 'copy'
	errorStrategy 'retry'

	input:
		path("")
		path("") 
		path bed
		val runname 
		path tasks

	output:
		tuple \
			val(runname), \
			path("exomedepth*"), emit: cnvs 

		tuple \
			val(runname), \
			path("exomedepth.toAnnotate.txt"), emit: toannotate
	
	script:
		if ( task.attempt == 1 )
			"""
			Rscript ${tasks}/exomeDepth.R -d . -o . -b ${bed} -n ${runname}
			"""
		else
			"""
			echo "No ExomeDepth" > exomedepth.txt
			"""
}		








process CONVADING {

	publishDir "${params.output}/cnvs/convading", mode: 'copy'
	errorStrategy 'retry'

	input:
		path("")
		path("") 
		path bed
		val runname 
		path tasks
		path fai


	output:
		tuple \
			val(runname), \
			path("CoNVaDING*"), emit: cnvs 

		tuple \
			val(runname), \
			path("CoNVaDING.toAnnotate.txt"), emit: toannotate
	
	script:
		if ( task.attempt == 1 )
			"""
			python ${tasks}/CoNVading_pipeline.py . ${bed} ./ ${runname} ${tasks} ${fai}
			"""
		else
			"""
			echo "No CoNVaDING" > CoNVaDING.txt
			"""
}





process PANELCNMOPS {

	publishDir "${params.output}/cnvs/panelmops", mode: 'copy'
	errorStrategy 'retry'

	input:
		path("")
		path("") 
		path bed
		val runname 
		path tasks

	output:
		tuple \
			val(runname), \
			path("panelcn.MOPS*"), emit: cnvs 
		
		tuple \
			val(runname), \
			path("panelcn.MOPS.toAnnotate.txt"), emit: toannotate
	
	script:
		if ( task.attempt == 1 )
			"""
			Rscript ${tasks}/panelcnMops.R -d . -o . -b ${bed} -n ${runname}
			"""
		else
			"""
			echo "No panelcn.MOPS" > panelcn.MOPS.txt
			"""
}





process CNV_RESULT_MIXER {

	publishDir "${params.output}/cnvs", mode: 'copy'
	
	input:
		tuple val(runname), path("exomedepth.toAnnotate.txt"), path("CoNVaDING.toAnnotate.txt"), path("panelcn.MOPS.toAnnotate.txt")
		path tasks
		path samples2analyce


	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.merged.bed"), emit: merged_bed 

		tuple \
			val(runname), \
			path("colnames.txt"), emit: colnames
	
	script:
		def samples2analyce_field   = samples2analyce ? "--samples ${samples2analyce} " : ""
		"""
		Rscript ${tasks}/CNV_result_mixer2.R \
		--inputdir . \
		--outputfile ${runname}.CNV.merged.bed \
		${samples2analyce_field}

		"""
}






process ANNOTSV {

	// publishDir "${params.output}/cnvs", mode: 'copy'
	
	input:
		tuple val(runname), path(merged_cnv)
		path genelist
		val annotsv_assembly
		path annotsv_path


	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.annotated.tsv"), emit: annotated_cnv 
	
	script:
		// def genelist_field   = genelist ? "-candidateGenesFile ${genelist} -candidateGenesFiltering " : ""

		"""
		export ANNOTSV=${annotsv_path}

		${annotsv_path}/bin/AnnotSV \
		-SVinputFile ${merged_cnv} \
		-outputFile ${runname}.CNV.annotated.tsv \
		-annotationMode both \
		-genomeBuild ${annotsv_assembly} \
		-svtBEDcol 4 \
		-samplesidBEDcol 5 \
		-SVminSize 20 
		
		mv *_AnnotSV/${runname}.CNV.annotated.tsv .
		"""
}





process PAM {

	publishDir "${params.output}/cnvs", mode: 'copy'
	
	input:
		tuple val(runname), path(annotated_cnv)
		tuple val(runname), path(colnames)
		path genefilter 
		path glowgenes
		path annotsv_path


	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.annotated.final.tsv"), emit: annotated_cnv 
	
	script:
		def genefilter_field = genefilter ? "--genefilter ${genefilter} " : ''
		def glowgenes_field  = glowgenes  ? "--glowgenes ${glowgenes} " : ''

		"""
		Rscript ${tasks}/post-AnnotSV_modification.R \\
		--input ${annotated_cnv} \\
		--outputfile ${runname}.CNV.annotated.final.tsv \\
		--extracolnames ${colnames} \\
		${genefilter_field} \\
		${glowgenes_field}
		"""
}




















process GVCF_HAPLOTYPECALLER {	

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val intervals
		val padding
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.g.vcf"), \
			path("${sample}.g.vcf.idx"), emit: gvcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_HaplotypeCaller" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_HaplotypeCaller" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}

		gatk HaplotypeCaller ${scratch_field} \
		--java-options "-Xmx${params.bigmem}g" \
		-R ${ref} \
		-I ${bam} \
		-ERC GVCF \
		-O ${sample}.g.vcf \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality \
		-A MappingQualityRankSumTest -A ReadPosRankSumTest \
		-A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}
		"""
}






process COMBINE_GVCF {	

	input:
		path("")
		path("")
		val runname
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("${runname}.g.vcf"), \
			path("${runname}.g.vcf.idx"), emit: gvcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_CombineGVCFs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_CombineGVCFs" : ""

		"""
		${scratch_mkdir}

		gatk CombineGVCFs ${scratch_field} \
		-R ${ref} \
		--variant *.g.vcf \
		-O ${runname}.g.vcf
		"""
}





process GENOTYPE_GVCF {	

	input:
		tuple val(runname), path(gvcf), path(idx)
		path ped
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("${runname}.vcf"), \
			path("${runname}.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_GenotypeGVCFs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_GenotypeGVCFs" : ""
		def ped_field = ped ? "-ped ${ped}" : ""

		"""
		${scratch_mkdir}

		gatk GenotypeGVCFs ${scratch_field} \
		-R ${ref} \
		-V ${gvcf} \
		-O ${runname}.vcf \
		${ped_field}
		"""
}




process CALCULATE_GENOTYPE_POSTERIORS {	

	input:
		tuple val(runname), path(vcf), path(idx)
		path ped
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("filtered_INDEL_SNP_data_GP_${runname}.vcf"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_CalculateGenotypePosteriors" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_CalculateGenotypePosteriors" : ""
		"""
		${scratch_mkdir}

		gatk CalculateGenotypePosteriors ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		-ped ${ped} \
		-O filtered_INDEL_SNP_data_GP_${runname}.vcf \
		--skip-population-priors 
		"""
}







process VARIANT_FILTRATION {	

	input:
		tuple val(runname), path(vcf)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("filtered_INDEL_SNP_data_GP_GQfiltered_${runname}.vcf"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_VariantFiltration" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_VariantFiltration" : ""
		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--filter-expression "GQ < 20" \
		--filter-name "lowGQ" \
		-O filtered_INDEL_SNP_data_GP_GQfiltered_${runname}.vcf
		"""
}







process VARIANT_ANNOTATOR {	

	input:
		tuple val(runname), path(vcf)
		path ped
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("${runname}.va.vcf"), 
			path("${runname}.va.vcf.idx"),emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_CalculateGenotypePosteriors" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_CalculateGenotypePosteriors" : ""
		"""
		${scratch_mkdir}

		gatk VariantAnnotator ${scratch_field} \
		-R ${ref} \
		-V ${vcf}  \
		-O ${runname}.va.vcf \
		-A PossibleDeNovo  \
		-A StrandBiasBySample  \
		-A AS_FisherStrand \
		-ped ${ped} 
		"""
}

