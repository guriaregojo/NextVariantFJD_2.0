#!/usr/bin/env nextflow
import java.util.regex.Matcher

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl=2



/*
 * Print the default parameters
 */ 

log.info """\
FJD-PIPELINE    -    NEXTFLOW  
=============================
samples      : $params.samples
"""


/* 
 * Import modules 
 */
include { BS_CHECK } from './modules/execution_modules'
include { LOCAL_CHECK } from './modules/execution_modules'
include { BS_COPY } from './modules/execution_modules'
include { FASTQ_CONCATENATION } from './modules/execution_modules'
include { BWA } from './modules/execution_modules'
include { FASTQTOSAM } from './modules/execution_modules'
include { MERGEBAMALIGNMENT } from './modules/execution_modules'
include { MARKDUPLICATESSPARK } from './modules/execution_modules'
include { SORTSAM } from './modules/execution_modules'
include { SETTAGS } from './modules/execution_modules'
include { BASERECALIBRATOR } from './modules/execution_modules'
include { APPLYBQSR } from './modules/execution_modules'
include { LOCALBAM } from './modules/execution_modules'

include { HAPLOTYPECALLER } from './modules/execution_modules'
include { SELECT_SNV } from './modules/execution_modules'
include { SELECT_INDEL } from './modules/execution_modules'
include { SELECT_MIX } from './modules/execution_modules'
include { FILTRATION_SNV } from './modules/execution_modules'
include { FILTRATION_INDEL } from './modules/execution_modules'
include { FILTRATION_MIX } from './modules/execution_modules'
include { MERGE_VCF } from './modules/execution_modules'
include { FILTER_VCF } from './modules/execution_modules'

include { GVCF_HAPLOTYPECALLER } from './modules/execution_modules'
include { COMBINE_GVCF } from './modules/execution_modules'
include { GENOTYPE_GVCF } from './modules/execution_modules'
include { CALCULATE_GENOTYPE_POSTERIORS } from './modules/execution_modules'
include { VARIANT_FILTRATION } from './modules/execution_modules'
include { VARIANT_ANNOTATOR } from './modules/execution_modules'

include { LOCALVCF } from './modules/execution_modules'
include { FORMAT2INFO } from './modules/execution_modules'
include { AUTOMAP } from './modules/execution_modules'
include { VEP } from './modules/execution_modules'
include { PVM } from './modules/execution_modules'

include { MOSDEPTH } from './modules/execution_modules'
include { MOSDEPTH_JOIN as MOSDEPTH_JOIN_SNV } from './modules/execution_modules'
include { MOSDEPTH_JOIN as MOSDEPTH_JOIN_CNV } from './modules/execution_modules'
include { MOSDEPTH_PLOT } from './modules/execution_modules'
include { MOSDEPTH_COV } from './modules/execution_modules'
include { GENOMECOV } from './modules/execution_modules'
include { SAMTOOLS_FLAGSTAT } from './modules/execution_modules'
include { READ_LENGTH_STATS } from './modules/execution_modules'
include { SEQUENCING_QUALITY_SCORES } from './modules/execution_modules'
include { SEQUENCING_CG_AT_CONTENT } from './modules/execution_modules'
include { NREADS_NONDUP_UNIQ } from './modules/execution_modules'
include { QC_SUMMARY } from './modules/execution_modules'
include { RUN_QC_CAT } from './modules/execution_modules'

include { BEDPROCCESING } from './modules/execution_modules'
include { EXOMEDEPTH } from './modules/execution_modules'
include { CONVADING } from './modules/execution_modules'
include { PANELCNMOPS } from './modules/execution_modules'
include { CNV_RESULT_MIXER } from './modules/execution_modules'
include { ANNOTSV } from './modules/execution_modules'
include { PAM } from './modules/execution_modules'


// def final_vcf  = Channel.fromPath(params.final_vcf)
// def omim       = params.omim       ? Channel.fromPath(params.omim)       : Channel.empty()
// def genefilter = params.genefilter ? Channel.fromPath(params.genefilter) : Channel.empty()
// def glowgenes  = params.glowgenes  ? Channel.fromPath(params.glowgenes)  : Channel.empty()

// println "final_vcf: ${final_vcf}"
// // println "${omim}"

// checkPathParamList = [
//     params.dbNSFP_gene, params.omim,
//     params.genefilter, params.glowgenes
// ]
// // for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// // Check mandatory parameters
// if (params.dbNSFP_gene) { ch_input = Channel.fromPath(params.dbNSFP_gene) } else { exit 1, 'dbNSFP_gene samplesheet not specified!' }
// if (genefilter) { println "HOLAAAAAAAA" } else { exit 1, 'genefilter samplesheet not specified!' }
// println "dbNSFP_gene: ${params.dbNSFP_gene}"
// println "omim: ${omim}"
// println omim.view()
// println "genefilter: ${genefilter}"
// println "glowgenes: ${params.glowgenes}"

// if (params.genefilter) { println "ADIOS" } else { exit 1, 'genefilter samplesheet not specified!' }

// println Channel.fromPath(params.dbNSFP_gene).type()
// println file(params.omim).type()

// exit 1

// aa = params.analysis.toUpperCase().split(",")

// println aa


// if(aa.contains("M"))  {println "HOLA"} else {println "ADIOS"}
// if(params.analysis.toUpperCase().contains("M"))  {println "1"} else {println "2"}

// channel.value(params.analysis).view()





  
workflow CHECK_PARAMS {

	take:


	main:
		println "Checking..."

		//Check the existance of reference genome
		if ( !params.reference_fasta ) {exit 1, "Error: Missing reference fasta file definition.\n"} 
		else {file(params.reference_fasta, type: "file", checkIfExists: true)}

		if ( !params.reference_index ) {exit 1, "Error: Missing index reference fasta file definition.\n"} 
		else {file(params.reference_index, type: "file", checkIfExists: true)}

		if ( !params.reference_dict  ) {exit 1, "Error: Missing dictionary reference fasta file definition.\n"} 
		else {file(params.reference_dict, type: "file", checkIfExists: true)}

		if ( !params.reference_gzi   ) {exit 1, "Error: Missing gzi reference fasta file definition.\n"} 
		else {file(params.reference_gzi, type: "file", checkIfExists: true)}
		println "Referencie file check"



		// Check the existance of the input file
		if ( !params.input ) {exit 1, "Error: Missing input file definition.\n"}
		if ( params.input && !params.basespace ) { file(params.input, type: "dir", checkIfExists: true) }
		println "Input folder check"
		


		// Define the run name
		if ( params.basespace )  { runname = params.input }
		else if (params.runname) { runname = params.runname } 
		else { runname = new Date().format("yyyy-MM-dd_HH-mm") }
		println "Run name: $runname" 



		// Check the existance of the output directory
		if (params.output) { 
			file(params.output, type: "dir").mkdir()
			outputdir = file(params.output, type: "dir", checkIfExists: true) 
		} else {
			exit 1, "Error: Missing output directory definition.\n"
		}
		println "outputdir check: $outputdir"



		// Check the existance of the sample file
		if ( params.samples ) { file(params.samples, type: "file", checkIfExists: true) }
		println "Sample file check"



		// Check the existance of the bed file
		m = params.analysis.toUpperCase() =~ /[CX]/
		assert m instanceof Matcher
		if ( !params.bed && ( m || params.intervals ) && params.capture != "G") {
			exit 1, "Error: No bedfile provided and CNV calling or '--intervals' was provided for Panels or WES analisis.\n"
		} else if ( params.bed ) {
			bedfile = file(params.bed, type: "file", checkIfExists: true)
			bedfilecontent = bedfile.readLines()[0].split("\t").length
			if ( bedfilecontent != 4 ) {exit 1, "Error: Bedfile must have 4 columns.\n"}
		}
		println "Bed file check"



		// Check the existance of the ped (pedigree) file
		if ( params.ped &&  params.analysis.toUpperCase().contains("G") ) { pedfile = file(params.ped, type: "file", checkIfExists: true) }
		if ( params.ped && !params.analysis.toUpperCase().contains("G") ) { println "WARN: Ped (pedigree) file provided, but not used" }
		println "Ped file check"



		// Check the type of the analysis
		m = params.analysis.toUpperCase() =~ /[MQSGCXAN]/
		assert m instanceof Matcher
		if ( !m ) {exit 1, "Error: Cannot recognice the any of the specified analisis analysis.\nThe available analysis are: M (Mapping), S (SNV individual), G (SNV GVCF),\nA (SNV annotation), C (CNV calling), N (CNV annotation), X (chrX CNV calling)\n"}
		println "Analysis type check"



		// BaseSpace must start with mapping
		if ( !params.analysis.toUpperCase().contains("M") && params.basespace ) {exit 1, "Error: If basespace parameter is specify, the mapping (M) analysis must be specify.\n"}
		println "Incompatibility check"



		if ( params.basespace ) {

			BS_CHECK(
				params.input,
				params.baseuser,
				params.samples,
				params.analysis.toUpperCase() )
			controlsamples  = BS_CHECK.out.controlsamples
			samples2analyce = BS_CHECK.out.samples2analyce

		} else {

			LOCAL_CHECK(
				params.input,
				params.samples,
				params.analysis.toUpperCase() )
			controlsamples  = LOCAL_CHECK.out.controlsamples
			samples2analyce = LOCAL_CHECK.out.samples2analyce
		}

		

	emit:
		controlsamples_file = controlsamples
		samples2analyce_file = samples2analyce
		controlsamples  = controlsamples.splitCsv()
		samples2analyce = samples2analyce.splitCsv()
		runname = runname
}




workflow MAPPING {
	take:
		mapping_sample

	main:
		
		if ( params.basespace ) {

			BS_COPY (
				params.input,
				mapping_sample,
				params.baseuser)

			fastq = BS_COPY.out.fastq

		} else {
			println  mapping_sample
			FASTQ_CONCATENATION (
					params.input,
					mapping_sample )

			fastq = FASTQ_CONCATENATION.out.fastq
		}

		fastq.view()

		BWA (
			fastq,
			params.reference_fasta,
			params.bwa_amb,
			params.bwa_ann,
			params.bwa_pac,
			params.bwa_bwt,
			params.bwa_sa,
			params.bwa_threads )

		FASTQTOSAM (
			fastq,
			params.scratch )

		
		MERGEBAMALIGNMENT (
			BWA.out.mapped_bam.join(FASTQTOSAM.out.unmapped_bam),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )
		
		if ( !params.ignoreduplicates ) {
			
			MARKDUPLICATESSPARK (
				MERGEBAMALIGNMENT.out.merged_bam )

			sorted_bam = MARKDUPLICATESSPARK.out.deduppedsorted_bam

		} else {

			SORTSAM (
				MERGEBAMALIGNMENT.out.merged_bam,
				params.scratch )

			sorted_bam = SORTSAM.out.sorted_bam
		}


		SETTAGS (
			sorted_bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )


		BASERECALIBRATOR (
			SETTAGS.out.tagged_bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.g1000_knownsites,
			params.g1000_knownsites_idx,
			params.mills_knownsites,
			params.mills_knownsites_idx,
			params.dbsnp_knownsites,
			params.dbsnp_knownsites_idx,
			params.scratch )


		APPLYBQSR (
			SETTAGS.out.tagged_bam.join(BASERECALIBRATOR.out.bqsr_table),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )


	emit:
		bam = APPLYBQSR.out.bam

}






workflow SNVCALLING {
	take:
		bam
	main:

		HAPLOTYPECALLER (
			bam,
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_SNV (
			HAPLOTYPECALLER.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_INDEL (
			HAPLOTYPECALLER.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_MIX (
			HAPLOTYPECALLER.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_SNV (
			SELECT_SNV.out.vcf,
			params.scratch )

		FILTRATION_INDEL (
			SELECT_INDEL.out.vcf,
			params.scratch )

		FILTRATION_MIX (
			SELECT_MIX.out.vcf,
			params.scratch )

		MERGE_VCF (
			FILTRATION_SNV.out.vcf.join(FILTRATION_INDEL.out.vcf).join(FILTRATION_MIX.out.vcf),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTER_VCF (
			MERGE_VCF.out.vcf )

	emit:
		finalvcf = FILTER_VCF.out.vcf
}









workflow ANNOTATION {

	take:
		final_vcf

	main:

		FORMAT2INFO( 
			final_vcf )


		AUTOMAP(
			final_vcf,
			params.automap_assembly )


		VEP(
			params.dbscSNV,
			params.dbscSNV_tbi,
			params.loFtool,
			params.exACpLI,
			params.dbNSFP,
			params.dbNSFP_tbi,
			params.maxEntScan,
			params.cADD_INDELS,
			params.cADD_INDELS_tbi,
			params.cADD_SNVS,
			params.cADD_SNVS_tbi,
			params.kaviar,
			params.kaviar_tbi,
			params.cCRS_DB,
			params.cCRS_DB_tbi,
			params.dENOVO_DB,
			params.dENOVO_DB_tbi,
			params.cLINVAR,
			params.cLINVAR_tbi,
			params.gNOMADg,
			params.gNOMADg_tbi,
			params.gNOMADe,
			params.gNOMADe_tbi,
			params.gNOMADg_cov,
			params.gNOMADg_cov_tbi,
			params.gNOMADe_cov,
			params.gNOMADe_cov_tbi,
			params.cSVS,
			params.cSVS_tbi,
			params.mutScore,
			params.mutScore_tbi,
			params.mAF_FJD_COHORT,
			params.mAF_FJD_COHORT_tbi,
			params.spliceAI_SNV,
			params.spliceAI_SNV_tbi,
			params.spliceAI_INDEL,
			params.spliceAI_INDEL_tbi,
			params.vep_threads,
			params.vep_cache,
			params.vep_plugins,
			params.vep_fasta,
			params.vep_fai,
			params.vep_gzi,
			params.vep_assembly,
			final_vcf.join(FORMAT2INFO.out.sample_info) )



		PVM(
			VEP.out.vep_tsv.join(AUTOMAP.out.roh_automap),
			params.dbNSFP_gene,
			params.omim,
			params.regiondict,
			params.maf,
			params.genelist,
			params.glowgenes,
			params.tasks )


	// emit:


}







workflow CNVCALLING {
	take:
		bam
		bai
		runname
		samples2analyce

	main:

		// Bed file transformation
		BEDPROCCESING(
			params.bed,
			params.min_target,
			params.window,
			params.cnv_chrx )


		// CNV calling
		EXOMEDEPTH(
			bam,
			bai,
			BEDPROCCESING.out.bed,
			runname,
			params.tasks )

		CONVADING(
			bam,
			bai,
			BEDPROCCESING.out.bed,
			runname,
			params.tasks,
			params.fai_convading )
		
		PANELCNMOPS(
			bam,
			bai,
			BEDPROCCESING.out.bed,
			runname,
			params.tasks )


		// CNV merge
		CNV_RESULT_MIXER(
			EXOMEDEPTH.out.toannotate.join(CONVADING.out.toannotate).join(PANELCNMOPS.out.toannotate),
			params.tasks,
			samples2analyce )


		// CNV ANNOTATION
		ANNOTSV ( 
			CNV_RESULT_MIXER.out.merged_bed,
			params.genelist,
			params.annotsv_assembly,
			params.annotsv_path )

		PAM (
			ANNOTSV.out.annotated_cnv,
			CNV_RESULT_MIXER.out.colnames,
			params.genelist,
			params.glowgenes,
			params.tasks)

	// emit:

}









workflow COMBINEDSNVCALLING {
	take:
		bam
		runname
	main:

		GVCF_HAPLOTYPECALLER (
			bam,
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		vcf = GVCF_HAPLOTYPECALLER.out.gvcf.collect().flatten().filter( ~/.*g.vcf$/ ).toList()
		idx = GVCF_HAPLOTYPECALLER.out.gvcf.collect().flatten().filter( ~/.*g.vcf.idx$/ ).toList()

		COMBINE_GVCF(
			vcf,
			idx,
			runname,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		GENOTYPE_GVCF(
			COMBINE_GVCF.out.gvcf,
			params.ped,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )




		SELECT_SNV (
			GENOTYPE_GVCF.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_INDEL (
			GENOTYPE_GVCF.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_MIX (
			GENOTYPE_GVCF.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_SNV (
			SELECT_SNV.out.vcf,
			params.scratch )

		FILTRATION_INDEL (
			SELECT_INDEL.out.vcf,
			params.scratch )

		FILTRATION_MIX (
			SELECT_MIX.out.vcf,
			params.scratch )

		MERGE_VCF (
			FILTRATION_SNV.out.vcf.join(FILTRATION_INDEL.out.vcf).join(FILTRATION_MIX.out.vcf),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )


		if (params.ped){

			CALCULATE_GENOTYPE_POSTERIORS(
				MERGE_VCF.out.vcf,
				params.ped,
				params.reference_fasta,
				params.reference_index,
				params.reference_dict,
				params.reference_gzi,
				params.scratch )

			VARIANT_FILTRATION(
				CALCULATE_GENOTYPE_POSTERIORS.out.vcf,
				params.reference_fasta,
				params.reference_index,
				params.reference_dict,
				params.reference_gzi,
				params.scratch )

			VARIANT_ANNOTATOR(
				VARIANT_FILTRATION.out.vcf,
				params.ped,
				params.reference_fasta,
				params.reference_index,
				params.reference_dict,
				params.reference_gzi,
				params.scratch )

			vcf2filt = VARIANT_ANNOTATOR.out.vcf
			
		} else {
				
			vcf2filt = MERGE_VCF.out.vcf
		}


		FILTER_VCF (
			vcf2filt )

	emit:
		finalvcf = FILTER_VCF.out.vcf
}










workflow QUALITYCHECK {
	take:
		bam
		runname
	
	main:

		MOSDEPTH(
			bam,
			params.bed )

		MOSDEPTH_PLOT(
			MOSDEPTH.out.mosdepth.collect(),
			params.tasks )

		MOSDEPTH_COV(
			bam,
			params.bed,
			params.padding )



		// // SNV quality check
		// // if ( params.analysis.toUpperCase().contains("S") ) {
		// if ( params.analysis.toUpperCase().contains("Q") ) {

		// 	MOSDEPTH_JOIN_SNV(
		// 		MOSDEPTH.out.mosdepth.collect(),
		// 		params.qsnvthreshold,
		// 		"snv",
		// 		runname )
		// }



		// // CNV quality check
		// // if ( params.analysis.toUpperCase().contains("C") ) {
		// if ( params.analysis.toUpperCase().contains("Q") ) {
			
		// 	MOSDEPTH_JOIN_CNV(
		// 		MOSDEPTH.out.mosdepth.collect(),
		// 		params.qcnvthreshold,
		// 		"cnv",
		// 		runname )
		// }




		GENOMECOV(
			bam,
			params.bed )

		SAMTOOLS_FLAGSTAT(
			bam )

		READ_LENGTH_STATS(
			bam,
			params.tasks )

		SEQUENCING_QUALITY_SCORES(
			bam )

		SEQUENCING_CG_AT_CONTENT(
			bam )

		NREADS_NONDUP_UNIQ(
			bam )

		QC_SUMMARY(
			GENOMECOV.out.genomecov.join(SAMTOOLS_FLAGSTAT.out.flagstat).join(READ_LENGTH_STATS.out.read_length_stats).join(SEQUENCING_QUALITY_SCORES.out.quality).join(SEQUENCING_CG_AT_CONTENT.out.cg_at).join(NREADS_NONDUP_UNIQ.out.nreads_nondup_uniq),
			params.bed,
			params.tasks	)


		// library_stats = QC_SUMMARY.out.quality_summary.collect().flatten().filter( ~/.*library.stats.txt$/ ).toList()
		quality_summary = QC_SUMMARY.out.quality_summary.collect().flatten().filter( ~/.*quality.summary.txt$/ ).toList()

		RUN_QC_CAT(
			quality_summary,
			runname )
	// emit:
}























workflow {
	

	// Check that all parameters are OK
	CHECK_PARAMS ()
	



	// Mapping
	if ( params.analysis.toUpperCase().contains("M") ) {
		
		MAPPING( CHECK_PARAMS.out.controlsamples )
	} 


	


	// Quality check
	if ( params.analysis.toUpperCase().contains("Q") ) {
		
		if ( params.analysis.toUpperCase().contains("M") ) {
				
			bam = MAPPING.out.bam

		} else {
			
			LOCALBAM (
				params.input,
				CHECK_PARAMS.out.samples2analyce )

			bam = LOCALBAM.out.bam
		}

		QUALITYCHECK(
			bam,
			CHECK_PARAMS.out.runname)
	}






	// SNV calling
	if ( params.analysis.toUpperCase().contains("S") ) {
		if ( params.analysis.toUpperCase().contains("M") ) {
			
			bam = MAPPING.out.bam

		} else {
			
			LOCALBAM (
				params.input,
				CHECK_PARAMS.out.samples2analyce )

			bam = LOCALBAM.out.bam
		}


		// SNV calling
		SNVCALLING ( bam )
	}




	// SNV calling GVCF mode
	if ( params.analysis.toUpperCase().contains("G") ) {
		if ( params.analysis.toUpperCase().contains("M") ) {
			
			bam = MAPPING.out.bam

		} else {
			
			LOCALBAM (
				params.input,
				CHECK_PARAMS.out.samples2analyce )

			bam = LOCALBAM.out.bam
		}


		// SNV calling
		COMBINEDSNVCALLING ( 
			bam,
			CHECK_PARAMS.out.runname )
	}




	// SNV annotation
	if ( params.analysis.toUpperCase().contains("A") ) {
		if ( params.analysis.toUpperCase().contains("S") ) {
			
			vcf = SNVCALLING.out.finalvcf

		} else {
			
			LOCALVCF (
				params.input,
				CHECK_PARAMS.out.samples2analyce )

			vcf = LOCALVCF.out.vcf
		}


		ANNOTATION ( vcf )
	}







	// CNV calling
	if ( params.analysis.toUpperCase().contains("C") ) {
		if ( params.analysis.toUpperCase().contains("M") ) {
			
			bam = MAPPING.out.bam.collect().flatten().filter( ~/.*bam$/ ).toList()
			bai = MAPPING.out.bam.collect().flatten().filter( ~/.*bai$/ ).toList()
			quality_bam = MAPPING.out.bam

		} else {

			LOCALBAM (
				params.input,
				CHECK_PARAMS.out.controlsamples )
			bam = LOCALBAM.out.bam.collect().flatten().filter( ~/.*bam$/ ).toList()
			bai = LOCALBAM.out.bam.collect().flatten().filter( ~/.*bai$/ ).toList()

		}


		// CNV calling and annotation
		CNVCALLING ( 
			bam,
			bai, 
			CHECK_PARAMS.out.runname,
			CHECK_PARAMS.out.samples2analyce_file )
	}





}


















// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/13-1509.final.vcf --sample 13-1509 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results
// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/18-0439_filtered.vcf --sample 18-0439 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results
// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/18-2109_filtered.vcf --sample 18-2109 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results
// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/21-0583.final.vcf --sample 21-0583 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results

// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --analysis M

// nextflow run ../annotation_nextflow/main.nf -profile hg19 --final_vcf head_21-4834-CES178.final.vcf --sample 21-4834-CES178 -with-report report.html -with-trace -with-timeline timeline.html -with-dag flowchart.png


// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --analysis A --input /mnt/genetica6/reanotacion/vcfs/ --output /mnt/genetica6/reanotacion/results/ -with-report report.7.html