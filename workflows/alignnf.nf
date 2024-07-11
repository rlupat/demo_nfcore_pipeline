/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { BWA_MEM 				 } from '../modules/nf-core/bwa/mem/main' 
include { PICARD_MERGESAMFILES 	 } from '../modules/nf-core/picard/mergesamfiles/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_alignnf_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ALIGNNF {

	take:
	samplesheet

	main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
	
    // Read in samplesheet
    ch_fastq = Channel.fromPath(params.input)
        .splitCsv( header:true )
        .map { row -> 
				meta = [:]
                meta.id = "$row.patient"+"-"+"$row.sample"+"-L"+"$row.lane"
                meta.patient = "$row.patient"
                meta.sample = "$row.sample"
                meta.lane = "$row.lane"
                if (!row.fastq_2) {
                    return [ meta + [ single_end:true ], [ row.fastq_1 ] ]
                } else {
                    return [ meta + [ single_end:false ], [ row.fastq_1, row.fastq_2 ] ]
                }
        }

	//println "ch_fastq: $ch_fastq"
	ch_fastq.view()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq
    )

    ch_bwa_index = Channel.fromPath(params.bwa_index).map {it -> [[id:'bwa'], it]}.collect()
    ch_genome_fasta = Channel.fromPath(params.fasta).map {it -> [[id:it.baseName], it]}.collect()

    // 
    // MODULE: Run bwa mem
    //
    BWA_MEM (
        ch_fastq,
		ch_bwa_index,
		ch_genome_fasta,
		true
    )

	// Set channel to merge bam
    ch_bams_merged = BWA_MEM.out.bam
        .map { meta, bam -> 
            meta = meta - meta.subMap('lane') 
            meta.id = "$meta.patient"+"-"+"$meta.sample"
            [ meta , bam ]
        }
        .groupTuple()


	//println "ch_bams_merged: $ch_bams_merged"
	ch_bams_merged.view()

    // 
    // MODULE: Run Picard MergeSamFiles
    //
    PICARD_MERGESAMFILES (
        ch_bams_merged
    )


    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

	emit:
	multiqc_report = MULTIQC.out.report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
