#!/usr/bin/env nextflow



/*
========================================================================================
                         AF3 Nextflow pipeline
========================================================================================

sh run.sh --fasta1 $PWD/fasta1.fa --fasta2 $PWD/fasta2.fa --outdir $PWD/development_output -resume

----------------------------------------------------------------------------------------
*/

/*
////////////////////////////////////////////////////////////////////
Set default params
////////////////////////////////////////////////////////////////////
*/
params.bin = "${params.base}/bin"
params.date = new Date().format( 'yyyyMMdd' )
if (params.outdir)   { ; } else { exit 1, 'Output directory path not specified!' }
params.results = "${params.outdir}/${params.date}"

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

/*
////////////////////////////////////////////////////////////////////
Import modules
////////////////////////////////////////////////////////////////////
*/
include { MAKE_JSONS } from './modules/af3/main.nf'
include { AF3_ALIGNMENT } from './modules/af3/main.nf'
include { AF3_INFERENCE } from './modules/af3/main.nf'
include { COMBINED_CONFIDENCE_SUMMARY } from './modules/af3/main.nf'

/*
////////////////////////////////////////////////////////////////////
Subworkflows
////////////////////////////////////////////////////////////////////
*/
workflow make_json {

    take : 
        fa1
        fa2

    main : 
        MAKE_JSONS(fa1, fa2)

    emit :
        jsons = MAKE_JSONS.out.jsons.flatten()

}

workflow af3_align {

    take : 
        json
        //msa_dir

    main : 
        AF3_ALIGNMENT(json)

    emit :
        alignments = AF3_ALIGNMENT.out.af3_alignment

}

workflow af3_inference {

    take : 
        data 

    main : 
        AF3_INFERENCE( data )

    emit : 
        inferences = AF3_INFERENCE.out.summary_conf

}

workflow combined_conf_summary { 

    take : 
        data 

    main : 
        COMBINED_CONFIDENCE_SUMMARY( data )

}

/*
////////////////////////////////////////////////////////////////////
Main Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    make_json(params.fasta1, params.fasta2)

    af3_align( make_json.out.jsons )

    af3_inference( af3_align.out.alignments )

    combined_conf_summary( af3_inference.out.inferences.collect() )

}