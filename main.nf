//=============================================================================
// PARAMS
//=============================================================================
println "Dev mode: ${params.dev_mode}"
//=============================================================================
// WORKFLOW
//=============================================================================
nextflow.enable.dsl=2

include { gand_bw } from './workflows/gand_bw.nf'
include { gand_sc } from './workflows/gand_sc.nf'

workflow {
    // Running ChIP analysis from BigWig Files
    gand_bw()

    // Partially related single Cell Analysis
    // This is here at the moment for simplicity
    // gand_sc()
}