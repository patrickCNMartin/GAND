//=============================================================================
// WORKFLOW
//=============================================================================
nextflow.enable.dsl=2

include { gand_bw } from './workflows/gand_bw.nf'
//include { gand_scpy } from './workflows/gand_scpy.nf'
include { gand_scR } from './workflows/gand_scR.nf'

workflow {
    // Running ChIP analysis from BigWig Files
    if (params.bigwig.run) {
        println "Run BigWig Analysis"
        gand_bw()

    }
    if (params.scrna_r.run) {
        println "Run Single Cell Analysis"
        gand_scR()
    }
    println "Completed Run"
}