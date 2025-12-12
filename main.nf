//=============================================================================
// WORKFLOW
//=============================================================================
nextflow.enable.dsl=2

include { dwl_data } from './workflows/dwl_data.nf'
include { gand_integration } from './workflows/gand_integration.nf'
include { gand_modelling} from './workflows/gand_modelling.nf'
include { build_report } from './workflows/build_report.nf'

workflow {
    
    
    if (params.run_download) {
        println "Check Data to be downloaded..."
        dwl = dwl_data(
            params.dwl.output_scrna,
            params.dwl.output_ref)
        dwl_status = dwl.status
    } else {
        dwl_status = Channel.of(true)
    }
    
    if (params.run_integration) {
        println "Integrating and annotating scRNA..."
        annotated = gand_integration(dwl_status)
        annotated_status = annotated.status
    } else {
        annotated_status = Channel.of(true)
    }

    if (params.run_modelling) {
        println "Modelling Cells..."
        models = gand_modelling(annotated_status)
        model_status = models.status
    } else {
        model_status = Channel.of(true)
    }

    if (params.run_report) {
        println "Building Report"
        def data_map = [
            "annotated": params.report.annotated ? file(params.report.annotated) : null,
            "mut_genes": params.report.mut_genes ? file(params.report.mut_genes) : null,
            "gene_sets": params.report.gene_sets ? file(params.report.gene_sets) : null
        ]
        
        build_report(data_map, model_status)
    }
    println "Completed Run"
}