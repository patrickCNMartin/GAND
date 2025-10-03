nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process seurat_integration {
    conda "${params.scrna_env_r}"
    publishDir params.tmp_scrna, 
        mode: 'copy',
        pattern: '*.rds' 

    input:
    path scrna_data        
    val manifest

    output:
    path "GAND_preprocessed.rds", emit: preprocessed
    path "GAND_seurat_integrated.rds", emit: integrated

    // Note that we are going to use conda for now
    // we could also parse renv.lock files
    script:
    """
    seurat_scRNA.R --input_dir ${scrna_data} --manifest ${manifest} 
    """
}


//=============================================================================
// MAIN WORKFLOW
//=============================================================================
include { RMARKDOWNNOTEBOOK } from  '../modules/nf-core/rmarkdownnotebook'
workflow gand_scR {
    // Create channels properly
    data_directory = Channel.fromPath(params.input_scrna)
    println "${params.input_scrna}"
    //println data_directory
    manifest = params.manifest
    // Call the process
    integrated = seurat_integration(
        data_directory,
        manifest)
    // building reports
    // if (params.build_report) {
    //     // Create a channel for the notebook template
    //     notebook_ch = Channel.fromPath(params.report_template)

    //     // Pass optional parameters into the Rmd
    //     params_map_ch = Channel.of([
    //         title: "Seurat scRNA Report",
    //         author: "Molly Easter - Patrick CN Martin",
    //         rds_file: "GAND_seurat_integrated.rds"
    //     ])

    //     // Pass the integrated RDS file to the module
    //     RMARKDOWNNOTEBOOK(
    //         meta: [ id: 'scrna_report' ],
    //         notebook: notebook_ch,
    //         parameters: params_map_ch,
    //         input_files: integrated.integrated,
    //         conda_loc: params.conda_loc
    //     )
    // }
}