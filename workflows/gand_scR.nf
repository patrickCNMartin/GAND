//=============================================================================
// PARAMS
//=============================================================================
println "Input Directory: ${params.input_scrna}"
println "Store Intermediate Objects: ${params.tmp_scrna}"
println "Export to: ${params.output_scrna}"
println "Manifest file: ${params.manifest}"
nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process seurat_integration {
    publishDir params.tmp_scrna, 
        mode: 'copy',
        pattern: '*.rds' 

    input:
    path scrna_data        
    val manifest
    val conda_loc          

    output:
    path "GAND_seurat_integrated.rds", emit: integrated

    script:
    """
    if [ "${params.dev_mode}" = "true" ]; then
        echo 'Development setup'
        source ${conda_loc}
        conda activate gand_scrna
    else
        echo 'Production setup'
        export DEBUG=0
    fi
    
    seurat_scRNA.R --input_dir ${scrna_data} --manifest ${manifest} 
    """
}


//=============================================================================
// MAIN WORKFLOW
//=============================================================================
include { RMARKDOWNNOTEBOOK } from  '../modules/nf-core/rmarkdownnotebook'
workflow gand_sc {
    // Create channels properly
    data_directory = file("data/scRNA")
    // Call the process
    integrated = seurat_integration(
        data_directory,
        params.manifest,      
        params.conda_loc
    )
    // building reports
    if (params.build_report) {
        // Create a channel for the notebook template
        notebook_ch = Channel.fromPath(params.report_template)

        // Pass optional parameters into the Rmd
        params_map_ch = Channel.of([
            title: "Seurat scRNA Report",
            author: "Molly Easter - Patrick CN Martin",
            rds_file: "GAND_seurat_integrated.rds"
        ])

        // Pass the integrated RDS file to the module
        RMARKDOWNNOTEBOOK(
            meta: [ id: 'scrna_report' ],
            notebook: notebook_ch,
            parameters: params_map_ch,
            input_files: integrated.integrated,
            conda_loc: params.conda_loc
        )
    }
}