//=============================================================================
// PARAMS
//=============================================================================
println "Input Directory: ${params.input_scrna}"
println "Export to: ${params.output_scrna}"
println "Manifest file: ${params.manifest}"
nextflow.enable.dsl=2
//=============================================================================
// HELPER FUNCTIONS
//=============================================================================


//=============================================================================
// PROCESSES
//=============================================================================

process scrna_integration {

   
    publishDir params.output_scrna, 
        mode: 'copy',
        pattern: '*.png',                    
        saveAs: { filename -> 
            "plots/${filename}" } 

    input:
    path scrna_data        
    val manifest
    val conda_loc          

    output:
    path "GAND_integrated.png", emit: plot
    path "*.log", emit: logs, optional: true  // Optional log files

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
    
    scRNA.py --path ${scrna_data} --manifest_name ${manifest} 2>&1 | tee analysis.log
    """
}

//=============================================================================
// MAIN WORKFLOW
//=============================================================================

workflow gand_sc {
    // Create channels properly
    data_directory = file("data/scRNA")
    
    // Call the process
    scrna_integration(
        data_directory,
        params.manifest,      
        params.conda_loc
    )
}