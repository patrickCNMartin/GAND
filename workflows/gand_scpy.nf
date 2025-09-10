nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================

process scrna_integration {
    conda "${params.scrna_env_py}"
    publishDir params.output_scrna, 
        mode: 'copy',
        pattern: '*.png',                    
        saveAs: { filename -> 
            "plots/${filename}" } 
    
    input:
    path scrna_data        
    val manifest       

    output:
    path "GAND_integrated.png", emit: plot
    path "*.log", emit: logs, optional: true  // Optional log files

    script:
    """
    scRNA.py --path ${scrna_data} --manifest_name ${manifest} 2>&1 | tee analysis.log
    """
}

//=============================================================================
// MAIN WORKFLOW
//=============================================================================

workflow gand_scpy {
    // Create channels properly
    data_directory = file("data/scRNA")
    
    // Call the process
    scrna_integration(
        data_directory,
        params.manifest
    )
}