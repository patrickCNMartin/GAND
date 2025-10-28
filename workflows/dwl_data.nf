nextflow.enable.dsl=2
include { isMac} from "${baseDir}/lib/utils.nf"

//=============================================================================
// PROCESSES
//=============================================================================

process download_scrna {
    publishDir "${params.scrna_r.input}", mode: 'copy', overwrite: true
    
    input:
    val url
    
    output:
    path "*"
    
    script:
    def useMac = isMac()
    """
    echo "Download VisiumHD"
    echo "Downloading from: ${url}"
    
    if [ "${useMac}" = "true" ]; then
        curl -L -J -O "${url}"
    else
        wget --content-disposition "${url}"
    fi
    
    downloaded_file=\$(ls -t *.tar | head -n 1)
    echo "Downloaded file: \$downloaded_file"
    
    echo "Extracting: \$downloaded_file"
    tar -xf "\$downloaded_file"
    
    echo "Extraction complete"
    """
}

process download_scrna_ref {
    publishDir "${params.scrna_r.ref}", mode: 'copy', overwrite: true
    
    input:
    val url
    
    output:
    path "*"
    
    script:
    def useMac = isMac()
    """
    echo "Downloading scRNA reference data"
    echo "Downloading from: ${url}"
    
    if [ "${useMac}" = "true" ]; then
        curl -L -O "${url}"
    else
        wget "${url}"
    fi
    
    downloaded_file=\$(ls -t | head -n 1)
    
    echo "Downloaded file: \$downloaded_file"
    """

}

//=============================================================================
// WORKFLOW
//=============================================================================

workflow dwl_data {
    if (params.dwl.dwl_scrna == true) {
        scrna_channel = Channel.fromList(params.dwl.scrna_url)
        download_scrna(scrna_channel)
    }
    if (params.dwl.dwl_ref == true) {
        ref_channel = Channel.fromList(params.dwl.ref_url)
        download_scrna_ref(ref_channel)
    }

}