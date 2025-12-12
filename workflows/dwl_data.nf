nextflow.enable.dsl=2
include { isMac} from "${baseDir}/lib/utils.nf"

//=============================================================================
// PROCESSES
//=============================================================================

process download_scrna {
    publishDir "${output_dir}", mode: 'copy', overwrite: true
    
    input:
    val url
    val output_dir
    
    output:
    path "*"
    
    script:
    def useMac = isMac()
    """
    echo "Download single cell RNA-seq data"
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
    publishDir "${output_dir}", mode: 'copy', overwrite: true
    
    input:
    val url
    val output_dir
    
    output:
    path "*"
    
    script:
    def useMac = isMac()
    """
    echo "Downloading single cell RNA-seq reference data"
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
    take:
    scrna_output_dir
    ref_output_dir

    main:
    scrna_channel = Channel.fromList(params.dwl.scrna_url)
        .combine(Channel.value(scrna_output_dir))
    scrna_done = download_scrna(scrna_channel)

    ref_channel = Channel.fromList(params.dwl.ref_url)
        .combine(Channel.value(ref_output_dir))
    ref_done = download_scrna_ref(ref_channel)

    all_done = scrna_done
        .combine(ref_done)
        .map { true }
        .first()
    
    emit:
    status = all_done
}