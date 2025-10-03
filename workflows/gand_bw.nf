nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process merge_by_tag {
    conda "${params.bigwig.env}"
    publishDir "${params.bigwig.tmp}", mode: 'copy', overwrite: true
    input:
    tuple val(tag), path(bw_files)
    output:
    tuple val(tag), path("${tag}_merged_hg38_uniqnorm_signal.bg"), emit: merged_bg
    script:
    """
    echo "Merging BigWig Files into BedGraph"
    bigWigMerge ${bw_files.join(' ')} ${tag}_merged_hg38_uniqnorm_signal.bg
    """
}

process peak_from_bdg {
    conda "${params.bigwig.env}"
    publishDir "${params.bigwig.tmp}", mode: 'copy', overwrite: true
    input:
    path bdg_file
    val min_l
    val cutoff
    output:
    path "${bdg_file.baseName}.narrowPeak"
    script:
    """
    echo "Calling peaks from BedGraph"
    macs3 bdgpeakcall -i ${bdg_file} -o ${bdg_file.baseName}.narrowPeak -l ${min_l} -c ${cutoff}
    """
}

process diff_peak {
    conda "${params.bigwig.env}"
    publishDir "${params.bigwig.tmp}", mode: 'copy', overwrite: true
    input:
    tuple path(cond), path(cond_input), path(ctr), path(ctr_input)
    val min_l
    val cutoff
    output:
    path "GAND_GATAD2B_diff_peaks*.*", emit: diff_files
    script:
    """
    echo "Running Diff peaks on merged files"
    macs3 bdgdiff --t1 ${cond} --c1 ${cond_input} --t2 ${ctr} --c2 ${ctr_input} --o-prefix GAND_GATAD2B_diff_peaks --cutoff ${cutoff} --min-len ${min_l}
    """
}


//=============================================================================
// WORKFLOW
//=============================================================================

workflow gand_bw {

    // Input tags
    tags = Channel.fromList(params.bigwig.tags)

    // All bigwig files
    bw_files = Channel
        .fromPath("${params.bigwig.input}", checkIfExists: true)
        .ifEmpty { error "No BigWig files found at ${params.bigwig.input}" }

    // Create [tag, files] grouped by tag
    tagged_files = tags
        .combine(bw_files)
        .filter { tag, file -> file.toString().toLowerCase().contains(tag.toLowerCase()) }
        .groupTuple()
        .map { tag, files -> tuple(tag, files) }
        .ifEmpty { error "No files matched any of the tags in ${params.bigwig.tags}" }

    // Merge bigwigs → emit (tag, merged_bg)
    merged = merge_by_tag(tagged_files)
    merged.merged_bg.view { tag, file -> "MERGED: ${tag} -> ${file}" }
    // Run peak calling
    peak_from_bdg(
        merged.merged_bg.map { it[1] },   // just the bg file, not the tag
        params.bigwig.min_peak_length,
        params.bigwig.peak_cutoff
    )

    // Filter merged outputs by tag name directly (no fromPath)
    cond_peak_ch       = merged.merged_bg.filter { tag, file -> tag == "GAND_GATAD2B" }.map { it[1] }
    cond_input_peak_ch = merged.merged_bg.filter { tag, file -> tag == "GAND_Input"   }.map { it[1] }
    ctr_peak_ch        = merged.merged_bg.filter { tag, file -> tag == "CTR_GATAD2B" }.map { it[1] }
    ctr_input_peak_ch  = merged.merged_bg.filter { tag, file -> tag == "CTR_Input"   }.map { it[1] }

    // Combine into one tuple for diff_peak
    combined_peak_files = cond_peak_ch
        .combine(cond_input_peak_ch)
        .combine(ctr_peak_ch)
        .combine(ctr_input_peak_ch)
        .map { cond, cond_input, ctr, ctr_input ->
            tuple(cond, cond_input, ctr, ctr_input)
        }

    diff_peak(combined_peak_files,
        params.bigwig.min_peak_length,
        params.bigwig.peak_cutoff)
}