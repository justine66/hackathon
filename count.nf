process count {

    publishDir params.resultdir, mode: 'copy'

    cpus threads

    input:
    file '.bam' from alignedReads.collect()
    file '.gtf' from human_genome

    output:
    file '.counts' into countData

    script:
    """
    featureCounts -T <CPUS> -t gene -g gene_id -s 0 -a input.gtf -o output.counts input.bam
    """
}