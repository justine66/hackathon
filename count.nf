process count {

    publishDir params.resultdir, mode: 'copy'

    cpus = $task.cpus

    input:
    file bam from alignedReads.collect()
    file gtf from human_genome

    output:
    file counts into countData

    script:
    """
    featureCounts -T $cpus -t gene -g gene_id -s 0 -a $gtf -o output.counts $bam
    """
}
