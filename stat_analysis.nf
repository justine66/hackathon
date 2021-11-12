process stat_analysis {

    publishDir params.resultdir, mode: 'copy'

    cpus = $task.cpus

    input:
    file '.counts' into countData counts

    output:
    file "*.pdf" into analysis
    file "*.RData" into analysis
    file "*pca.vals.txt" into analysis
    file "*pca.vals_mqc.tsv" into analysis
    file "*sample.dists.txt" into analysis
    file "*sample.dists_mqc.tsv" into analysis
    file "*.log" into analysis
    file "size_factors" into analysis

    script:
    def label_lower = params.multiqc_label.toLowerCase()
    def label_upper = params.multiqc_label.toUpperCase()
    """
    deseq2_qc.r \\
        --count_file $counts \\
        --outdir ./ \\
        --cores $task.cpus \\
        $options.args
    if [ -f "R_sessionInfo.log" ]; then
        sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" <$pca_header_multiqc >tmp.txt
        sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" tmp.txt
        cat tmp.txt *.pca.vals.txt > ${label_lower}.pca.vals_mqc.tsv
        sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" <$clustering_header_multiqc >tmp.txt
        sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" tmp.txt
        cat tmp.txt *.sample.dists.txt > ${label_lower}.sample.dists_mqc.tsv
    fi
    """
}
