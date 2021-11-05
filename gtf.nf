process gtf {
    publishDir params.resultdir, mode: 'copy'
    
    cpus threads

    output:
    file '*' into human_genome

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
    """
}
