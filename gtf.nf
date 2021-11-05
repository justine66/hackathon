process gtf {
    publishDir params.resultdir, mode: 'copy'
    
    cpus=1

    output:
    file 'annot.gtf' into human_genome

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
    gunzip -c Homo_sapiens.GRCh38.104.chr.gtf.gz > annot.gtf
    """
}
