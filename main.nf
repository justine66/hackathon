params.project = "SRA062359" #numeo SRA fournis par l'article
params.resultdir = 'results' #repertoire de sortie des reultats
params.list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT"] #liste de tous les chromosomes humain
projectSRId = params.project
list = params.list


process getSRAIDs {
	
	publishDir params.resultdir, mode: 'copy'

	input:
	val projectID from projectSRId
	
	output:
	file 'sra.txt' into sraIDs
	
	script:
	"""
	esearch -db sra -query $projectID  | efetch --format runinfo | grep SRR | cut -d ',' -f 1 > sra.txt
	"""
}

sraIDs.splitText().map { it -> it.trim() }.filter(  ~/^SRR62858.*/ ).set { singleSRAId}

process fastqDump {
	
	publishDir params.resultdir, mode: 'copy'

	input:
	val id from singleSRAId

	output:
	tuple file('*1.fastq.gz'),file('*2.fastq.gz') into reads #associe les read 1 et 2 dans un tuple
	val id into idmapping #recupere les SRR associe aux reads

	script:
	"""
	parallel-fastq-dump --sra-id $id --threads ${task.cpus} --split-files --gzip;
	"""	
}
readss = idmapping.join(reads) #ajoute les SRR aux tuples des reads pour former des tuble reads1 reads2 SRR

process chromosome {

    input:
    val chr from list

    output:
    file 'Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz' into chrfasta #recuperes tous les chromosomes telecharges dans 1 channel

    script: 
    """
    wget -o ${chr}.fa.gz "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz"
    """
}

process mergechr {
	
    publishDir params.resultdir, mode: 'copy'

    input:
    file allchr from chrfasta.collect() #attend que tous les chromosomes soit telecharges

    output:
    file 'ref.fa' into fasta #unique fichier contenant toutes les chromosomes

    script:
    """
    gunzip -c ${allchr}> ref.fa
    """
}

process gtf {
    
    output:
    file 'annot.gtf' into human_genome

    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
    gunzip -c Homo_sapiens.GRCh38.104.chr.gtf.gz > annot.gtf
    """
}



process index{
	publishDir params.resultdir, mode: 'copy'

	input:
	file c from fasta
	file annot from human_genome

	output:
	file 'ref/' into index #renvoie un unique repertoire contenant tous les fichiers de l'index de reference

	script: 
	"""
	mkdir ref
	chmod +x ref
	STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ${c} --sjdbGTFfile ${annot}
	"""
}

process mapping {
	publishDir params.resultdir, mode: 'copy'

	input:
	tuple file (r1), file (r2), val id  from reads
	file ref from index

	output:
	file '*.bam' into lbam		#recupere les fichier bam pour l'indexation samtools
	file '*.bam' into alignedReads	#recupere les fichier bam pour le comptage

	script :
	"""
	STAR --outSAMstrandField intronMotif \
	--outFilterMismatchNmax 4 \
	--outFilterMultimapNmax 10 \
	--genomeDir ${ref}\
	--readFilesIn <(gunzip -c ${r1}) <(gunzip -c ${r2}) \
	--runThreadN ${task.cpus} \
	--outSAMunmapped None \
	--outSAMtype BAM SortedByCoordinate \
	--outStd BAM_SortedByCoordinate \
	--genomeLoad NoSharedMemory \
	--limitBAMsortRAM 50000000000 \
	> ${id}.bam
	"""
}

process mapping2 {
	publishDir params.resultdir, mode: 'copy'
	
	input:
	file bam from lbam

	output:
	file '*.bai' into map

	script:
	"""
	samtools index ${bam}
	"""
}

process count {

    publishDir params.resultdir, mode: 'copy'

    input:
    file bam from alignedReads.collect() #attend que tous les fichiers bam soit disponible
    file gtf from human_genome

    output:
    tuple file ('output.counts'), file ('output.counts.summary') into countData

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a ${gtf} -o output.counts ${bam}
    """
}

process stat_analysis {
    publishDir params.resultdir, mode: 'copy'

    input:
    tuple file ('output.counts'), file (summary) from countData

    output:
    tuple file ('heatmap.pdf'), file ('d.vst_1000genes.txt') into ana_stat

    script:
    """
    Rscript /ifb/data/mydatalocal/stat_analysis.R\
        --count_file ${'output.counts'} \
        --outdir ./ \
        --cores ${task.cpus}

    """
}
