params.project = "SRA062359"

params.resultdir = 'results'

params.list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT"]

projectSRId = params.project

list = params.list


process getSRAIDs {
	
	cpus 1

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

sraIDs.splitText().map { it -> it.trim() }.filter(  ~/^SRR62858.*/ ).set { singleSRAId }

/*process fastqDump {
	
	publishDir params.resultdir, mode: 'copy'

	input:
	val id from singleSRAId

	output:
	file '*.fastq.gz' into reads

	script:
	"""
	parallel-fastq-dump --sra-id $id --threads ${task.cpus} --split-files --gzip
	"""	
}*/

process chromosome {

    input:
    val chr from list

    output:
    file 'Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz' into chrfasta

    script: 
    """
    wget -o ${chr}.fa.gz "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz"
    """
}

process mergechr {

    input:
    file allchr from chrfasta.collect()

    output:
    file 'ref.fa' into fasta

    script:
    """
    gunzip -c ${allchr}> ref.fa
    """
}

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

process index{


    input:
    file c from fasta
    file annot from human_genome

    output:
    file '*.txt' into index
 
    script: 
    """
    mkdir ref
    STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ${c} --sjdbGTFfile ${annot}
    """
}

/*process star {

    publishDir params.resultdir, mode: 'copy'

    input:
    file read from reads
    file index from index

    output:
    file '.bam' into alignedReads

    script:
    readName = read.toString() - ~/(.fastq.gz)?$/

    """
    STAR --genomeDir $index --readFilesIn $read --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $readName
    """
}*/



