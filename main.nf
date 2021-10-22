params.project = "SRA062359"

params.resultdir = 'results'

params.list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","Mt"]

projectSRId = params.project

list = params.list



process getSRAIDs {
	
	cpus 1

	input:
	val projectID from projectSRId
	
	output:
	file 'sra.txt' into sraIDs
	
	script:
	"""
	esearch -db sra -query $projectID  | efetch --format runinfo | grep SRR | cut -d ',' -f 1 > sra.txt
	"""
}

sraIDs.splitText().map { it -> it.trim() }.set { singleSRAId }

process fastqDump {
	
	publishDir params.resultdir, mode: 'copy'


	input:
	val id from singleSRAId

	output:
	file '*.fastq.gz' into reads

	script:
	"""
	parallel-fastq-dump --sra-id $id --threads ${task.cpus} --gzip --split-files
	"""	
}
process chromosome{

	input:
	val chr from list

	output:
	file 'ref.fa' into fasta

	script: 
	"""
	wget -o $chr.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.!${chr}.fa.gz | gunzip -c *.fa.gz > ref.fa    
	""""
}
process index{
	input:
	file c from fasta

	output:
	file "format a trouver" into index 
	script: 
	"""
	STAR --runThreadN <nb cpus> --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ref.fa
	"""
}




