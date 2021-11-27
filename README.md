<div align="center"><h1>Reproduction d'analyses RNA-seq à l'aide d'un workflow</h1></div>

<div align="justify">
  <p>
    Dans le cadre de notre formation en M2 AMI2B à l'Université Paris-Saclay, nous avons été amenés à réaliser un workflow d'analyses RNA-Seq pour l'UE Hackathon     Reproductible. L'objectif de ce Projet conssite à reproduire les résultats des analyses décrites dans ces deux articles :
    
   * [Recurrent mutations at codon 625 of the splicing factor SF3B1 in uveal melanoma](https://pubmed.ncbi.nlm.nih.gov/23313955), Harbour et al. (2013)
   * [SF3B1 mutations are associated with alternative splicing in uveal melanoma](https://pubmed.ncbi.nlm.nih.gov/23861464), Marais et al. (2013)
  </p>
</div>

<div align="left"><h2>Utilisation du workflow</h2></div>

<div align="justify">
  <p>
    Le workflow s'exécute dans Nextflow et fait appel à Docker pour les conteneurs. Nextflow est lancé depuis Conda (Bioconda). Pour exécuter le workflow, il faut     donc au préalable avoir installé Conda, Nextflow et Docker sur sa machine. La configuration du workflow proposée nécessite d'avoir au minimum 16 CPUs et 50 GB     de mémoire vive.  <br>
    Procédure à suivre pour lancer le workflow : 
    
    ```
    $ conda activate nextflow
    $ nextflow run main.nf -resume
    ```
  </p>
</div>


<div align="left"><h2>Ressources utilisées</h2></div>

<div align="justify">
  <p>
    
  Comme expliqué précédemment, le workflow que nous proposons utilise [Nextflow](https://nextflow.io/) comme Workflow Management System. Nextflow peut être lancé   depuis [Conda](https://conda.io).
    
  Les [Dockers](https://www.docker.com/en) utilisés sont : 
     
   * [pprietob/star-nf](https://hub.docker.com/r/pprietob/star-nf) : utilisé pour récuperer les fichiers fastq, associer les reads 1 et 2, créer l'index de référence et réaliser le mapping des reads
   * [evolbioinfo/samtools](https://hub.docker.com/r/evolbioinfo/samtools) (version 1.14) : utilisé pour indexer le mapping
   * [evolbioinfo/subread:v2.0.1](https://hub.docker.com/r/evolbioinfo/subread) (version 2.0.1) : utilisé pour réaliser la matrice de comptage
   * [evolbioinfo/deseq2:v1.28.1](https://hub.docker.com/r/evolbioinfo/deseq2) (version 1.28.1) : utilisé pour l'analyse statistique R
  
  Les données biologiques utilisées dans notre étude sont celles des SRR contenues dans le fichier [SraAccList_SRA062359.txt](https://github.com/justine66/hackathon/blob/main/SraAccList_SRA062359.txt).
    
  </p>
</div>


<div align="left"><h2>Auteurs</h2></div>

<div align="justify">
  <p>
CARO Hugo, MICHAU Thomas-Sylvestre, NOURRY Justine, PATAT Cedric
  </p>
</div>
