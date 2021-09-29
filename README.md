# TTVnucleoTranslations
Is a nucleotide to protein translator for torque teno viruses (TTVs). The pipeline is a Nextflow DSL1 based and translates sequences into 6 frame protein sequences and using a TTV database to help with selection and cutting of ORF (Open Reading Frame) regions.

 ## Software requirements 
 All versions of the softwares should be compatible with the pipeline. Currently I don't have any conflicts between softwares. 
 - [Nextflow DSL1](https://www.nextflow.io/)
 - [Python3](https://www.python.org/downloads/)
 - [seqtk](https://github.com/lh3/seqtk)
 - [BLAST/blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 - [gotranseq](https://github.com/feliixx/gotranseq)


#### Latest software versions used in testing the pipeline
Softwares
| Software     | Version      |
| --------     | -------      |
| nextflow     | 19.07.0.5106 |
| gotranseq    | 0.2          |
| BLAST/blastp | 2.9.0+       |
| python3      | 3.7.8        |
| seqtk        | 1.3-r106     |

## Database
In the folder called "TTV_prot_DB" we have the TTV database containing 289 sequences. If you need to update the database with more sequences you can use the following command bellow.
```
makeblastdb -in <reference.fasta> -dbtype prot -parse_seqids -out TTV_prot_ORFs -title "TTV_prot_ORFs"
```
When running this command, replace <reference.fasta> with the file name of your TTV sequences and when finished, replace the content in "TTV_prot_DB".

## Running the pipeline
The user should create two folders one called 'input_fasta' store all samples with nucleotide sequences in input_fasta. All files must be in fasta format with the file extension ".fasta". 

To run the pipeline in command line:
```
nextflow -C TTVnucleoTranslations.config run TTVnucleoTranslations.nf
```
To run the pipeline in command line and resume from cache memory:
```
nextflow -C TTVnucleoTranslations.config run TTVnucleoTranslations.nf -resume
```

### Output folder
In the output folder "prot_output" you can find your final results with protein sequences in "final_ORFs_results" and "multiORF_fasta". In "final_ORFs_results" sequences are divided by ORF region in separate files and "multiORF_fasta" with all ORFs in one file.
