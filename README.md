# scEDIT

Single cell Edit Detection and Identification Tool (scEDIT) for analysis of single cell amplicon DNA-seq data generated on Tapestri platform.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Implementation](#implementation)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

## Introduction

CRISPR Cas9 system has been revolutionizing the gene therapy research and paving the way for personalized medicine. Although CRISPR can silence/correct/modify the gene with greater precision compared to other nuclease off-target and impreseaice activity remains a major safety concern. Bulk DNA sequencing provides a population level view of the CRISPR however details of rare and large deletions are diluted. Single cell DNA-seq Tapestri platform provides a tool to study the CRISPR activity at single cell resolution. Currently, analysis of Tapestri generated can only be done using the custom Tapestri bioinformatics pipeline with restricted availability. scEDIT is lightweight and portable software tool for analysis of CRISPR editing at single cell resolution. 
scEDIT starts with raw sequence data from Tapestry library, extracts and identifies cell barcode and amplicon count in each cell. Furthermore, scEDIT detects and identifies CRISPR induced small indels as well as large deletions.    

## Features

 Main features of scEDIT:

- Processing of raw sequence data (Only other alternative is Tapestri Bioinformatics Pipeline with restricted access).   
- Extraction of cell barcode and generation of amplicon count matrix
- Detection and identification indels and large deletions
- Identifies zygosities of CRISPR-induced mutations
- Lightweight and portable can be run efficiently on inexpensive desktop with multicore/threaded CPU (AMD/Intel) 
- Easy to install and use with conda environment
- Also include tool for separating reads by cell barcode and Amplicon ID for deeper downstream analysis

## Prerequisites

List any prerequisites, for example:

- Python 3.9
- Dependencies in src/scEDIT_conda_env_packageList.yml

## Installation
Download the folder with all the 
Create conda environment using the file src/scEDIT_conda_env_packageList.yml

## Implementation (run in scEDIT conda environment)
 
Usage: python scEDIT.py --mode BaseEdit/InDel --samples path/SampleList.xlsx --outDir path/output_directory

    Single cell Edit detection and Identification Tool (scEDIT) runs two modes of analysis: InDel and BaseEdit
        -InDel mode is used for single cell analysis of CRISPR-Cas induced indels 
        -BaseEdit mode is used for single cell analysis of CRISPR-Cas induced base edits
    
    Each mode required different set of excel files containing amplicon primer information and gRNA details.
    
    The excel files need to have mode specific format.

Above command generates bash scripts in the /full-directory-path/output-dir for each sample in the SampleSheet (Example SampleSheet and Amplicon sheets are provided in scCRISPR_suppFiles)

Run each bash_(*SAMPLEID).sh file using following command
$sh bash_(*SAMPLEID).sh


Final fastq files are stored in /SubjectID/SampleName/result_out/fastqs/
Reads did not qualify the criteria are stored in dump fastq files stored in /SubjectID/SampleName/result_out/fastqs/

Final counts and edit details files stored in /SubjectID/SampleName/result_out/final_count/
final_count folder also contains plots showing cell rank and barcode count. 

## Usage
Single cell CRISPR data analyzed using scEDIT is from work of Elisa ten Hacken et al published in Genome Biology
# High throughput single-cell detection of multiplex CRISPR-edited gene modifications
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02174-1
Single cell CRISPR base editing data analyzed using scEDIT is from work of Jorge D Martin-Rufino et al published in Cell
2023 May 25;186(11):2456-2474.e24. doi: 10.1016/j.cell.2023.03.035. Epub 2023 May 2. 
# Massively parallel base editing to map variant eXects in human hematopoiesis. Cell 186 (2023)
https://www.cell.com/cell/fulltext/S0092-8674(23)00332-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS009286742300332X%3Fshowall%3Dtrue

For converting sra file to fastq file use following command
parallel-fastq-dump --sra-id $srr --threads 12 --split-files --origfmt --gzip --defline-seq '@$ac.$sn  $sn length=$rl' --include-technical --outdir /your_fastq_dir

## License
scEDIT is avaible under MIT licence.

## Contact
Gajendra Suryawanshi
gw.suryawanshi@gmail.com

### Important note
Author has solely created all the work provided here and preliminary results presented here are authors own interpretation of the data. Author's current or previous employers are in any way not associated with work or results presented here. 
