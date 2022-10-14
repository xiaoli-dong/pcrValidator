# pcrValidator
pcrValidator is a stand-alone and automated pipleine developed in Python3. It was designed to assess the DNA oligonucleotide sequences of primers and probes used in diagnostic qPCR or conventional PCR assays in ProvLab, Alberta. It firstly searches user supplied primer and probe sequences against the target DNA sequences, which can be genome or gene sequences, to identify the target regions of the DNA sequences. Secondly, the identified target regions, amplicons, were then run through Mafft program to generate multiple sequence alignments. Based on the amplicon multiple sequence alignments, pcrValidator produces reports on assay performance and sequence variants in the primer and probe sites.

![pcr_validator_workflow](https://user-images.githubusercontent.com/52679027/194174235-b1fe2c9e-2cd1-4ac3-a8d7-73db311fa13f.png)

# Required Third-party software
* [ThermonucleotideBLAST](https://github.com/jgans/thermonucleotideBLAST): a software program for searching a target database of nucleic acid sequences using an assay-specific query.
* [mafft](https://mafft.cbrc.jp/alignment/software/): a multiple sequence alignment program

# Installation

* Firstly, you need to make sure that tntblast and mafft programs are installed and added into your PATH
* Download pcrValidator program and move into the pcrValidator directory
* Run the following commands to create a virtual environtment, activate the virtual environment, and install the dependency packages into the newly created virtual environment:

```
#create a virtual environment and activate it
python3 -m venv .venv

#activate the virtual environment
source .venv/bin/activate

# install required packages contained in requirement.product.txt, then you are ready to run the program
pip install -r requirements.product.txt
```

Now, you are ready to run the program

# Quick start
 pcrValidator allows users to choose blastn or tnbblast to search PCR template sequences or primer and probe sequences against the user defined fasta format reference sequences to identity target amplicons. When the tntblast is chosen, the program only requires user to provide a assay.csv file, in which each line describes a PCR assays and define the location of the fasta format reference file. When the blastn is chosen, the program requires both the assay.csv file and a fasta format PCR template file, which contains one or multiple PCR template sequences. The assay.csv file has multiple lines and each line has eight columns (for taqman assay) or six columns (for conventional PCR assay) to describe the assay. assay.csv file is in the format of:
 
 ```
 assay_name,forward_primer_name,forward_primer_seq,reverse_primer_name,reverse_primer_seq,probe_name,probe_seq 
```

See below for an example assay csv file

```
test1,test1_Fwd,GGAAAATGTAAAGACAACGAATACAG,test1_Rev,GCTATCACATAATCTGGAAGCGTA,test1_Probe,AAGCCGTAATCTATGTTGTCTATCGTGTCC
test2,test2_Fwd,CACACCGTCTCTTCCACAGA,test2_Rev,GATACAGGTTAATTTCCACATCG,test2_Probe,AACCCGTCGTAACCAGCAATACATTT
```

Notes:

* all oligo sequences should be written in the 5' to 3' orientation
* degenerate nucleotides are permitted in the assay oligo sequences
* the probe name and probe sequence can be omitted for a conventional PCR

Example commands to run the program

```
# check the script input options:
python fetch_gbk.py -h
python pcrValidator.py -h

# download the genbank format nucleotide sequences from NCBI
python fetch_gbk.py -q "txid138948[Organism]" -d nucleotide -e test@example.com -p mpx -r gb 

# run the pcrValidator program
python pcrValidator.py -a assay.csv -g mpx.fasta -o results -p mpx -c 0 
```

The pcrValidator program generates five reports for each assay:

* pcr report (*_pcr.tsv)
* assay report (*_assay_report.tsv)
* forward primer variant report (*_fwd_variants.tsv)
* reverse primer variant report (*_rev_variants.tsv)
* proble variant report (*_proble_variants.tsv)

In addition, the program outputs also include the identified fasta format amplicon sequence file and its fasta format multiple sequence alignment file, which can be visualized by some other software such as [JalView](https://www.jalview.org/)

The example outputs can be accessed from [test directory](https://github.com/xiaoli-dong/pcrValidator/tree/main/test/results)