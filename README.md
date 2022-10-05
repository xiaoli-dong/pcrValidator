# pcrValidator
pcrValidator is a stand-alone and semi-automated pipleine. It was designed to assess the primer and probe oligo nucleotides of diagnostic qPCR or PCR assays in ProvLab. It searches user supplied primer and probes sequences against the target nucleotide sequences file to identify the target amplicons. The identified amplicons were then run through Mafft to generate multiple sequence alignment. Based on the multiple sequene alignments, pcrValidator generates a list of reports:

![pcr_validator_workflow](https://user-images.githubusercontent.com/52679027/194162582-c5b064b6-eaab-4f10-9e50-8194bce3534f.png)

# Required Third-party software
* [tntblast](https://github.com/jgans/thermonucleotideBLAST): a software program for searching a target database of nucleic acid sequences using an assay-specific query.
* [mafft](https://mafft.cbrc.jp/alignment/software/): a multiple sequence alignment program

# Installation
Firstly download the program and then create a virtual environtment and then run it
```
#create a virtual environment and activate it
python3 -m venv .vent
source .venv/bin/activate
# install required packages contained in requirement.product.txt, then you are ready to run the program
pip install -r requirements.product.txt
```

# Quick start
 pcrValidator takes an csv file, in which each line describles a PCR assays and a fasta format sequences file, which contains the target sequences. The assay csv file can contain multiple lines and each line describe an assay and  uses the following format: assay_name, forward_primer_name, forward_primer_seq, reverse_primer_name, reverse_primer_seq, probe_name, probe_seq. 

See below an example assay csv file

```
test1,test1_Fwd,GGAAAATGTAAAGACAACGAATACAG,test1_Rev,GCTATCACATAATCTGGAAGCGTA,test1_Probe,AAGCCGTAATCTATGTTGTCTATCGTGTCC
test2,test2_Fwd,CACACCGTCTCTTCCACAGA,test2_Rev,GATACAGGTTAATTTCCACATCG,test2_Probe,AACCCGTCGTAACCAGCAATACATTT
```

Notes

* all oligo sequences should be writen in the 5' to 3' orientation
* degenerate nucleotides are permitted in the assay oligo sequences
* the probe name and probe sequence can be omitted for a conventional PCR

See belwo on how to run the program
```
# check the script input options:
python fetch_gbk.py -h
python pcrValidator.py -h

# download the genbank format nucleotide sequences belong to the organism, whose taxonomic id is 138943 from ncbi
python fetch_gbk.py -q "txid138948[Organism]" -p mpx -r gb -d nucleotide -e test@example.com

# validate the assay
python pcrValidator.py -a assay.csv -g mpx.fasta -o results -p mpx -c 0 
```

The pcrValidator program generates following reports
* pcr report (*_pcr.tsv)
* assay report (*_assay_report.tsv)
* forward primer variant report (*_fwd_variants.tsv)
* reverse primer variant report (*_rev_variants.tsv)
* proble variant report (*_proble_variants.tsv)
