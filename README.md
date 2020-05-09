![alt text](https://static1.squarespace.com/static/591d9c8cbebafbf01b1e28f9/t/5eae1f2c2bffcf20e605efea/1588469550974/image+%2811%29.png?format=2500w "LABRAT")

# LABRAT <br/> <br/>Lightweight Alignment-Based Reckoning of Alternative Three-prime ends

## Overview

LABRAT is designed to quantify the usage of alternative polyadenylation and cleavage sites in RNAseq data and identify genes whose usage of these sites varies across experimental conditions. <br/>

It takes advantage of the kmer-based, quasi-mapping of reads to transcripts performed by [salmon](https://combine-lab.github.io/salmon/). When it comes to quantifying the usage of alternative 3' UTRs, this strategy has many advantages over classical methods that count the number of counts contained within transcript regions after alignment to the transcriptome. Since alternative 3' UTRs often contain large amounts of sequence in common between them, many reads will map to multiple 3' UTRs, reducing the discriminatory power they contain when align-then-count methods are used.  Transcript abundance quantification with salmon's "lightweight alignments" circumvents this issue to allow accurate quantification of alternative 3' ends.

LABRAT quantifies alternative polyadenylation (APA) site usage by assigning a "psi" value (apologies to the alternative splicing field) for each gene in each sample. Psi values of 0 indicate exclusive usage of the most upstream APA site while psi values of 1 indicate exclusive usage of the most downstream APA site.  When comparing psi values across conditions, an increase in psi reflects increased usage of downstream APA sites while a decrease in psi reflects increased usage of upstream APA sites.  LABRAT uses psi values in experimental replicates to identify genes whose APA site usage significantly differs between experimental conditions.

## Installation

LABRAT is purely python-based (python3), but requires a number of non-standard python modules.  These are most easily installed with [conda](https://docs.conda.io/projects/conda/en/latest/index.html). They are listed below. Versions of each module that are known to be supported are listed, but other versions may work as well.

- python 3.6
- gffutils 0.9
- numpy 1.18.1
- biopython 1.69
- pandas 1.0.3
- statsmodels 0.10.2
- scipy 1.3.1
- salmon (not a python module, but installable using conda)

A safe and easy way to install all of these would be to create a conda environment containing them. We have provided a configuration file that contains all the information needed to setup a LABRAT-ready environment.

```
conda env create -f labratenv.yml
```

This will create a conda environment called 'labrat' that contains all the necessary modules.  To activate the environment, type

```
source activate labrat
```

Then, to make sure you are ready to go, ask for the help options in the LABRAT script by typing

```
python LABRAT.py -h
```

If you see something similar to

```
usage: LABRAT.py [-h] [--mode {makeTFfasta,runSalmon,calculatepsi}]
                 [--gff GFF] [--genomefasta GENOMEFASTA] [--lasttwoexons]
                 [--txfasta TXFASTA] [--reads1 READS1] [--reads2 READS2]
                 [--samplename SAMPLENAME] [--threads THREADS]
                 [--salmondir SALMONDIR] [--sampconds SAMPCONDS]

optional arguments:
  -h, --help            show this help message and exit
  --mode {makeTFfasta,runSalmon,calculatepsi}
  --gff GFF             GFF of transcript annotation. Needed for makeTFfasta
                        and calculatepsi.
```

then you are good to go.  If you get an **ImportError**, one or more of the modules did not install correctly.  In that case, using an alternative package manager like **pip** may help get all these modules installed.

## Annotations

LABRAT requires a genome annotation in gff3 format. We **strongly** recommend using [Gencode](www.gencodegenes.org) annotations. LABRAT is expecting certain tags and attributes to be associated with transcripts in the GFF (e.g. 'protein_coding', etc.) that are not reliably present in GFFs from other sources. Unfortunately, Gencode annotations are only available for human and mouse genomes.

For *Drosophila* annoatations, a separate LABRAT, LABRAT_dm6annotation.py, is provided. This version will accept GFF files from [Ensembl](ftp://ftp.ensembl.org/pub/release-99/gff3/drosophila_melanogaster/). It *should* also work with Ensembl GFF files for other species, but this has not been rigorously tested.

## Running LABRAT

Running LABRAT consists of three steps:</br>
- makeTFfasta
- runSalmon
- calculatepsi

### makeTFfasta

The first step consists of making a fasta file of transcripts that will later be quantified by salmon.  This is done using the following command.

```
python LABRAT.py --mode makeTFfasta --gff <genomegff> --genomefasta <genome sequence in fasta format> --lasttwoexons
```

This will create a database that stores information about the gff using [gffutils](https://daler.github.io/gffutils/index.html). Initial creation of this database can take up to an hour, but it is written to disk so that it does not have to be created in future runs.

The option ```--lasttwoexons``` is optional, but recommended. If included, it tells LABRAT to only consider the last two exons of transcripts for future quantification.  This may be important because it removes assumptions present in the GFF file about relationships between alternative splicing outcomes happening early in a transcript and polyA site choice.

Only protein-coding transcripts and those with confidently defined ends (i.e. those that are not tagged with 'mRNA_end_NF') will end up in this fasta. The output of this step will be a file named either 'TFseqs.fa' or 'wholetranscriptseqs.fa', depending on whether or not the ```--lasttwoexons``` flag was used.

### runSalmon

After creating a fasta file, transcript abundances are calculated using salmon. Reads can be either fastq or fasta, and either gzipped or not. If gzipped, filenames must end in '.gz'. Single and paired end reads are supported. **This step should be performed in a clean, empty directory.**

```
python LABRAT.py --mode runSalmon --txfasta <output of makeTFfasta> --reads1 <comma separated list of forward read files> --reads2 <Optional, comma separated list of reverse read files> --samplename <comma separated list of sample names> --threads <number of threads to use>
```

As an example:

```
python LABRAT.py --mode runSalmon --txfasta TFseqs.fa --reads1 CondARep1_1.fq.gz,CondARep2_1.fq.gz,CondARep3_1.fq.gz,CondBRep1_1.fq.gz,CondBRep2_1.fq.gz,CondBRep3_1.fq.gz --reads2 CondARep1_2.fq.gz,CondARep2_2.fq.gz,CondARep3_2.fq.gz,CondBRep1_2.fq.gz,CondBRep2_2.fq.gz,CondBRep3_2.fq.gz --samplename CondARep1,CondARep2,CondARep3,CondBRep1,CondBRep2,CondBRep3 --threads 8
```

The output of this step will be directories containing salmon quantification files, with one directory per sample.

### calculatepsi

The final step is the calculation of psi values for each gene in each sample, and the identification of genes that show significantly different psi values across conditions.  Transcripts whose 3' ends are less than 25 nt from each other are grouped together during this step and counted as using the same polyadenylation site.

sampconds is a two column, tab-delimited file with sample names split by condition.  ConditionA samples should be in the first column and ConditionB samples in the second column. Reported delta psi values will be ConditionB - ConditionA. These sample names should match those given to runSalmon with ```--samplename```.

salmondir should be a directory that contains **ALL** of the salmon output directories created by runSalmon. **It should contain no other directories besides, optionally, one called 'txfasta.idx'.**

```
python LABRAT.py --mode calculatepsi --salmondir <directory of salmon outputs> --sampconds <sampconds file> --gff <genomegff>
```

## Expected output files

The main output file is named '**LABRAT.psis.pval**'. It contains the psi value for every gene in every condition.  Genes that did not meet an expression cutoff (TPM < 5) in a condition have psi values of NA.  This file also contains delta psi values (ConditionB - ConditionA), a raw pvalue, and a Benjamini-Hochberg corrected pvalue (FDR).

Additionally, this file contains a column called 'genetype'. This contains information about whether the alternative polyadenylation sites in this gene are contained within the same exon (TUTR) or within different exons (ALE). If there are only 2 APA sites for this gene, the gene must be labeled as either TUTR or ALE.  If there are more than 2, if *all* sites are contained within either the same exon, the gene is labeled TUTR.  If *all* sites are contained within different exons, it is labeled ALE. If neither of these is true, the gene is labeled 'mixed'.

A secondary output file is name **'numberofposfactors.txt'**.  This file contains information about the transcripts that were assigned to each APA site. The column 'numberofposfactors' indicates the number of distinct (separated by at least 25 nt) APA sites found for this gene. The column 'txids' has the IDs of the transcripts assigned to each APA site.  APA sites are separated by underscores and transcripts assigned to the same site are separated by commas. In this column, APA sites are ordered from most upstream to most downstream.  The final column contains the distance between the APA sites, but only if the gene contains just 2 sites.

## Advanced Usage

LABRAT identifies genes with significantly different psi values across conditions using a linear mixed effects model. When you are simply comparing across conditions, this is essentially a t-test. However, an advantage of this approach is that it allows for the inclusion of covariates when testing for significance. This has been tested and does seem to work, although currently, incorporating covariates requires editing the code directly (see the model formulae in the getdpsis() function).

Future releases aim to incorporate this in a more user-friendly fashion.

