![alt text](https://images.squarespace-cdn.com/content/v1/591d9c8cbebafbf01b1e28f9/1588470361245-PJCCCJTEQ8JXCSIM173R/ke17ZwdGBToddI8pDm48kGcCOoVw7tMeq96q09YNvTBZw-zPPgdn4jUwVcJE1ZvWQUxwkmyExglNqGp0IvTJZamWLI2zvYWH8K3-s_4yszcp2ryTI0HqTOaaUohrI8PIveFD_g_aea0KYeAi1GiEsiYHubgbp9p27L_1ORGVwOoKMshLAGzx4R3EDFOm1kBS/image%2B%252811%2529.jpg?format=2500w "LABRAT")

# LABRAT <br/> <br/>Lightweight Alignment-Based Reckoning of Alternative Three-prime ends

## Overview

LABRAT is designed to quantify the usage of alternative polyadenylation and cleavage sites in RNAseq data and identify genes whose usage of these sites varies across experimental conditions. <br/>

It takes advantage of the kmer-based, quasi-mapping of reads to transcripts performed by [salmon](https://combine-lab.github.io/salmon/). When it comes to quantifying the usage of alternative 3' UTRs, this strategy has many advantages over classical methods that count the number of counts contained within transcript regions after alignment to the transcriptome. Since alternative 3' UTRs often contain large amounts of sequence in common between them, many reads will map to multiple 3' UTRs, reducing the discriminatory power they contain when align-then-count methods are used.  Transcript abundance quantification with salmon's "lightweight alignments" circumvents this issue to allow accurate quantification of alternative 3' ends.

LABRAT quantifies alternative polyadenylation (APA) site usage by assigning a "psi" (𝜓) value (apologies to the alternative splicing field) for each gene in each sample. Psi values of 0 indicate exclusive usage of the most upstream APA site while 𝜓 values of 1 indicate exclusive usage of the most downstream APA site.  When comparing 𝜓 values across conditions, an increase in 𝜓 reflects increased usage of downstream APA sites while a decrease in 𝜓 reflects increased usage of upstream APA sites.  LABRAT uses 𝜓 values in experimental replicates to identify genes whose APA site usage significantly differs between experimental conditions.

## Installation

### Option 1: conda

The *easiest* and *best* way to install LABRAT is using the [conda](https://docs.conda.io/projects/conda/en/latest/index.html) package manager. LABRAT is available on the [bioconda](https://bioconda.github.io/) channel.

```
conda install -c bioconda labrat
```

### Option 2: manual installation

Alternatively, you can download LABRAT directly from this repository. LABRAT is purely python-based (python3), but requires a number of non-standard python modules.  These are most easily installed with [conda](https://docs.conda.io/projects/conda/en/latest/index.html). They are listed below. Versions of each module that are known to be supported are listed, but other versions may work as well.

- python 3.6
- gffutils 0.9
- numpy 1.18.1
- biopython 1.69
- pandas 1.0.3
- statsmodels 0.10.2
- scipy 1.3.1
- salmon 0.14 (not a python module, but installable using conda)

A safe and easy way to install all of these would be to create a conda environment containing them. We have provided a configuration file that contains all the information needed to setup a LABRAT-ready environment.

*Note:* The most recent version of [salmon](https://github.com/COMBINE-lab/salmon/releases) has changed the way that transcriptome indexes are made. LABRAT does not currently support these versions (salmon >= 1.0.0), but will in the future.

```
conda env create -f labratenv.yml
```

This will create a conda environment called 'labrat' that contains all the necessary modules.  To activate the environment, type

```
source activate labrat
```

Uncompress the repository and move into the uncompressed directory. Install LABRAT by typing 

```
python setup.py install
```

Then, to make sure you are ready to go, ask for the help options in the LABRAT script by typing

```
LABRAT.py -h
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

When given to LABRAT, these gff annotation files must **not** be compressed.

## Library types

LABRAT accepts high-throughput RNA sequencing data of two types: RNAseq and 3 prime end sequencing. 3 prime end sequencing is very informative for the study of alternative polyadenylation since the majority of reads are relatively easily assignable to the usage of one particular polyadenylation site. This may lead to more accurate quantification of polyadenylation site usage. However, the majority of publicly available data is of the RNAseq variety.

LABRAT uses a parameter `--librarytype` parameter to distinguish between these two possibilities. The possible values for this parameter are `RNAseq` and `3pseq`. If the value is `RNAseq`, the last two exons of every valid transcript are quantified, and the quantifications used for further calculation of 𝜓 are the length-normalized TPM values. Gene-level expression thresholds are set at 5 TPM.

If the value is `3pseq`, then only the last 300 nt of each transcript are quantified (this is due to the average insert length of many 3' end libraries being 150-300 nt), and the quantifications used for further calculation of 𝜓 are read count values (since length normalization of 3' end data is not necessary or wanted). Gene-level expression thresholds are set at 100 counts.

## Running LABRAT

Running LABRAT consists of three steps:</br>
- makeTFfasta
- runSalmon
- calculatepsi

### makeTFfasta

The first step consists of making a fasta file of transcripts that will later be quantified by salmon.  This is done using the following command.

```
LABRAT.py --mode makeTFfasta --gff <genomegff> --genomefasta <genome sequence in fasta format> --lasttwoexons --librarytype <librarytype>
```

This will create a database that stores information about the gff using [gffutils](https://daler.github.io/gffutils/index.html). Initial creation of this database can take up to several hours, but it is written to disk so that it does not have to be created in future runs. Compressed gff files are not currently supported.

*Important*: If you kill the database creation process before it is finished, you will still have a .db file written. The next time you run LABRAT, it will see that file and think it is complete. This will obviously lead to problems. If this happens, simply delete the .db file. This will force LABRAT to create the .db file again.

If you would like to skip this, gff annotations and pre-built db files for hg38, mm10, and dm6 are available [here](https://www.dropbox.com/sh/qy3jzd00k00w3ga/AACkNE2q3d68sr3wKQDLboORa?dl=0). **Importantly**, LABRAT will expect a gff file and its corresponding db to be located in the same directory.

The option ```--lasttwoexons``` is optional, but recommended. If included, it tells LABRAT to only consider the last two exons of transcripts for future quantification.  This may be important because it removes assumptions present in the GFF file about relationships between alternative splicing outcomes happening early in a transcript and polyA site choice.

Only protein-coding transcripts and those with confidently defined ends (i.e. those that are not tagged with 'mRNA_end_NF') will end up in this fasta. The output of this step will be a file named either 'TFseqs.fa' or 'wholetranscriptseqs.fa', depending on whether or not the ```--lasttwoexons``` flag was used.

### runSalmon

After creating a fasta file, transcript abundances are calculated using salmon. Reads can be either fastq or fasta, and either gzipped or not. If gzipped, filenames must end in '.gz'. Single and paired end reads are supported. **This step should be performed in a clean, empty directory.**

```
LABRAT.py --mode runSalmon --librarytype <librarytype> --txfasta <output of makeTFfasta> --reads1 <comma separated list of forward read files> --reads2 <Optional, comma separated list of reverse read files> --samplename <comma separated list of sample names> --threads <number of threads to use>
```

As an example:

```
LABRAT.py --mode runSalmon --txfasta TFseqs.fa --reads1 CondARep1_1.fq.gz,CondARep2_1.fq.gz,CondARep3_1.fq.gz,CondBRep1_1.fq.gz,CondBRep2_1.fq.gz,CondBRep3_1.fq.gz --reads2 CondARep1_2.fq.gz,CondARep2_2.fq.gz,CondARep3_2.fq.gz,CondBRep1_2.fq.gz,CondBRep2_2.fq.gz,CondBRep3_2.fq.gz --samplename CondARep1,CondARep2,CondARep3,CondBRep1,CondBRep2,CondBRep3 --threads 8
```

The output of this step will be directories containing salmon quantification files, with one directory per sample.

### calculatepsi

The final step is the calculation of 𝜓 values for each gene in each sample, and the identification of genes that show significantly different 𝜓 values across conditions.  Transcripts whose 3' ends are less than 25 nt from each other are grouped together during this step and counted as using the same polyadenylation site.

salmondir should be a directory that contains **ALL** of the salmon output directories created by runSalmon. **It should contain no other directories besides, optionally, one called 'txfasta.idx'.**

sampconds is a tab-delimited file with a column header row that gives information about the samples. Column names are strict, and the first two columns are required. The first column should be named '**sample**' and contain the names of every sample to be compared. These sample names should match those given to runSalmon with ```--samplename```. The second column should be named '**condition**' and should contain two factors identifying the grouping of the samples for the comparison of interest.  These factors must also be given as the arguments **--conditionA** and **--conditionB**.  Delta 𝜓 values will be reported as B-A.  Additional columns representing covariates can be included, but not are not requried. Covariate column names must contain 'covariate' within them. A sample sampconds file is provided below.

| sample | condition | covariate1 |
---------|-----------|--------
| Brain_M1 | Brain| M |
| Brain_F1 | Brain | F |
| Liver_M1 | Liver | M |
| Liver_F1 | Liver | F |
| Liver_F2 | Liver | F |

LABRAT compares 𝜓 values of experimental replicates across experimental conditions to identify genes with statistically significantly different 𝜓 values between conditions.  This is done using a mixed linear effects model that tests the relationship between 𝜓 values and experimental condition. A null model is also created in which the term denoting the experimental condition has been removed.  A likelihood ratio test compares the goodness of fit of these two models to the observed data and assigns a p value for the probability that the real model is a better fit than the null model. In simple comparisons between two conditions, this approach mimics a t-test.  However, this technique has the advantage of being able to easily incorporate covariates into significance testing.  After performing this test on all genes, the raw p values are corrected for multiple hypothesis testing using a Benjamini-Hochsberg correction.

```
LABRAT.py --mode calculatepsi --salmondir <directory of salmon outputs> --sampconds <sampconds file> --conditionA <conditionA> --conditionB <conditionB> --gff <genomegff> --librarytype <librarytype>
```

## Expected output files

The main output file is named '**LABRAT.psis.pval**'. It contains the 𝜓 value for every gene in every condition.  Genes that did not meet an expression cutoff (TPM < 5) in a condition have 𝜓 values of NA.  This file also contains delta 𝜓 values (ConditionB - ConditionA), a raw pvalue, and a Benjamini-Hochberg corrected pvalue (FDR).

Additionally, this file contains a column called 'genetype'. This contains information about whether the alternative polyadenylation sites in this gene are contained within the same exon (TUTR) or within different exons (ALE). If there are only 2 APA sites for this gene, the gene must be labeled as either TUTR or ALE.  If there are more than 2, if *all* sites are contained within either the same exon, the gene is labeled TUTR.  If *all* sites are contained within different exons, it is labeled ALE. If neither of these is true, the gene is labeled 'mixed'.

A secondary output file is name **'numberofposfactors.txt'**.  This file contains information about the transcripts that were assigned to each APA site. The column 'numberofposfactors' indicates the number of distinct (separated by at least 25 nt) APA sites found for this gene. The column 'txids' has the IDs of the transcripts assigned to each APA site.  APA sites are separated by underscores and transcripts assigned to the same site are separated by commas. In this column, APA sites are ordered from most upstream to most downstream.  The final column contains the distance between the APA sites, but only if the gene contains just 2 sites.

## Runtime

If LABRAT is encountering a gff genome annotation file for the first time, it indexes this file using [gffutils](https://github.com/daler/gffutils/). This process can take a few hours, depending on the size of the annotation. However, it only needs to be completed once.  All future runs will automatically make use of a database file written after the indexing completes. Importantly, if indexing is interrupted, this file will still be written, and LABRAT will attempt to use this truncated file in the next run. This will cause problems. To prevent this, if indexing is interrupted, be sure to delete the resulting database file. It can be found at the location of the gff annotation, and ends with '.db' 

To test the runtime requirement of LABRAT, we focused on a dataset that considered two conditions with two replicates per condition. Each sample contained approximately 25 million paired end reads. Using a modern Intel Mac laptop running OSX 10.15 with 12 cores, LABRAT analysis of this data took approximately 25 minutes.  This does not include the time taken to index the genome annotation as described above.

## Using LABRAT with single cell RNAseq data

LABRAT uses [salmon](https://github.com/COMBINE-lab/salmon/releases) to quantify APA from RNAseq data.  LABRATsc uses an analogous tool, [alevin](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y), for quantification of APA from single cell RNAseq data. An example of the use of LABRATsc can be found in section 4 of [this paper](https://www.sciencedirect.com/science/article/pii/S0076687921001312).


LABRATsc compares 𝜓 values between **predefined** groups or clusters of cells. These clusters can be defined using standard approaches such as tSNE or UMAP.  Importantly, quantification with alevin and designation of clusters is not performed by LABRATsc and must be done beforehand.  Following quantification with alevin, transcripts are filtered to keep only those with confidently defined 3' ends exactly as they are in LABRAT-based quantification. Analysis then proceeds upon one of two paths defined by the `--mode` parameter.

If the `--mode` parameter is set to `cellbycell`, then a 𝜓 value is calculated for every (applicable) gene in every cell. In practice, the coverage for most genes in most cells is low or nonexistant. Genes that do not pass a read coverage threshold (indicated by the `--readcount` parameter) have 𝜓 values of NA in that cell. For each gene, 𝜓 values are then compared across cell clusters using the 𝜓 value of individual cells as independent observations.

Alternatively, if the `--mode` parameter is set to `subsampleClusters`, then read counts for each transcript are first summed across all the cells within a cluster. This has the advantage of raising the number of reads associated with each gene, but single cell resolution is lost. Tests to identify genes with regulated APA across cell clusters are performed by creating a distribution of 𝜓 alues for each gene in each cluster through bootstrapping resampling.

### Important considerations for LABRATsc

Many single cell RNAseq libraries are well suited to the quantification of APA because the reads they produce are at or near polyA sites. Still, several important caveats must be considered. 

First, low read depth and dropouts limit the reliable detection and quantification of APA in single cell data. These limiations make the selection of appropriate minimum read thresholds (using the `--readcountfilter` argument) critical. We suggest thresholds of at least 100 counts per gene for cluster-level 𝜓 quantification (`--mode subsampleClusters`) and at least 5 counts per cell for cell-level 𝜓 quantification (`--mode cellbycell`) as starting points, but these thresholds (particularly in the `cellbycell` case) vary considerably between experiments due to a range of technical and biological factors.

Second, while scRNAseq libraries generally capture the 3' ends of mRNAs, they also contain reads that arise due to internal priming on genomically encoded A-rich regions. If these internal priming events occur in close proximity to bonafide polyA sites, they may skew raw 𝜓 values substantially. However, while this may impact the accuracy of raw 𝜓 values for some genes, the relative 𝜓 values between cells should be less affected as the rate of internal priming should be largely consistent across cells.

### Generating input matrices for LABRATsc using alevin

Prior to running LABRATsc, cell-by-isoform count matrices must be produced. This relies on the single-cell transcriptome quantification tool [alevin](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1670-y). For general information about alevin, see its [documentation](https://salmon.readthedocs.io/en/latest/alevin.html). 

The following arguments must be passed to **alevin** for use with LABRATsc, ideally in this order: </br>

- `-l`: Library type. For most single-cell libraries, this will be "ISR".
- `-1`: A list of files containing the forward sequencing reads.
- `-2`: A list of files containing the reverse sequencing reads in the same order as `-1`.
- `--dropseq --chromium --chromiumV3`: One of these depending on the sequencing platform used.
- `-i`: A salmon index, generated with LABRAT using the `--librarytype 3pseq` argument (see above).
- `-p`: Number of threads to be used by alevin.
- `-o`: Output path for each count matrix and metadata
- `--tgMap`: A transcript-to-gene map file, which consists of each transcript ID in the salmon index listed twice per line, separated by tab
- `--fldMean 250`: Expected mean fragment length (250 for consistency with LABRAT's execution of salmon).
- `--fldSD 20`: Expected standard deviation of fragment length (20 for consistency with LABRAT's execution of salmon).
- `--validateMappings`: Enables selective alignment of reads
- `--whitelist`: A whitelist of cell barcodes from a previous analysis to restrict quantification to previous identified valid barcodes (Optional)

### Calculating 𝜓 values with LABRATsc

Transcript counts from one or more single-cell libraries are used to calculate 𝜓 values for every gene with at least two APA sites. As with LABRAT, gene-level 𝜓 values are then compared across conditions to identify genes with significant 𝜓 value changes. As described above, LABRATsc provides two diferent approaches for 𝜓 calculation and testing: per-cell or using subsampling within clusters. The relevant options for quantification of 𝜓 values with LABRATsc are as follows: </br>

- `--mode`: `cellbycell` or `subsampleClusters`, as described above
- `--gff`: path to the gff annotation to be used. It should be the same annotation used to generate the salmon index provided to alevin.
- `--alevindir`: A directory containing alevin quantification subdirectories with one for each sample. The names of these subdirectories will be appended to the cell names in each sample matrix to form a "sample_barcode" cell ID for each cell. An example `alevindir` can be found [here](https://github.com/TaliaferroLab/LABRAT/tree/singlecell/testdata/alevin_example/alevin_out).
- `--conditions`: A tab delimited text file with column names "sample" and "condition". The first column contains cell IDs and the second column contains cell condition or cluster. The cell IDs in the sample column must follow the "sample_barcode" structure described above. Note that unlike LABRAT, LABRATsc does not currently support covariates.  An example `conditions` file can be found [here](https://github.com/TaliaferroLab/LABRAT/blob/singlecell/testdata/alevin_example/conditions.tsv).
- `--readcountfilter`: Minimum read count necessary for calculation of 𝜓 values. Genes that do not pass this threshold will have 𝜓 values of NA. If in `cellbycell` mode, this is the number of reads mapping to a gene in that single cell. If in `subsampleClusters` mode, this is the summed numer of reads mapping to a gene across all cells in a predefined cluster.
- `--conditionA` and `--conditionB`: In order to define a difference in 𝜓 across conditions, the direction of comparison must be defined. Delta 𝜓 for each gene is defined as the mean 𝜓 value in condition B minus the mean 𝜓 value in condition A. Both `conditionA` and `conditionB` must be found in the condition column of the `conditions` file.

### Expected output of LABRATsc

Following quantification, 𝜓 values for all genes in all conditions as well as raw and Benjamini-Hochberg corrected p-values are reported in files named "results.subsampleclusters.txt" (`subsampleClusters` mode) or "results.cellbycell.txt" (`cellbycell` mode). Differences in mean 𝜓 values across conditions are also reported. In `cellbycell` mode, the results file additionally includes the number of cells in each condition passing read depth filters for each gene. Finally, per-cell psi values are reported for each gene when run in `cellbycell` mode in a file called "psis.cellbycell.txt.gz". These results can be sued to annotate existing single cell analyses.

