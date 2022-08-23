## TransMiCron: Accurate prediction of insertion probabilities improves detection of cancer driver genes from transposon mutagenesis screens

## Background: Transmicron
Transmicron is an easy to run pipeline to identify common insertion sites (CIS) from transposon mutagenesis screens. It comprises two steps: A) The mutagenesis model to correct for transposon integration biases, and B) the selection model to identify CIS. In the following, we provide a guide on how to run Transmicron. Details on the method can be found in the manuscript [link]. 

![Overview of the Transmicron method](transmicron_method.png)

##Running transmicron

### Installation & Dependencies
1. Make sure anaconda is installed on your system (https://docs.anaconda.com/anaconda/install/).
2. @ata: please help, something like: Clone the repository / Download the folder, open a terminal and navigate to the Transmicron folder.
3. Create and activate the [Conda environment](environment.yml) using the following comands:

```
conda env create --name transmicron --file=environment.yml
conda activate transmicron
```

### Prepare insertion data
YOU CAN SKIP THIS STEP IF YOU RUN TRANSMICRON ON THE TEST DATA PROVIDED BY US.

4. Please prepare a BED file specifying the genomic locations of your transposon insertions. 5 columns are mandatory:
Column 1 should specify the chromosome (Format: chr1, chr2...chrX, chrY).
Column 2 specifies the START of the insertion sites (a number).
Column 3 specifies the END of the insertion sites (a number). PLEASE ENSURE THAT START AND END ARE CONSECUTIVE NUCLEOTIDES (e.g. Start: 1234, End: 1235).
Column 4 specifies the orientation of transposon insertions ("-" or "+").
Column 5 supplies sample / Tumour identifiers.

Column headers for columns 4 and 5 should be "orientation" and "TumorID" (other column names do not matter). Additional columns will be disregarded. An example can be found under Input/TestInsertionBed/. 

### Executing the code
5. Execute the [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to run Transmicron using the following commands:
```
snakemake -n #dry run
snakemake --cores 8
```
You can add several options to these commands, depending on how you want to run Transmicron. If you want to use default parameters and data, you do not need to specify any additional parameters:
-d Please supply the filepath to the BED file containing the locations of your insertions. Do not specify if you want to run Transmicron on the testing dataset.
-t Which transposon system was used? Options: 1. PB [the default, PiggyBac] 2. SB [Sleeping Beauty]
-m Which version of the mutagenesis model do you want to apply? Options: 1. pretrained [the default; a pretrained version of the mutagenesis model on mESC insertions is used] 2. retrain [the mutagenesis model is retrained on your data using our predefined features] 3. null [no mutagenesis model is applied, Transmicron contols only for the distribution of TA / TTAA nucleotides] 4. filepath [please supply one or more filepaths to BED files containing feature information. The mutagenesis model is retrained using the distance of insertions to these features as input.]
-a Please Specify the target annotation used to identify CIS. Options: 1. genes [the default] 2. 10kb, 20kb...[any binlength] 3. filepath [specify a filepath to a BED file file custom features of interest, e.g. regulatory elementes]


It is possible to also invoke single workflows explicitly e.g. for prepareAnnoation with:
```
snakemake prepareAnnoation --cores 10 # run with 10 cores
```

The code is written in R and python and tested on python 3.9 and a Linux OS but should run on any OS. 

TODO: make a script to take do this, instead of directly changing config and snakemake file.

## Required Ressources
The ressources required strongly depend on which parameters you choose above. On the test data, using default parameters, roghly XX mb of memory are reuiqred. The pipeline takes roughly xx minutes to run. 

## Features and Annotations
By default the following datasets are used:

Insertion Sites: 
* Input/TestInsertionBed/DLBCLPB.BED
Input features for the mutagenesis model:
* Input/FeatureBeds/ATACmesc_SRX2514792.bed
* Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed

Gene annotation:
* Input/Annotations/genes/mm10.refGene.gtf.gz 

All files are provided in Input folder.

## Acknowledgements and Funding
