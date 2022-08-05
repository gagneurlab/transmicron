## TransMiCron: Accurate prediction of insertion probabilities improves detection of cancer driver genes from transposon mutagenesis screens
test transposon


![Overview of the Transmicron method](transmicron_method.png)


## Installation & Dependencies
The code is written in R and python and tested on python 3.9 and a Linux OS but should run on any OS. A [Conda environment](environment.yml) is provided which contains all the required libraries to run the code and needs to be activated before attempting to run. 
Using the following conda command, you can create an environment named "transmicron" with specified libraries in environment.yml and then activate the environment:

```
conda env create --name transmicron --file=environment.yml
conda activate transmicron
```


## Running the Code
The project can be run via [Snakemake](https://snakemake.readthedocs.io/en/stable/). The Snakemake workflow management system is a tool to create reproducible and scalable data analyses. To run the whole workflow, some variables in [Snakemake Config file](config.yaml) and [config file](config.yaml) need to be changed. Depending on how you wanna run the code, you can run the code with one of the following settings:




### Running the whole pipeline
Set the following configurations.
 * In [Snakemake Config file](config.yaml):
     * Set mutagenesis_method: ["predefinedFeatures"]
     * Set annotation: ["genes", "10kb", "20kb"...] # genes or any bins
     * Supply the Dataset and InsertionFile
 

### Running with precomputed model
Set the following configurations.
 * In [Snakemake Config file](config.yaml):
     * Set mutagenesis_method: ["pretrainedModel"]
After setting the proper configurations and activating previously created conda environment, start running the code by the following commands:
```
snakemake -n #dry run
snakemake --cores 128
```


lkhjdslkjflk

It is possible to also invoke single workflows explicitly e.g. for prepareAnnoation with:
```
snakemake prepareAnnoation --cores 10 # run with 10 cores
```
TODO: make a script to take do this, instead of directly changing config and snakemake file.

## Datasets
Following datasets are used in the pipeline and provided in Input folder.
* Input/FeatureBeds/ATACmesc_SRX2514792.bed
* Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed
* Input/Annotations/genes/mm10.refGene.gtf.gz 

## Acknowledgements and Funding
