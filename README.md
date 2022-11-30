## TransMiCron: Accurate prediction of insertion probabilities improves detection of cancer driver genes from transposon mutagenesis screens

## Background: Transmicron
Transmicron is an easy to run pipeline to identify common insertion sites (CIS) from transposon mutagenesis screens. It comprises two steps: A) The mutagenesis model to correct for transposon integration biases, and B) the selection model to identify CIS. In the following, we provide a guide on how to run Transmicron. Details on the method can be found in the manuscript [link]. 

![Overview of the Transmicron method](transmicron_method.png)

##Running transmicron

### Installation & Input Data
1. Make sure anaconda is installed on your system (https://docs.anaconda.com/anaconda/install/).
2. Clone or download the repository, then open a terminal and navigate to the Transmicron folder.
```
git clone https://github.com/gagneurlab/transmicron.git
cd transmicron
```
3. Create and activate the [Conda environment](environment.yml) using the following comands:

```
conda env create --name transmicron --file=environment.yml
conda activate transmicron
```
4. (Optional, Recommended) If you want to use our precomputed features, download all Input and features from projects [zenodo](https://zenodo.org/record/7373066) using the following code. The files are about 15GB, so this might take some time.
```
wget https://zenodo.org/record/7373066/files/Input.zip?download=1
unzip -o Input.zip?download=1 && rm Input.zip?download=1.zip
snakemake --cores 1 --touch
```
If you don't download the data from zenodo, you can still run the program with the Input folder provided in this github repository, but all features will be calculated from scratch.
### Prepare insertion data
YOU CAN SKIP THIS STEP IF YOU RUN TRANSMICRON ON THE TEST DATA PROVIDED BY US.

4. Please prepare a [BED](https://www.genomatix.de/online_help/help_regionminer/bedformat_help.html) file specifying the genomic locations of your transposon insertions. 5 columns are mandatory:
* Column 1 specifies the chromosome (Format: chr1, chr2...chrX, chrY).
* Column 2 specifies the position [START] of the insertion site (a number).
* Column 3 specifies the END of the insertion site (a number). Alywas: END = START + 1  (e.g. START: 1234, END: 1235).
* Column 4 specifies the orientation of transposon insertions ("-" or "+").
* Column 5 supplies Sample / Tumour identifiers.

Column headers for columns 4 and 5 should be "orientation" and "TumorID" (other column names do not matter). Additional columns will be disregarded. An example can be found under Input/TestInsertionBed/. 

### Executing the code
5. Run the whole Transmicron pipeline using [Snakemake](https://snakemake.readthedocs.io/en/stable/)  by the following commands:
```
python run.py --snakemakeRule="-n" #dry run
python run.py --snakemakeRule="--cores 8"
```

You can add several parameters to these commands, depending on how you want to run Transmicron. All the parameters are listed below, with their default value. The program sets the default value to any parameter that is not defined by the user.
```
python run.py   
		--insertionFile="Input/BEDInsertionTesting/DLBCLPB.BED"  # path to comma seperated name of BED insertion files, the final results will be saved to [outputDir/datset] for each dataset provided here, e.g. Output/DLBCLPB
		--transposonSystem="PB"  # transposon systems of given datasets, the same order as insertion BED files
		--mutagenesisMethod="predefinedFeatures" 
		--customFeatures="Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed,Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed" 
		--annotation="genes" 
		--multestCorrection="bonferroni" 
		--snakemakeRule="--cores 8" 
		--outputDir="Output" 
				
```
All paths can be given relative to the root of the program (i.e. where the snakefile is) or as absolute paths. If you want to use default parameters and data, you do not need to specify any additional parameters:
* insertionFile: Please supply the filepath to the BED file containing the locations of your insertions. Do not specify if you want to run Transmicron on the testing dataset.
* transposonSystem: Options: 1. PB [PiggyBac, default] 2. SB [Sleeping Beauty]
* mutagenesisMethod: Which version of the mutagenesis model do you want to apply? Options:
  * 1. "predefinedFeatures": default option; a pretrained version of the mutagenesis model on mESC insertions is used.
  * 2. "pretrainedModel": the mutagenesis model is retrained on your data using our predefined features.
  * 3. "null": no mutagenesis model is applied, Transmicron contols only for the distribution of TA / TTAA nucleotides.
  * 4. filepath: please supply one or more filepaths to BED files containing feature information. The mutagenesis model is retrained using the distance of insertions to these features as input.
* annotation: Please Specify the target annotation used to identify CIS. Options: 1. "genes" [default] 2. "10kb", "20kb"...[any binlength] 3. filepath [specify a filepath to a BED file file custom features of interest, e.g. regulatory elementes]
* snakemakeRule: Specify snakemake rule here
* outputDir: Path to where the results of analysis will be written to. 


It is also possible to also invoke single workflows explicitly e.g. for defineMutagenesisFeatures with:
```
python run.py --snakemakeRule="-n" #dry run
python run.py --snakemakeRule="--until defineMutagenesisFeatures --cores 8" # run with 8 cores
```

The code is written in R and python and tested on python 3.9 and a Linux OS but should run on any OS. 


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
This work was supported by the German Bundesministerium für Bildung und Forschung (BMBF) supported the study through the VALE (Entdeckung und Vorhersage der Wirkung von genetischen Varianten durch Artifizielle Intelligenz für LEukämie Diagnose und Subtyp-Identifizierung) project (031L0203B to CB, XQ and JG)
