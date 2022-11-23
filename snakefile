configfile: "./user_config.yaml"


#####################################################################################
# memory management @ata: please feel free to optimize ressource management in whatever way works best
#####################################################################################

#dictionary to map dataset name to insertion BED-file
insertionFilesZIP = zip(config["datasets"], config["insertionFile"])
insertionFilesDict = dict(insertionFilesZIP)
#add unselected datasets to dictionary (we use these datasets to calculate precomputed insertion rates)
insertionFilesDict["unselectedPB"] = "Input/BEDInsertionTesting/PBmESC.BED"
insertionFilesDict["unselectedSB"] = "Input/BEDInsertionTesting/SBmESC.BED"


#ressources needed to build the feature matrices for all TTAA / TA sites in the genome
model_resources = {
    "PB": {"mem_mb": 60000},
    "SB": {"mem_mb": 250000}}
    
glm_resources = {
    "AML": {"mem_mb": 60000},
    "TALL": {"mem_mb": 60000},
    "DLBCLPB": {"mem_mb": 30000},
    "DLBCLSB":{"mem_mb": 30000},
    "Pancreas": {"mem_mb": 250000},
    "Melanoma": {"mem_mb": 150000},
    "colonRAD": {"mem_mb": 60000},
    "Cuka": {"mem_mb": 60000},
    "SBdriver_ML": {"mem_mb": 160000},
    "wholeBodyHealthy":{"mem_mb": 60000},
    "DLBCL_PBtail": {"mem_mb": 60000},
    "DLBCL_SBtail": {"mem_mb": 60000},
    "unselected": {"mem_mb": 60000},
    "PancreasNoCancer": {"mem_mb": 260000},
    "unselectedPB": {"mem_mb": 60000},
    "unselectedSB": {"mem_mb": 60000},
    "Liver": {"mem_mb": 260000},
    "ControlsSB": {"mem_mb": 60000},
    "ControlsPB": {"mem_mb": 60000},
    "PRIM_dataset": {"mem_mb": 160000},
    "Intestinal": {"mem_mb": 250000},
    "unselectedYohsibaPBIllumina": {"mem_mb": 60000},
    "cuSCC": {"mem_mb": 60000}}
    

transposon_system = {
    "AML": {"transposon_system": "PB"},
    "TALL": {"transposon_system": "PB"},
    "DLBCLPB": {"transposon_system": "PB"},
    "DLBCLSB":{"transposon_system": "SB"},
    "Pancreas": {"transposon_system": "SB"},
    "Melanoma": {"transposon_system": "SB"},
    "colonRAD": {"transposon_system": "PB"},
    "Cuka": {"transposon_system": "SB"},
    "SBdriver_ML": {"transposon_system": "SB"},
    "wholeBodyHealthy":{"transposon_system": "PB"},
    "DLBCL_PBtail": {"transposon_system": "PB"},
    "DLBCL_SBtail": {"transposon_system": "SB"},
    "unselected": {"transposon_system": "PB"},
    "PancreasNoCancer": {"transposon_system": "SB"},
    "unselectedSB": {"transposon_system": "SB"},
    "unselectedPB": {"transposon_system": "PB"},
    "ControlsSB": {"transposon_system": "SB"},
    "Liver": {"transposon_system": "SB"},
    "ControlsPB": {"transposon_system": "PB"},
    "PRIM_dataset": {"transposon_system": "SB"},
    "Intestinal": {"transposon_system": "SB"},
    "unselectedYohsibaPBIllumina":  {"transposon_system": "PB"},
    "cuSCC":  {"transposon_system": "SB"}}

 
#####################################################################################
# create target rule
#####################################################################################

#loop through "datasets" and "transposonSystem" in config-file
input_target_rule=list()
for i in range(0, len(config["datasets"])):
    
    if("pretrainedModel" in config["mutagenesis_method"]):
        input_target_rule.append(expand("Output/{dataset}/{transposonSystem}/results/results_{annotation}_{insertionRates}_{mutaFeatures}.csv", dataset=config["datasets"][i], annotation=config["annotation"], transposonSystem=config["transposonSystem"][i], insertionRates = "precompInsRates",mutaFeatures = "precompFeatures" ))

    if("predefinedFeatures" in config["mutagenesis_method"]):
        input_target_rule.append(expand("Output/{dataset}/{transposonSystem}/results/results_{annotation}_{insertionRates}_{mutaFeatures}.csv", dataset=config["datasets"][i], annotation=config["annotation"], transposonSystem=config["transposonSystem"][i], insertionRates = "NewInsRates",mutaFeatures = "precompFeatures" ))

    if("noMutagenesis" in config["mutagenesis_method"]):
        input_target_rule.append(expand("Output/{dataset}/{transposonSystem}/results/results_{annotation}_{insertionRates}_{mutaFeatures}.csv", dataset=config["datasets"][i], annotation=config["annotation"], transposonSystem=config["transposonSystem"][i], insertionRates = "EqualWeightInsRates",mutaFeatures = "precompFeatures"))

    if("OnlyCustomFeatures" in config["mutagenesis_method"]):
        input_target_rule.append(expand("Output/{dataset}/{transposonSystem}/results/results_{annotation}_{insertionRates}_{mutaFeatures}.csv", dataset=config["datasets"][i], annotation=config["annotation"], transposonSystem=config["transposonSystem"][i], insertionRates = "NewInsRates",mutaFeatures = "OnlyCustomFeatures"))

    if("AddCustomFeatures" in config["mutagenesis_method"]):
        input_target_rule.append(expand("Output/{dataset}/{transposonSystem}/results/results_{annotation}_{insertionRates}_{mutaFeatures}.csv", dataset=config["datasets"][i], annotation=config["annotation"], transposonSystem=config["transposonSystem"][i], insertionRates = "NewInsRates",mutaFeatures = "AddCustomFeatures"))

rule target:
    input: 
        input_target_rule

    
##########################################################################
#Step 0: Preprocess user input
##########################################################################
    
# 0.1.: preprocess insertion data (BED-file)



rule processInsertionFile:
    params:
        insertionFile = lambda wildcards: insertionFilesDict[wildcards.dataset]
    output:
        Granges_insertions = "Output/{dataset}/{transposonSystem}/Granges_insertions.RData"
    script:
        "Scripts/processInsertionBED.R"
        
# 0.2.: create table of features for Step 1 (mutagenesis). If users supply custom features these are included here.

#dictionary to assign path for the feature tables (we supply precomputed files for our predefined features)
pathNegatives = {
    "precompFeatures" : "Input", #if predefined features: use precomputed files
    "OnlyCustomFeatures" : "{dataset}", #else: recompute tables
    "AddCustomFeatures" : "{dataset}"
}


rule defineMutagenesisFeatures:
    params:
        REFSEQgtf = "Input/Annotations/genes/mm10.refGene.gtf.gz",
        DNAse_mESC = "Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed", #pre-supplied DNAse-seq data
        ATAC_mESC = "Input/FeatureBeds/ATACmesc_SRX2514792.bed",
        customFeatures = config["customFeatures"],#list of Bed-files, user might supply to control for
    output:
        featureList = "Input/{transposonSystem}/InputFeaturesMutagenesis/featureList_{mutaFeatures}.RData"
    script:
        "Scripts/MutagenesisModel/createInputFeatures.R"
    

# 0.3.: prepare annotation according to users choice (creates GRanges Objects)

rule prepareAnnoation:
    params:
        CustomAnnotation = config["customAnnotation"],# if the user has supplied a custom annotation, this is fed in here; @ata: maybe you find a better way to supply user defined files as inputs?
        REFSEQgtf = "Input/Annotations/genes/mm10.refGene.gtf.gz"
    output:
        AnnotationGRanges = "Input/Annotations/{annotation}/AnnotationGRanges.RData"
    script:
        "Scripts/Annotations/PrepareAnnotationGRanges.R"
        
        
##########################################################################
#Step 1: Mutagenesis Model
##########################################################################

# 1.1.: Create Model Input

#create PRECALCULATED feature table with all TTAA / TA sites in the genome
rule GenomeTargetSitesPrecomp: 
    input:
        featureList = "Input/{transposonSystem}/InputFeaturesMutagenesis/featureList_{mutaFeatures}.RData"
    params:
        InputFunctions = workflow.source_path("Scripts/MutagenesisModel/define_functions.R"), #script containing functions to create input matrices 
        mutagenesis_method = config["mutagenesis_method"],
        REFSEQgtf = "Input/Annotations/genes/mm10.refGene.gtf.gz"
    resources:
        mem_mb = lambda wildcards: model_resources[wildcards.transposonSystem]['mem_mb']
    output:
        one_hot_encoder = "Input/{transposonSystem}/InputFeaturesMutagenesis/one_hot_encoder_{mutaFeatures}.RData",
        granges_negatives_controls = "Input/{transposonSystem}/InputFeaturesMutagenesis/granges_negatives_controls_{mutaFeatures}.RData",
        step1_input_negatives_controls = "Input/{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_controls_{mutaFeatures}.RData",
        one_bp_ranges_negatives ="Input/{transposonSystem}/InputFeaturesMutagenesis/one_bp_ranges_negatives_{mutaFeatures}.RData",
        step1_input_negatives_whole_genome = "Input/{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_whole_genome_{mutaFeatures}.h5"
    script:
        "Scripts/MutagenesisModel/GenomeTargetSites.R"

        
#reaclculate the feature table with all TTAA / TA sites in the genome (needed if user supplies their own features)
rule GenomeTargetSitesNew: 
    input:
        featureList ="Input/{transposonSystem}/InputFeaturesMutagenesis/featureList_{mutaFeatures}.RData"
    params:
        InputFunctions = workflow.source_path("Scripts/MutagenesisModel/define_functions.R"), #script containing functions to create input matrices 
        mutagenesis_method = config["mutagenesis_method"],
        REFSEQgtf = "Input/Annotations/genes/mm10.refGene.gtf.gz"
    resources:
        mem_mb = lambda wildcards: model_resources[wildcards.transposonSystem]['mem_mb']
    output:
        one_hot_encoder = "Output/{dataset}/{transposonSystem}/InputFeaturesMutagenesis/one_hot_encoder_{mutaFeatures}.RData",
        granges_negatives_controls = "Output/{dataset}/{transposonSystem}/InputFeaturesMutagenesis/granges_negatives_controls_{mutaFeatures}.RData",
        step1_input_negatives_controls = "Output/{dataset}/{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_controls_{mutaFeatures}.RData",
        one_bp_ranges_negatives ="Output/{dataset}/{transposonSystem}/InputFeaturesMutagenesis/one_bp_ranges_negatives_{mutaFeatures}.RData",
        step1_input_negatives_whole_genome = "Output/{dataset}/{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_whole_genome_{mutaFeatures}.h5"
    script:
        "Scripts/MutagenesisModel/GenomeTargetSites.R"
        
ruleorder: GenomeTargetSitesPrecomp > GenomeTargetSitesNew

#create input table to train mutagenesis model:

rule MutagenesisModelInput:
    input:
        featureList ="Input/{transposonSystem}/InputFeaturesMutagenesis/featureList_{mutaFeatures}.RData",
        Granges_insertions = "Output/{dataset}/{transposonSystem}/Granges_insertions.RData",
        one_hot_encoder = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures],"{transposonSystem}/InputFeaturesMutagenesis/one_hot_encoder_{mutaFeatures}.RData"]),
         granges_negatives_controls = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures],"{transposonSystem}/InputFeaturesMutagenesis/granges_negatives_controls_{mutaFeatures}.RData"]),
         step1_input_negatives_controls = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures], "{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_controls_{mutaFeatures}.RData"])
    params:
        InputFunctions = workflow.source_path("Scripts/MutagenesisModel/define_functions.R"), #script containing functions to create input matrices 
        transposon_system = lambda wildcards: transposon_system[wildcards.dataset]['transposon_system'],
        REFSEQgtf = "Input/Annotations/genes/mm10.refGene.gtf.gz"
    resources:
        mem_mb=160000
    output:
        step1_input_positives_RAW ="Output/{dataset}/{transposonSystem}/step1_input_positives_RAW_{mutaFeatures}.RData",
        step1_input_positives = "Output/{dataset}/{transposonSystem}/step1_input_positives_{mutaFeatures}.RData",
        MutagenesisTrainingData = "Output/{dataset}/{transposonSystem}/MutagenesisTrainingData_{mutaFeatures}.RData",
        MutagenesisTestData = "Output/{dataset}/{transposonSystem}/MutagenesisTestData_{mutaFeatures}.RData"
    script:
        "Scripts/MutagenesisModel/MutagenesisModelInput.R"
        
        
# 1.2.: Train Mutagenesis Model

rule trainMutagenesisModel:
    input:
        MutagenesisTrainingData = "Output/{dataset}/{transposonSystem}/MutagenesisTrainingData_{mutaFeatures}.RData" # @ata: we could change how we transfer files between R and python from .RData files to .csv; alternatively, we could dimplement the model in r
    output:
        MutagenesisModel = "Output/{dataset}/{transposonSystem}/MutagenesisModel_{mutaFeatures}.pickle"
    resources:
        mem_mb=60000
    script:
        "Scripts/MutagenesisModel/MutagenesisModelTaining.py"
        
#1.3.: Use Mutagenesis Model to predict insertion rates for all TTAA / TA sites in the genome

rule predictInsertionRates:
    input:
        predict_insertion_model = "Output/{dataset}/{transposonSystem}/MutagenesisModel_{mutaFeatures}.pickle",
        step1_input_negatives_whole_genome = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures], "{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_whole_genome_{mutaFeatures}.h5"]),
        step1_input_positives = "Output/{dataset}/{transposonSystem}/step1_input_positives_{mutaFeatures}.RData"
    resources:
        mem_mb=100000
    output:               
        InsertionRatesGenome = "Output/{dataset}/{transposonSystem}/InsertionRates/Genomewide_{insertionRates}_{mutaFeatures}.csv"
    script:
        "Scripts/MutagenesisModel/predictInsertionRatesGenome.py"

# if not mutagenesis correction is used, this rule creates a vector with equals weights for all TTAA / TA sites in the genome:




rule NULL_model_predictions:
    input:
        step1_input_negatives_whole_genome = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures], "{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_whole_genome_{mutaFeatures}.h5"])
    resources:
        mem_mb=100000
    output:               
        InsertionRatesGenome = "Output/{dataset}/{transposonSystem}/InsertionRates/Genomewide_EqualWeightInsRates_{mutaFeatures}.csv"
    script:
        "Scripts/MutagenesisModel/NoMutagenesisInsertionRates.py"
     
#create PRECOMPUTED insertion rates and annotation matrices; this rule will not be executed if our precomputed files are present
#the precomputed rates are calculated based on unselected datasets

#create dictionary to map transposon system to the appropriate unselected dataset
unselectedData = {
    "PB": "unselectedPB",
    "SB": "unselectedSB"}
    

rule predictInsertionRatesPRECOMPUTED:
    input:
        predict_insertion_model = lambda wildcards: "/".join([unselectedData[wildcards.transposonSystem], "{transposonSystem}/MutagenesisModel_precompFeatures.pickle"]),
        step1_input_negatives_whole_genome = "Input/{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_whole_genome_precompFeatures.h5",
        step1_input_positives = lambda wildcards: "/".join([unselectedData[wildcards.transposonSystem], "{transposonSystem}/step1_input_positives_precompFeatures.RData"])
    resources:
        mem_mb=100000
    output:               
        InsertionRatesGenome = "Input/{transposonSystem}/InsertionRates/Genomewide_precompInsRates_precompFeatures.csv"
    script:
        "Scripts/MutagenesisModel/predictInsertionRatesGenome.py"

rule NULL_model_predictionsPRECOMPUTED:
    input:
        step1_input_negatives_whole_genome = "Input/{transposonSystem}/InputFeaturesMutagenesis/step1_input_negatives_whole_genome_precompFeatures.h5"
    resources:
        mem_mb=100000
    output:               
        InsertionRatesGenome = "Input/{transposonSystem}/InsertionRates/Genomewide_EqualWeightInsRates_precompFeatures.csv"
    script:
        "Scripts/MutagenesisModel/NoMutagenesisInsertionRates.py"

ruleorder: NULL_model_predictionsPRECOMPUTED > predictInsertionRatesPRECOMPUTED > NULL_model_predictions > predictInsertionRates

# for debugging
#rule precompAnnotations:
 #   output:
  #      AnnotationEqualWeight = "Input/{transposonSystem}/{annotation}/AnnotationMatrix_EqualWeightInsRates_precompFeatures.RData",
   #     AnnotationPrecomputed = "Input/{transposonSystem}/{annotation}/AnnotationMatrix_precompInsRates_precompFeatures.RData"
        
########################################################################
#Step 2: Selection Model
########################################################################

# 2.1.: Create input for the selection model

#dictionatry to match insertion rates to file paths depending on wildcard, i.e. user choice (we supply precomputed insertion rates)

pathInsertion = {
    "precompInsRates": "Input",
    "NewInsRates": "{dataset}",
    "EqualWeightInsRates" : "Input"
}

rule matrixInsRatesPrecomp:
    input:
        InsertionRatesGenome = "Input/{transposonSystem}/InsertionRates/Genomewide_precompInsRates_precompFeatures.csv",
        #InsertionRatesGenome = "Input/{transposonSystem}/InsertionRates/Genomewide_precompInsRates_precompFeatures.csv",
        one_bp_ranges_negatives = "Input/{transposonSystem}/InputFeaturesMutagenesis/one_bp_ranges_negatives_precompFeatures.RData",
        AnnotationGRanges = "Input/Annotations/{annotation}/AnnotationGRanges.RData"
    output:
        AnnotationMatrix= "Input/{transposonSystem}/{annotation}/AnnotationMatrix_precompInsRates_precompFeatures.RData"
    resources:
        mem_mb=60000
    script:
        "Scripts/Annotations/AnnotationInsertionRates.R"


        
rule matrixInsRatesNew:
    input:
        InsertionRatesGenome = "Output/{dataset}/{transposonSystem}/InsertionRates/Genomewide_{insertionRates}_{mutaFeatures}.csv",
        one_bp_ranges_negatives = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures],"{transposonSystem}/InputFeaturesMutagenesis/one_bp_ranges_negatives_{mutaFeatures}.RData"]),
        AnnotationGRanges = "Input/Annotations/{annotation}/AnnotationGRanges.RData"
    output:
        AnnotationMatrix= "Output/{dataset}/{transposonSystem}/{annotation}/AnnotationMatrix_{insertionRates}_{mutaFeatures}.RData"
    resources:
        mem_mb=60000
    script:
        "Scripts/Annotations/AnnotationInsertionRates.R"

ruleorder: matrixInsRatesPrecomp > matrixInsRatesNew
        
        
rule createInputSelectionModel:
    input:
        annotationMatrix = lambda wildcards: "Output/" + "/".join([pathInsertion[wildcards.insertionRates],"{transposonSystem}/{annotation}/AnnotationMatrix_{insertionRates}_{mutaFeatures}.RData"]),
        AnnotationGRanges = "Input/Annotations/{annotation}/AnnotationGRanges.RData",
        InsertionRatesGenome =  lambda wildcards: "Output/" + "/".join([pathInsertion[wildcards.insertionRates],"{transposonSystem}/InsertionRates/Genomewide_{insertionRates}_{mutaFeatures}.csv"]),
        Granges_insertions = "Output/{dataset}/{transposonSystem}/Granges_insertions.RData",
        one_bp_ranges_negatives = lambda wildcards: "/".join([pathNegatives[wildcards.mutaFeatures],"{transposonSystem}/InputFeaturesMutagenesis/one_bp_ranges_negatives_{mutaFeatures}.RData"])
    params:
        transposon_system = lambda wildcards: transposon_system[wildcards.dataset]['transposon_system']
    output:
        SelectionModelInput = "Output/{dataset}/{transposonSystem}/{annotation}/SelectionModelInput_{insertionRates}_{mutaFeatures}.RData"
    resources:
        mem_mb=30000
    script:
        "Scripts/SelectionModel/createSelectionModelInput.R"

# 2.2.: Run selection model and generate output
rule runSelectionModel:
    input:
        SelectionModelInput50kb = "Output/{dataset}/{transposonSystem}/50kb/SelectionModelInput_{insertionRates}_{mutaFeatures}.RData", #used to fit the first glm 
        SelectionModelInputTargetAnnotation = "Output/{dataset}/{transposonSystem}/{annotation}/SelectionModelInput_{insertionRates}_{mutaFeatures}.RData" #used to call CIS
    output:
        results = "Output/{dataset}/{transposonSystem}/results/results_{annotation}_{insertionRates}_{mutaFeatures}.csv"
    params:
        multest_correction = config["multest_correction"]
    resources:
        mem_mb = lambda wildcards: glm_resources[wildcards.dataset]['mem_mb']
    script:
        "Scripts/SelectionModel/runSelectionModel.R"
