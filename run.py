import yaml
import os
from argparse import ArgumentParser, Namespace

def get_args() -> Namespace:
    """Parses given arguments
    Returns:
        Namespace: parsed arguments
    """

    parser = ArgumentParser(description="parameters")
    parser.add_argument("--datasets", default='DLBCLPB', type=str,)
    parser.add_argument("--insertionFile", default='Input/TestInsertionBed/DLBCLPB.BED', type=str,)
    parser.add_argument("--transposonSystem", default='PB', type=str,)
    parser.add_argument("--mutagenesisMethod", default='predefinedFeatures', type=str,)
    parser.add_argument("--annotation", default='genes', type=str,)
    parser.add_argument("--CustomAnnotation", default='Input/testAnnot.BED', type=str,)
    parser.add_argument("--multestCorrection", default='bonferroni', type=str,)
    parser.add_argument("--promoters", default=0, type=int, )
    parser.add_argument("--snakemakeRule", default="-n", type=str)
    parser.add_argument("--customFeatures", default="Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed,Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed", type=str) 
    parser.add_argument("--outputDir", default="./Output", type=str,)
    return parser.parse_args()


args = get_args()
args = vars(args)
args['datasets'] = args['datasets'].split(',')
args['insertionFile'] = args['insertionFile'].split(',')
args['transposonSystem'] = args['transposonSystem'].split(',')
args['multest_correction'] = args['multestCorrection'].split(',')
args['mutagenesis_method'] = args['mutagenesisMethod'].split(',')
args['annotation'] = args['annotation'].split(',')
args['customFeatures'] = args['customFeatures'].split(',')



with open('./user_config.yaml', 'w') as file:
    documents = yaml.dump(args, file, default_flow_style=False)


#os.system("./run_slurm.sh" + " " + args['snakemake_rule'])
os.system("snakemake " + args['snakemakeRule'])
