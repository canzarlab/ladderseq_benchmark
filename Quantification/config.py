# Config file for entire ladder kallisto simulation and analysis
# Authour: Shounak Chakraborty
# 07.08.2019

N_THREADS = 32

# config data path from simulation
CONFIG_DATA = '../configData/'

ANNO_PREFIX = 'Homo_sapiens.GRCh38.91'
ANOTATION_GTF_FILE = CONFIG_DATA+ANNO_PREFIX+'.gtf'
ANOTATION_FASTA_FILE = CONFIG_DATA+ANNO_PREFIX+'.fa'
ANNOTATION_RSEM = CONFIG_DATA+'RSEMRef/rsemReference'

GENOME_FASTA_FILE = CONFIG_DATA+'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
TRANS_LENGTH_FILE = CONFIG_DATA+'transLengthFile.txt'

REAL_DATA_BASE_PATH = CONFIG_DATA+'realData/'
REAL_DATA_FILE_NAME = 'NA12716_7'
REAL_DATA_FILE_TYPE = '.fastq.gz'
#SIM_DATA_BASE_PATH = '../simulations/sim_30000000'
BAND_PROB_FILE_PATH = CONFIG_DATA+'bandProbabilities/'


GENOME_BWT = CONFIG_DATA+'bowtie/bowtieGenome'
ANNOTATION_BWT = CONFIG_DATA+'bowtie/bowtieAnnotation'

KAL_IDX = CONFIG_DATA+'kallisto_index.idx'


## Softwares
KALLISTO_LS_BINARY = './ext/linuxBinaries/kallisto-ls'


# total number of simulations
SIMS_LIST = [str(i) for i in range(1, 3)]
NUM_SIM = 2
# Number of reads
N_READS = "30000000"
SIM_DATA_BASE_PATH = '../simulations/sim_'+N_READS
BETA_ESTIMATION = 'TRUE'
BANDS = ["band1","band2","band3","band4","band5","band6","band7"]
BAND_PROB = ["betas_realistic","betas_perfect"]
ANALYSIS_TYPE=["geneComp"]
