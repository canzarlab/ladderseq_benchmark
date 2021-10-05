### Snakemake config file for Running Stringtie_original and Stringtie_modified on samples and analysing the results. ###
### Author: Shounak ###
### 11.05.2018 ###



# softwares
KALLISTO_LS_BINARY = './ext/linuxBinaries/kallisto-ls'
STAR_BINARY = './ext/linuxBinaries/STAR'
RSEM_BINARY_PATH = './ext/linuxBinaries/rsemBinaries'
BOWTIE2_BINARY_PATH = './ext/linuxBinaries/bowtie2Binaries'
SAMTOOLS_BINARY = './ext/linuxBinaries/samtools'
EXONREFINE_BINARY = './ext/linuxBinaries/exonRefine'
GROUPGENES_BINARY = './ext/linuxBinaries/groupGenes'
GTFFILTER_BINARY = './ext/linuxBinaries/gtfFilter'
GFFREAD_BINARY = './ext/linuxBinaries/gffread'
BEDTOOLS_BINARY = './ext/linuxBinaries/bedtools'
STRINGTIE2_BINARY = './ext/linuxBinaries/stringtie'
GFFCOMPARE_BINARY = './ext/linuxBinaries/gffcompare'




CONFIG_DATA = '../configData/'


# simulation configurations
ANNO_PREFIX = 'Homo_sapiens.GRCh38.91'
ANOTATION_GTF_FILE = CONFIG_DATA+ANNO_PREFIX+'.gtf'
ANOTATION_FASTA_FILE = CONFIG_DATA+ANNO_PREFIX+'.fa'
ANNOTATION_RSEM = CONFIG_DATA+'RSEMRef/rsemReference'

GENOME_FASTA_FILE = CONFIG_DATA+'Homo_sapiens.GRCh38.dna.primary_assembly.fa'
GENEID_TRANSID_FILE = CONFIG_DATA+'geneId_transId.txt'
TRANS_LENGTH_FILE = CONFIG_DATA+'transLengthFile.txt'

REAL_DATA_BASE_PATH = CONFIG_DATA+'realData/'
REAL_DATA_FILE_NAME = 'NA12716_7'
REAL_DATA_FILE_TYPE = '.fastq.gz'

BAND_PROB_FILE_PATH = CONFIG_DATA+'bandProbabilities/'

GENOME_BWT = CONFIG_DATA+'bowtie/bowtieGenome'
ANNOTATION_BWT = CONFIG_DATA+'bowtie/bowtieAnnotation'

STAR_GENOME_DIR = CONFIG_DATA+"starGenomeDirectory"

N_READS = "30000000"
SIM_DATA_BASE_PATH = '../simulations/sim_'+N_READS
SIMS_LIST = [str(i) for i in range(1, 2)]
NUM_SIM = 1
BAND_PROB = ["betas_realistic","betas_increase_1", "betas_increase_4", "betas_increase_7","betas_perfect"]
#BAND_PROB = ["betas_realistic"]

KAL_IDX = CONFIG_DATA+'kallisto_index.idx'

# Ground truth configurations
ABUNDANCE_FILE = REAL_DATA_BASE_PATH+"rsem/out.isoforms.results"
GROUND_TRUTH_BASE_PATH =  '../GroundTruth'
## the complexity divided ground truth files
GROUND_TRUTH_PATH = GROUND_TRUTH_BASE_PATH+'/tpmFilter0.1/CoverageFilter'
COMPLETE_GROUND_TRUTH_FILE_PATH = GROUND_TRUTH_BASE_PATH+'/complexityFilter'
SIMULATION_TYPE = 'RSEM'
READ_LENGTH = 75
OUTPUT_TYPE_READtoBED = 'NotOnly_transcript'
WITHOUT_ORIGINAL_TRANSCRIPTS = 'F'
DESIRED_TPM = "0.1"
ALL_COMPLEXITIES = list(range(1,12))
LAST_COMPLEXITY = 11





## This has to range from 1 to the total number of Transcript complexity in the ground truth file
COMPLEXITIES = list(range(1,11))

# Ref based assembly configurations
INDIVIDUAL_BAND = ["band1","band2","band3","band4","band5","band6","band7"]
COMBINED_BAND = ["comb12","comb23","comb34","comb45","comb56","comb67"]
N_THREADS = 8
DATASET = "30Mil"
FILTER_TYPE_INDIVIDUAL="CORRECT_PLUS_ONE"
FILTER_TYPE_COMBINED="CORRECT_PLUS_ONE"
PAIRWISE = 'TRUE'
STRINGTIE_LOW_THRESHOLD = 'FALSE'
LADDER_RESOLUTION = "HIGH"
LADDER_HIGH_FILTER_THRESHOLD = 1
POLY_A = 200
INDIVIDUAL_STRINGTIE_THRESHOLD = 1
MOD_KALLISTO_FILTER = 0.1
SUBEXON_LENGTH_THRESHOLD = 2
RESTORE = 'TRUE'
RESTORED_TPM_FILTER = 1
