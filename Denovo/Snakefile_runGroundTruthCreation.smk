### Snakemake file for ###
### 1. Concatenating the output files from stringtie ###

### Author: Shounak ###
### 23.07.2018 ###



rule exonrefine:
	input:
		ANOTATION_GTF_FILE
	output:
		GROUND_TRUTH_BASE_PATH+"/initialExonRefined.gtf",
	shell:
		"""
			{EXONREFINE_BINARY} -p {GROUND_TRUTH_BASE_PATH}/initialExonRefined {input[0]}
		"""
rule groupGenes:
	input:
		GROUND_TRUTH_BASE_PATH+"/initialExonRefined.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_annotation.gtf"
	shell:
		"""
			{GROUPGENES_BINARY} {input[0]} {output[0]}
		"""

rule filterByExpressed:
	input:
		ABUNDANCE_FILE,
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_annotation.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/expressedTranscripts.txt",
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_expressed.gtf"
	shell:
		"""
			awk 'NR!=1 && $6>0 {{print $1}}' {input[0]} > {output[0]}
			{GTFFILTER_BINARY} --list {output[0]} -m whiteT {input[1]} {output[1]}
		"""


rule run_create_FullyCovered_GroupedGenesFiles:
	input:
		GROUND_TRUTH_BASE_PATH+"/fullyCoveredTranscripts.txt",
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_expressed.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_fullyCovered.gtf"
	shell:
		"""
		{GTFFILTER_BINARY} -m whiteT -l {input[0]} {input[1]} {output[0]}
		"""


rule run_ReadsToBedFiles:
	input:
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_expressed.gtf",
		SIM_DATA_BASE_PATH+"/sim_1/R_1.fq",
		SIM_DATA_BASE_PATH+"/sim_1/R_2.fq",
	output:
		GROUND_TRUTH_BASE_PATH+"/readPosFile.bed",
		GROUND_TRUTH_BASE_PATH+"/transcriptPosFile.bed"
	shell:
		"""
		python ReadsToBedFiles.py {input[0]} {ABUNDANCE_FILE} {input[1]} {input[2]} {GROUND_TRUTH_BASE_PATH}/ {SIMULATION_TYPE} {READ_LENGTH} {OUTPUT_TYPE_READtoBED}
		"""


rule run_BEDTOOLS_Operations:
	input:
		readPosFile = GROUND_TRUTH_BASE_PATH+"/readPosFile.bed",
		transcriptPosFile = GROUND_TRUTH_BASE_PATH+"/transcriptPosFile.bed"
	output:
		GROUND_TRUTH_BASE_PATH+"/sortedReadPos.bed",
		GROUND_TRUTH_BASE_PATH+"/mergedSortedReadPos.bed",
		GROUND_TRUTH_BASE_PATH+"/intersectionFile.bed",
		GROUND_TRUTH_BASE_PATH+"/fullyCoveredTranscripts.txt"
	shell:
		"""
		{BEDTOOLS_BINARY} sort -i {input[0]} > {output[0]}
		{BEDTOOLS_BINARY} merge -i {output[0]} > {output[1]}
		{BEDTOOLS_BINARY} intersect -a {input[1]} -b {output[1]} -wa -f 1 > {output[2]}
		awk '{{print $1}}' {output[2]} > {output[3]}
		"""

rule run_Prepare_Initial_GTF_file:
	input:
		gtfFile = GROUND_TRUTH_BASE_PATH+"/groupedGenes_expressed.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/InitialGtf_withTransId.gtf"
	shell:
		"""
		grep transcript_id {input[0]} > {output[0]}
		"""


rule run_filterByComplexity:
	input:
		infile = GROUND_TRUTH_BASE_PATH+"/InitialGtf_withTransId.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/complexityFilter/{complexity}_iso_ground_truth.gtf"
	shell:
		"""
		mkdir -p {GROUND_TRUTH_BASE_PATH}/complexityFilter
		{GTFFILTER_BINARY} -li {wildcards.complexity} -ui {wildcards.complexity} {input.infile} {output[0]}
		"""

rule run_filterByComplexity_lastComplexity:
	input:
		infile = GROUND_TRUTH_BASE_PATH+"/InitialGtf_withTransId.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/complexityFilter/"+str(LAST_COMPLEXITY)+"_iso_ground_truth.gtf"
	shell:
		"""
		mkdir -p {GROUND_TRUTH_BASE_PATH}/complexityFilter
		{GTFFILTER_BINARY} -li {LAST_COMPLEXITY} -ui {LAST_COMPLEXITY} {input.infile} {output[0]}
		"""


rule run_generateTPMFilteredListofTranscripts:
	input:
		abundance_file = {ABUNDANCE_FILE}
	output:
		GROUND_TRUTH_BASE_PATH+"/filteredTranscriptList.txt"
	shell:
		"""
		awk '$6>{DESIRED_TPM} {{print $1}}' {input.abundance_file} > {output[0]}
		"""


rule run_create_tpmFiltered_GroupedGenesFiles:
	input:
		GROUND_TRUTH_BASE_PATH+"/filteredTranscriptList.txt",
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_fullyCovered.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_fullyCov_"+DESIRED_TPM+"Tpm.gtf",
		GROUND_TRUTH_BASE_PATH+"/groupedGenes_fullyCov_"+DESIRED_TPM+"Tpm.fa"
	shell:
		"""
		{GTFFILTER_BINARY} -m whiteT -l {input[0]} {input[1]} {output[0]}
		{GFFREAD_BINARY} {output[0]} -g {GENOME_FASTA_FILE} -w {output[1]}
		"""

rule run_GtfFilter_filterByTpm:
	input:
		infile = GROUND_TRUTH_BASE_PATH+"/complexityFilter/{complexity}_iso_ground_truth.gtf",
		filterBy = GROUND_TRUTH_BASE_PATH+"/filteredTranscriptList.txt"
	output:
		temp(GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/{complexity}_tpm_ground_truth_cleaned.gtf"),
		GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/{complexity}_tpm_ground_truth.gtf"
	shell:
		"""
		mkdir -p {GROUND_TRUTH_BASE_PATH}"/tpmFilter"{DESIRED_TPM}
		grep transcript_id {input.infile} > {output[0]}
		{GTFFILTER_BINARY}  --list {input.filterBy} -m whiteT {output[0]} {output[1]}
		"""

rule run_GtfFilter_filterByCoverage:
	input:
		infile = GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/{complexity}_tpm_ground_truth.gtf",
		filterBy = GROUND_TRUTH_BASE_PATH+"/fullyCoveredTranscripts.txt"
	output:
		temp(GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/CoverageFilter/{complexity}_iso_ground_truth_cleaned.gtf"),
		GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/CoverageFilter/{complexity}_iso_ground_truth.gtf"
	shell:
		"""
		mkdir -p {GROUND_TRUTH_BASE_PATH}"/tpmFilter"{DESIRED_TPM}"/CoverageFilter"
		grep transcript_id {input.infile} > {output[0]}
		{GTFFILTER_BINARY}  --list {input.filterBy} -m whiteT {output[0]} {output[1]}
		"""

rule run_GffRead_createFastaFiles:
	input:
		GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/CoverageFilter/{complexity}_iso_ground_truth.gtf"
	output:
		GROUND_TRUTH_BASE_PATH+"/tpmFilter"+DESIRED_TPM+"/CoverageFilter/{complexity}_iso_ground_truth.fa"
	shell:
		"""
		{GFFREAD_BINARY} {input[0]} -g {GENOME_FASTA_FILE} -w {output[0]}
		"""
