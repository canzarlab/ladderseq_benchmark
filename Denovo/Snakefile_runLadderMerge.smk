### Snakemake file for ###
### 1. Denovo assembly using trinity ###

### Author: Shounak ###
### 24.0.2019 ###




rule ladder_seq_trinity_individual:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}/R1" + READS_SFX + "." + READS_TYPE,
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}/R2" + READS_SFX + "." + READS_TYPE
    threads:
    	N_THREADS
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/{band}_trinity/Trinity.fasta"
    message: "Running Trinity on {wildcards.band} of {wildcards.simNum} of {wildcards.bandProb}..."
    shell:
        """
        mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.band}_trinity/
        {TRINITY_BINARY} --seqType {READS_TYPE} --left {input[0]} --right {input[1]} --CPU 16 --max_memory 20G --output {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.band}_trinity
        """

rule filter_trinity_individual:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}_trinity/Trinity.fasta",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt"
    threads:
    	N_THREADS
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/{band}/Trinity_filtered.fasta"
    message: "Filtering Trinity output on {wildcards.band} of {wildcards.simNum} of {wildcards.bandProb}..."
    shell:
        """
        python FilterLadderFasta.py {input[0]} {output[0]} length {wildcards.band} {FILTER_EXTENT} {FILTER_RES} {input[1]}
        """


rule concat_contigs:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/"+"Trinity_filtered.fasta", band=BANDS)
    threads:
        N_THREADS
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/concat/concatenated.fasta"
    shell:
        """
        mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/concat
        cat {input} | grep -v '^#' > {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/concat/concatenated.fasta
        """


rule run_kallistoMod_index:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/concat/concatenated.fasta"
	threads:
	        N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/ladderMerge/modKallistoIndex.idx"
	shell:
		"""
		mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/ladderMerge
		{KALLISTO_LS_BINARY} index --make-unique -i {output[0]} {input[0]}
		"""


rule run_kallisto_modified_denovo:
	input:
		kalIndex = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/ladderMerge/modKallistoIndex.idx",
		bandProbFile = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt",
		ladderReads = expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/R{pair}_sim.fq",band=BANDS,pair=[1,2]),
		mergedReads = expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/merged/R{pair}_sim.fq",pair=[1,2])
	threads:
	        N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/ladderMerge/abundance.tsv"
	shell:
		"""
		{KALLISTO_LS_BINARY} quant -t {threads} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/ladderMerge/ -i {input.kalIndex} -p {input.bandProbFile} {input.ladderReads} {input.mergedReads}
		"""


rule run_filterModKallistoGff:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/concat/concatenated.fasta",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/ladderMerge/abundance.tsv"
    threads:
        N_THREADS
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/ladderMerge/transcripts_tobeKept.txt",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/ladderMerge/ladderMergeFinal.fasta"
    shell:
    ## the fourth and 5th argument is filterByTpm and filterByCount
        """
        awk '$5 >= {MOD_KALLISTO_FILTER} {{print $1}}' {input[1]} > {output[0]}
        python FilterLadderFasta.py {input[0]} {output[1]} id {output[0]}
        """

# Run Trinity on full set of reads
rule original_trinity:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R1" + READS_SFX + "." + READS_TYPE,
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R2" + READS_SFX + "." + READS_TYPE
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/original_trinity/Trinity.fasta"
    message: "Running Trinity on full set..."
    shell:
        """
        mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/original_trinity
		Trinity --seqType {READS_TYPE} --left {input[0]} --right {input[1]} --CPU 16 --max_memory 20G --output {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/original_trinity/
		"""

# Run BLAT using the Ladder-Seq results
rule blat_ladder_seq:
    input:
        DENOVO_FASTA,
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/ladderMerge/ladderMergeFinal.fasta"
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/blat/ladder_seq.psl"
    shell:
        """
        mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/blat
        {BLAT_BINARY} {input} -t=dna -q=dna -minIdentity={BLAT_MIN_IDENT} {output}
        """


# Run BLAT using the full set results
rule blat_full:
    input:
        DENOVO_FASTA,
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original_trinity/Trinity.fasta"
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/blat/full.psl"
    shell:
        """
        mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/blat
        {BLAT_BINARY} {input} -t=dna -q=dna -minIdentity={BLAT_MIN_IDENT} {output}
        """


# Run blat for ladder seq per complexity
rule blat_ladder_seq_complexity:
	input:
		COMPLETE_GROUND_TRUTH_FILE_PATH+"/{complexity}_iso_ground_truth.fa",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/ladderMerge/ladderMergeFinal.fasta"
	threads:
		N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/blat/complexity{complexity}/ladder_seq.psl"
	run:
		shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/blat/complexity{wildcards.complexity}")
		shell("{BLAT_BINARY} {input[0]} {input[1]} -t=dna -q=dna -minIdentity={BLAT_MIN_IDENT} {output}")

# Run BLAT using the full set per complexity
rule blat_full_complexity:
	input:
		COMPLETE_GROUND_TRUTH_FILE_PATH + "/{complexity}_iso_ground_truth.fa",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original_trinity/Trinity.fasta"
	threads:
		N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/blat/complexity{complexity}/full.psl"
	shell:
		"{BLAT_BINARY} {input[0]} {input[1]} -t=dna -q=dna -minIdentity={BLAT_MIN_IDENT} {output}"



# Calculate precision and recall of both BLAT results using R script
rule blat_analysis_orig_byComplexity:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/ladderMerge/ladderMergeFinal.fasta",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original_trinity/Trinity.fasta",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/blat/ladder_seq.psl",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/blat/full.psl",
        DENOVO_FASTA
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/denovoResults.txt"
    shell:
        """
        TS=$(grep '>' {input[4]} | wc -l)
        S1=$(grep '>' {input[0]} | wc -l)
        S2=$(grep '>' {input[1]} | wc -l)
        Rscript trinityAnalysis.R {input[2]} {input[3]} ${{S1}} ${{S2}} ${{TS}} {output} {wildcards.bandProb}
        """

# Calculate precision and recall of complexity wise BLAT results using R script
rule blat_analysis_ladder_byComplexity:
	input:
		COMPLETE_GROUND_TRUTH_FILE_PATH + "/{complexity}_iso_ground_truth.fa",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/ladderMerge/ladderMergeFinal.fasta",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original_trinity/Trinity.fasta",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/blat/complexity{complexity}/ladder_seq.psl",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/blat/complexity{complexity}/full.psl"
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/blat/complexity{complexity}/results.txt"
	shell:
		"""
			TS=$(grep '>' {input[0]} | wc -l)
			S1=$(grep '>' {input[1]} | wc -l)
			S2=$(grep '>' {input[2]} | wc -l)
			Rscript trinityAnalysis.R {input[3]} {input[4]} ${{S1}} ${{S2}} ${{TS}} {output} {wildcards.bandProb}
		"""

# Assimilate results
rule blat_assimilateResults:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/blat/complexity{complexity}/results.txt", complexity=COMPLEXITIES),
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/denovoResults.txt"
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/CompRe_ladder.txt",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/CompRe_original.txt"
    shell:
        """
            bash assimilateResults.sh {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/blat/ {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/
        """


# Assimilate results across betas for ladder seq by complexity
rule blat_assimilateResults_ladderComplexities_acrossBetas:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{bandProb}/CompRe_ladder.txt", bandProb=BAND_PROB)
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/CompRe_ladder_acrossBetas.txt"
    shell:
        """
            cat {input} > {output}
        """

# Assimilate results across betas for conventional by complexity
rule blat_assimilateResults_originalComplexities_acrossBetas:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{bandProb}/CompRe_original.txt",  bandProb=["betas_realistic"])
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/CompRe_original_acrossBetas.txt"
    shell:
        """
            cat {input} > {output}
        """


# Assimilate results for both ladder and original by all complexities together
rule blat_assimilateResults_acrossBetas:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{bandProb}/denovoResults.txt", bandProb=BAND_PROB)
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/denovoResults_acrossBetas.txt"
    shell:
        """
            cat {input} > {output}
        """




# Plot results
rule blat_Plot_individialBetas:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/CompRe_ladder.txt",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/CompRe_original.txt",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/denovoResults.txt"
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/Complexity_Re_80.png",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/Complexity_Re_90.png",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/Complexity_Re_95.png",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/PrecisionComparison.png",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/RecallComparison.png"
    shell:
        """
            Rscript AnalyseResults.R {input[1]} {input[0]} {input[2]} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/ ./ {DATASET}
        """





# Plot results
rule blat_Plot_acrossBetas:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/CompRe_ladder_acrossBetas.txt",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/CompRe_original_acrossBetas.txt",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/denovoResults_acrossBetas.txt"
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/Complexity_Re_80_acrossBetas.svg",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/Complexity_Re_90_acrossBetas.svg",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/Complexity_Re_95_acrossBetas.svg",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/PrecisionComparison_acrossBetas.svg",
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/RecallComparison_acrossBetas.svg"
    shell:
        """
            Rscript AnalyseResults_denovoAssembly.R {input[1]} {input[0]} {input[2]} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/ ./ {DATASET}
        """
