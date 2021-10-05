### Snakemake file for Running {GFFCOMPARE_BINARY} comparing the stringtie outputs to the provided ground truth

### Author: Shounak ###
### 23.05.2018 ###


rule run_gffcompare_original:
	input:
		groundTruthGff = GROUND_TRUTH_PATH + "/{complexity}_iso_ground_truth.gtf",
		completeGroundTruthGff = COMPLETE_GROUND_TRUTH_FILE_PATH + "/{complexity}_iso_ground_truth.gtf",
		stringtieOutputOriginal = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/stringtie_out_original.gtf",
		splitAssemblyOutput = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/finalAssembly_{complexity}.gtf"
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/gffComp_{complexity}.stats",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/CompleteGround_{complexity}.stats"
	shell:
		"""
		{GFFCOMPARE_BINARY} -Q -r {input.groundTruthGff} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/originalStringtie/gffComp_{wildcards.complexity} {input.stringtieOutputOriginal}
		{GFFCOMPARE_BINARY} -Q -r {input.completeGroundTruthGff} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/originalStringtie/CompleteGround_{wildcards.complexity} {input.splitAssemblyOutput}
		"""

rule run_gffcompare_singleExonsRemoved:
	input:
		groundTruthGff=GROUND_TRUTH_PATH + "/{complexity}_iso_ground_truth.gtf",
		completeGroundTruthGff = COMPLETE_GROUND_TRUTH_FILE_PATH + "/{complexity}_iso_ground_truth.gtf",
		singleExonRemovedAssembly=SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedAssembly.gtf",
		splitAssemblyOutput = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/finalAssembly_{complexity}.gtf"
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/gffComp_{complexity}.stats",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/CompleteGround_{complexity}.stats"
	shell:
		"""
		{GFFCOMPARE_BINARY} -Q -r {input.groundTruthGff} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/gffComp_{wildcards.complexity} {input.singleExonRemovedAssembly}
		{GFFCOMPARE_BINARY} -Q -r {input.completeGroundTruthGff} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/CompleteGround_{wildcards.complexity} {input.splitAssemblyOutput}
		"""



## Assimilating gffcompare results complexity divided
rule assimilate_gffcompare_original:
	input:
		expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/originalStringtie/gffComp_{complexity}.stats",complexity=COMPLEXITIES),
		expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/singleExonsRemovedAssembly/gffComp_{complexity}.stats",complexity=COMPLEXITIES)
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffOriginalResults.txt"
	shell:
		"""
		bash assimilateResults.sh gff {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/originalStringtie/ {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/gffOriginalResults {wildcards.bandProb}
		"""
rule assimilate_singleExonsRemovedAssembly:
	input:
		expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/originalStringtie/gffComp_{complexity}.stats",complexity=COMPLEXITIES),
		expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/singleExonsRemovedAssembly/gffComp_{complexity}.stats",complexity=COMPLEXITIES)
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffSingleExonsRemovedAssemblyResults.txt"
	shell:
		"""
		bash assimilateResults.sh gff {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/ {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/gffSingleExonsRemovedAssemblyResults {wildcards.bandProb}
		"""

# Assimilate results across betas for original Stringtie complexity divided
rule assimilateResults_gffcompare_original_acrossBetas:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{bandProb}/gffOriginalResults.txt",  bandProb=["betas_realistic"])
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/gffOriginalResults_acrossBetas.txt"
    shell:
        """
            cat {input} > {output}
        """
# Assimilate results across betas for Stringtie-ls complexity divided
rule assimilateResults_gffcompare_ladder_acrossBetas:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{bandProb}/gffSingleExonsRemovedAssemblyResults.txt", bandProb=BAND_PROB)
    output:
        "{SIM_DATA_BASE_PATH}/sim_{simNum}/gffSingleExonsRemovedAssemblyResults_acrossBetas.txt"
    shell:
        """
            cat {input} > {output}
        """



rule plot_gffcompareResults_acrossBetas:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/gffOriginalResults_acrossBetas.txt",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/gffSingleExonsRemovedAssemblyResults_acrossBetas.txt"
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/Stringtie_Pre_gffCompare_acrossBetas.svg",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/Stringtie_F1_gffCompare_acrossBetas.svg",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/Stringtie_Re_gffCompare_acrossBetas.svg"
	shell:
		"""
		Rscript AnalyseResults_refBased_acrossBetas.R {input} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/ {DATASET} gffCompare ShounaksPipeline
		"""
