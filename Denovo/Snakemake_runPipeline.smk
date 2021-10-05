### Snakemake file for ###


### Author: Shounak ###
### 30.07.2018 ###

include: "config.py"
include: "Snakefile_prepareData.smk"
include: "Snakefile_runLadderMerge.smk"
include: "Snakefile_runGroundTruthCreation.smk"


rule all:
	input:
		## Plot Stringtie results acrossBetas
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/Complexity_Re_80_acrossBetas.svg", simNum=SIMS_LIST),
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/Complexity_Re_90_acrossBetas.svg", simNum=SIMS_LIST),
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/Complexity_Re_95_acrossBetas.svg", simNum=SIMS_LIST),
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/PrecisionComparison_acrossBetas.svg", simNum=SIMS_LIST),
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/RecallComparison_acrossBetas.svg", simNum=SIMS_LIST)
