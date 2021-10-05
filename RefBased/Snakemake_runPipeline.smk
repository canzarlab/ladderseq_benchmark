### Snakemake file for ###


### Author: Shounak ###
### 30.07.2018 ###

include: "config.py"
include: "Snakefile_prepareData.smk"
include: "Snakefile_runAssembly.smk"
include: "Snakefile_runLadderAssembly.smk"
include: "Snakefile_runEvaluation.smk"
include: "Snakefile_runGroundTruthCreation.smk"


rule all:
	input:
		## Plot Stringtie results acrossBetas
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/Stringtie_Pre_gffCompare_acrossBetas.svg", simNum=SIMS_LIST),
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/Stringtie_F1_gffCompare_acrossBetas.svg", simNum=SIMS_LIST),
		expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/Stringtie_Re_gffCompare_acrossBetas.svg", simNum=SIMS_LIST),
