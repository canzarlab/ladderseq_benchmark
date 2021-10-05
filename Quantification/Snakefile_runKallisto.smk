## Runs kallisto original and modified
rule run_kallisto_original:
	input:
		expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/merged/R{pair}_sim.fq", pair = [1,2])
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original/abundance.h5",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original/abundance.tsv",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original/run_info.json"
	run:
		outDir = "{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/original/"
		shell("mkdir -p {0}".format( outDir ))
		shell("kallisto quant -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/original/ -i {KAL_IDX} {input}")



def determineBeta(wildcards):
	otherInputs = [SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band1/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band1/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band2/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band2/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band3/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band3/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band4/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band4/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band5/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band5/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band6/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band6/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band7/R1_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band7/R2_sim.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/R_1.fq",
	SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/R_2.fq"]

	beta = wildcards.bandProb

	if beta=="betas_perfect":
		returnList = [str(BAND_PROB_FILE_PATH+beta+".txt")]
	elif beta!="betas_perfect" and BETA_ESTIMATION=="FALSE":
		returnList = [str(BAND_PROB_FILE_PATH+beta+".txt")]
	else:
		returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+beta+"/estimatedBeta_high_text.txt")]
	returnList.extend(otherInputs)
	return returnList



rule run_kallisto_modified:
	input:
		determineBeta,
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/outputMod/abundance.h5",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/outputMod/abundance.tsv",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/outputMod/run_info.json"
	run:
		outDir = "{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/outputMod/"
		shell("mkdir -p {0}".format( outDir ))
		shell("{KALLISTO_LS_BINARY} quant -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/outputMod/ -i {KAL_IDX} -p {input}")
