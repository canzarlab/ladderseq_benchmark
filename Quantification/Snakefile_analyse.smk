
rule run_analysis:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/outputMod/abundance.tsv",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/original/abundance.tsv",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/geneId_transId.txt"
	threads:
	        N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/geneComp/allIncluded_MARDS_Mod_File.txt",
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/geneComp/allIncluded_MARDS_Orig_File.txt",
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/geneComp/allIncluded_tpm_Mod_File.txt",
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/geneComp/allIncluded_tpm_Orig_File.txt"

	run:
		shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/geneEffectiveComp")
		shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/geneComp")
		shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/difference")
		cmd = " ".join([
            "Rscript WriteComparisonFiles.R",
			"{input[0]}",
			"{input[1]}",
			"{REAL_DATA_BASE_PATH}rsem/out.isoforms.results",
			"{wildcards.bandProb}",
			"{input[2]}",
			"{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/",
			"RSEM"
            ])
		shell(cmd)




rule run_assimilate_MARDSMod_allIncluded:
        input:
                expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/{{bandProb}}/geneComp/allIncluded_MARDS_Mod_File.txt", simNum=SIMS_LIST),
        threads:
                N_THREADS
        output:
                SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_MARDS_Mod_File_Aggregated.txt",
        run:
                shell("mkdir -p {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp")
                shell("Rscript AggregateResults.R {NUM_SIM} geneComp {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp/ allIncluded_MARDS_Mod_File_Aggregated.txt {wildcards.bandProb} {input}")



rule run_assimilate_MARDSOrig_allIncluded:
        input:
                expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/{{bandProb}}/geneComp/allIncluded_MARDS_Orig_File.txt", simNum=SIMS_LIST),
        threads:
                N_THREADS
        output:
                SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_MARDS_Orig_File_Aggregated.txt",
        run:
                shell("mkdir -p {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp")
                shell("Rscript AggregateResults.R {NUM_SIM} geneComp {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp/ allIncluded_MARDS_Orig_File_Aggregated.txt {wildcards.bandProb} {input}")




rule run_assimilate_tpmMod_allIncluded:
        input:
                expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/{{bandProb}}/geneComp/allIncluded_tpm_Mod_File.txt", simNum=SIMS_LIST),
        threads:
                N_THREADS
        output:
                SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_tpm_Mod_File_Aggregated.txt",
        run:
                shell("mkdir -p {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp")
                shell("Rscript AggregateResults.R {NUM_SIM} geneComp {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp/ allIncluded_tpm_Mod_File_Aggregated.txt {wildcards.bandProb} {input}")



rule run_assimilate_tpmOrig_allIncluded:
        input:
                expand(SIM_DATA_BASE_PATH+"/sim_{simNum}/{{bandProb}}/geneComp/allIncluded_tpm_Orig_File.txt", simNum=SIMS_LIST),
        threads:
                N_THREADS
        output:
                SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_tpm_Orig_File_Aggregated.txt",
        run:
                shell("mkdir -p {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp")
                shell("Rscript AggregateResults.R {NUM_SIM} geneComp {SIM_DATA_BASE_PATH}/{wildcards.bandProb}/geneComp/ allIncluded_tpm_Orig_File_Aggregated.txt {wildcards.bandProb} {input}")











rule run_assimilate_MARDSMod_acrossBetas_allIncluded:
	input:
		expand(SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_MARDS_Mod_File_Aggregated.txt", bandProb=BAND_PROB),
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_MARDS_Mod_File_acrossBetas.txt",
	run:
		shell("mkdir -p {SIM_DATA_BASE_PATH}/geneComp")
		shell("cat {input} > {output}")



rule run_assimilate_MARDSOrig_acrossBetas_allIncluded:
	input:
		expand(SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_MARDS_Orig_File_Aggregated.txt", bandProb=["betas_realistic"]),
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_MARDS_Orig_File_acrossBetas.txt",
	run:
		shell("mkdir -p {SIM_DATA_BASE_PATH}/geneComp")
		shell("cat {input} > {output}")


rule run_assimilate_tpmMod_acrossBetas_allIncluded:
	input:
		expand(SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_tpm_Mod_File_Aggregated.txt", bandProb=BAND_PROB)
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_tpm_Mod_File_acrossBetas.txt",
	run:
		shell("mkdir -p {SIM_DATA_BASE_PATH}/geneComp")
		shell("cat {input} > {output}")



rule run_assimilate_tpmOrig_acrossBetas_allIncluded:
	input:
		expand(SIM_DATA_BASE_PATH+"/{bandProb}/geneComp/allIncluded_tpm_Orig_File_Aggregated.txt", bandProb=["betas_realistic"])
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_tpm_Orig_File_acrossBetas.txt",
	run:
		shell("mkdir -p {SIM_DATA_BASE_PATH}/geneComp")
		shell("cat {input} > {output}")
