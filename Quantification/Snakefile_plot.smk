## Runs kallisto original and modified


rule run_plot_acrossBetas_allIncluded:
	input:
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_MARDS_Mod_File_acrossBetas.txt",
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_MARDS_Orig_File_acrossBetas.txt",
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_tpm_Mod_File_acrossBetas.txt",
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_tpm_Orig_File_acrossBetas.txt"
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_geneCompTpm.svg",
		SIM_DATA_BASE_PATH+"/geneComp/allIncluded_geneCompMARDS.svg"
	run:
		shell("Rscript PlotQuantResults.R {SIM_DATA_BASE_PATH}/ {N_READS}Reads")
