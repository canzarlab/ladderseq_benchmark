########################################################################
# Generate simulations
########################################################################

# simulate from RSEM using errors
rule rsem_simulate:
    input:
        REAL_DATA_BASE_PATH+"rsem/out.isoforms.results"
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/R_1.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/R_2.fq"
    threads: 1
    run:
        shell("mkdir -p ../simulations")
        shell("mkdir -p {SIM_DATA_BASE_PATH}/")
        res = "{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/R"
        cmd = " ".join([
            "rsem-simulate-reads",
            ANNOTATION_RSEM,
            "{REAL_DATA_BASE_PATH}rsem/out.stat/out.model",
            "{REAL_DATA_BASE_PATH}rsem/out.isoforms.results",
            "0.0",
            "{N_READS}",
            "{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/R",
            "--seed", "{wildcards.simNum}"
            ])
        shell(cmd)
        #shell("gzip -f {0}*.fq".format( res ))


########################################################################
# Split reads to bands
########################################################################

rule split_reads_to_bands:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/R_1.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/R_2.fq"
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band1_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band1_R2_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band2_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band2_R2_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band3_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band3_R2_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band4_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band4_R2_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band5_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band5_R2_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band6_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band6_R2_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band7_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band7_R2_sim.fq",
    run:
        shell('mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}')
        if ({wildcards.bandProb} == {'betas_perfect'}):
            cmd = " ".join([
                "python split_reads_into_bands.py",
                ANOTATION_GTF_FILE,
                "{REAL_DATA_BASE_PATH}rsem/out.isoforms.results",
                "{input}",
                "{CONFIG_DATA}{BAND_PROB_FILE_PATH}{wildcards.bandProb}_csv.txt",
                "{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/",
                "RSEM",
                "LOW",
                "INITIAL"
                ])
        else:
            cmd = " ".join([
                "python split_reads_into_bands.py",
                ANOTATION_GTF_FILE,
                "{REAL_DATA_BASE_PATH}rsem/out.isoforms.results",
                "{input}",
                "{CONFIG_DATA}{BAND_PROB_FILE_PATH}{wildcards.bandProb}_csv.txt",
                "{SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/",
                "RSEM",
                "HIGH",
                "INITIAL"
                ])
        shell(cmd)

rule organiseReadsInFolders:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}_R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}_R2_sim.fq"
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}/R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}/R2_sim.fq"
    run:
        shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.band}")
        shell("mv {input[0]} {output[0]}")
        shell("mv {input[1]} {output[1]}")


rule mergeReads:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band1/R{pair}_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band2/R{pair}_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band3/R{pair}_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band4/R{pair}_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band5/R{pair}_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band6/R{pair}_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/band7/R{pair}_sim.fq",
    threads:
        N_THREADS
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R{pair}_sim.fq"
    run:
        shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/merged/")
        shell("python interleave_fqs.py {input} | sed '/^$/d' > {output}")




## Runs kallisto original and modified
rule run_kallisto_pseudoalignment:
    input:
        KAL_IDX,
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R1_sim.fq",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R2_sim.fq"
    threads:
        N_THREADS
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/pseudoalignments.bam"
    run:
        shell("{KALLISTO_LS_BINARY} pseudo --pseudobam -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/ -i {input}")


rule getUniqueMappings:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/pseudoalignments.bam"
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/uniqueMappings.txt"
	run:
		shell("samtools view {input} | awk 'NR%2 == 0 && $12==\"NH:i:1\" {{print $1,\"\t\",$3,\"\t\",$12}}' > {output[0]}")


rule estimateBetas:
    input:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/uniqueMappings.txt"
    threads:
        N_THREADS
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/geneId_transId.txt",
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt"
    run:
        shell("Rscript EstimateBetas.R {input[0]} 0 {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/")
