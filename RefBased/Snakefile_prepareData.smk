
rule get_geuvadis:
    output:
        "../configData/realData/NA12716_7_{pair}.fastq.gz"
    shell:
        "mkdir -p ../configData/realData;"
        "cd ../configData/realData; "
        "wget -O NA12716_7_{wildcards.pair}.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188140/ERR188140_{wildcards.pair}.fastq.gz; "


rule get_genome:
    output:
        "../configData/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    shell:
        "cd ../configData; "
        "wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; "
        "gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"


rule get_anno:
    input:
        "../configData/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        "../configData/Homo_sapiens.GRCh38.91.gtf",
        "../configData/Homo_sapiens.GRCh38.91.fa"
    shell:
        "cd ../configData; "
        "wget ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz; "
        "gunzip Homo_sapiens.GRCh38.91.gtf.gz; "
        "gffread Homo_sapiens.GRCh38.91.gtf -g Homo_sapiens.GRCh38.dna.primary_assembly.fa -w Homo_sapiens.GRCh38.91.fa; "


rule runSTARGenomeGenerate:
    input:
        ANOTATION_FASTA_FILE
    output:
        STAR_GENOME_DIR
    shell:
        "mkdir -p {STAR_GENOME_DIR};"
        "{STAR_BINARY} --runThreadN 4 --limitGenomeGenerateRAM 140026839306 --runMode genomeGenerate --genomeDir {STAR_GENOME_DIR} --genomeFastaFiles {GENOME_FASTA_FILE}"



# Prepare the main annotation to run RSEM on it
rule rsem_prepare:
    input:
        ANOTATION_FASTA_FILE
    output:
        ANNOTATION_RSEM+'.grp'
    run:
        shell('mkdir -p {CONFIG_DATA}/RSEMRef')
        shell('{RSEM_BINARY_PATH}/rsem-prepare-reference {ANOTATION_FASTA_FILE} {ANNOTATION_RSEM}')


rule kallisto_index:
    input:
        ANOTATION_FASTA_FILE
    output:
        KAL_IDX
    threads: 1
    shell:
        KALLISTO_LS_BINARY + ' index '
        '-i ' + KAL_IDX + ' ' +
        ANOTATION_FASTA_FILE


rule genome_bwt2_index:
    input:
        GENOME_FASTA_FILE
    output:
        GENOME_BWT + ".1.bt2",
        GENOME_BWT + ".2.bt2",
        GENOME_BWT + ".3.bt2",
        GENOME_BWT + ".4.bt2",
        GENOME_BWT + ".rev.1.bt2",
        GENOME_BWT + ".rev.2.bt2"
    threads: 32
    run:
        shell('mkdir -p {CONFIG_DATA}/bowtie')
        shell('{BOWTIE2_BINARY_PATH}/bowtie2-build --seed 42 --offrate 1 {input} {GENOME_BWT}')


rule anno_bwt2_index:
    input:
        ANOTATION_FASTA_FILE
    output:
        ANNOTATION_BWT + ".1.bt2",
        ANNOTATION_BWT + ".2.bt2",
        ANNOTATION_BWT + ".3.bt2",
        ANNOTATION_BWT + ".4.bt2",
        ANNOTATION_BWT + ".rev.1.bt2",
        ANNOTATION_BWT + ".rev.2.bt2"
    threads: 32
    run:
        shell('mkdir -p {CONFIG_DATA}/bowtie')
        shell('{BOWTIE2_BINARY_PATH}/bowtie2-build --seed 42 --offrate 1 {input} {ANNOTATION_BWT}')



########################################################################
# Real data to act as seed for simulation
########################################################################

# align for RSEM on "real" data
rule rsem_align:
    input:
        REAL_DATA_BASE_PATH+REAL_DATA_FILE_NAME+"_1"+REAL_DATA_FILE_TYPE,
        REAL_DATA_BASE_PATH+REAL_DATA_FILE_NAME+"_2"+REAL_DATA_FILE_TYPE,
        ANNOTATION_BWT + ".1.bt2"
    output:
        REAL_DATA_BASE_PATH+REAL_DATA_FILE_NAME+"_rsem.bam"
    threads:
        N_THREADS
    message:
        'Running bowtie2 (thr = {threads}) on {input}'
    shell:
        'bowtie2 '
        '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 '
        '--no-mixed --no-discordant '
        '--seed 42 '
        '-k 200 '
        '-x ' + ANNOTATION_BWT + ' '
        '-p {threads} '
        '-1 {input[0]} '
        '-2 {input[1]} | '
        '{SAMTOOLS_BINARY} view -Sb - > {output}'


# Run RSEM on a "real" sample
rule rsem:
    input:
        REAL_DATA_BASE_PATH+REAL_DATA_FILE_NAME+"_rsem.bam",
        ANNOTATION_RSEM + ".grp"
    output:
        REAL_DATA_BASE_PATH+"rsem/out.isoforms.results"
    threads:
        N_THREADS
    run:
        shell("mkdir -p {REAL_DATA_BASE_PATH}/rsem")
        out_prefix = "{REAL_DATA_BASE_PATH}/rsem/out"
        cmd = " ".join(["rsem-calculate-expression", "--bam", "--paired-end",
            "--no-bam-output", "-p {threads}", "{input[0]}", ANNOTATION_RSEM, out_prefix])
        shell(cmd)



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
                "{BAND_PROB_FILE_PATH}{wildcards.bandProb}_csv.txt",
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
                "{BAND_PROB_FILE_PATH}{wildcards.bandProb}_csv.txt",
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
		#shell("mkdir -p {0}".format( outDir ))
		shell("{KALLISTO_LS_BINARY} pseudo --pseudobam -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/ -i {input}")


rule getUniqueMappings:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/pseudoalignments.bam"
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/uniqueMappings.txt"
	run:
		shell("{SAMTOOLS_BINARY} view {input} | awk 'NR%2 == 0 && $12==\"NH:i:1\" {{print $1,\"\t\",$3,\"\t\",$12}}' > {output[0]}")


rule estimateBetas:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/uniqueMappings.txt"
	threads:
		N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt"
	run:
		shell("Rscript EstimateBetas.R {input[0]} 0 {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/")


def determineBands(wildcards):
    bandCombination = wildcards.combBand
    pair = wildcards.pair
    if bandCombination=="comb12":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band1/R"+pair+"_sim.fq"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band2/R"+pair+"_sim.fq")]
    elif bandCombination=="comb23":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band2/R"+pair+"_sim.fq"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band3/R"+pair+"_sim.fq")]
    elif bandCombination=="comb34":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band3/R"+pair+"_sim.fq"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band4/R"+pair+"_sim.fq")]
    elif bandCombination=="comb45":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band4/R"+pair+"_sim.fq"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band5/R"+pair+"_sim.fq")]
    elif bandCombination=="comb56":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band5/R"+pair+"_sim.fq"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band6/R"+pair+"_sim.fq")]
    elif bandCombination=="comb67":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band6/R"+pair+"_sim.fq"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band7/R"+pair+"_sim.fq")]

    print(returnList)
    return returnList



rule combineReads:
    input:
        determineBands
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/Comb_R{pair}_sim.fq"
    shell:
        """
            mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}
            cat {input} > {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}/Comb_R{wildcards.pair}_sim.fq
        """


rule run_STAR_Band_INDIVIDUAL:
	input:
		genomeDir = STAR_GENOME_DIR,
		inputFilePair1 = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}/R1_sim.fq",
		inputFilePair2 = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{band}/R2_sim.fq"
	threads:
	        N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/{band}/Aligned.sortedByCoord.out.bam"
	shell:
		"""
		mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.band}
		{STAR_BINARY} --runThreadN 8 --genomeDir {input.genomeDir} --readFilesIn {input.inputFilePair1} {input.inputFilePair2} --outSAMstrandField intronMotif --outFileNamePrefix {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.band}/ --outSAMtype BAM SortedByCoordinate
		"""

rule run_STAR_Merged:
	input:
		genomeDir = STAR_GENOME_DIR,
		inputFilePair1 = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R1_sim.fq",
		inputFilePair2 = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/R2_sims.fq"
	threads:
	        N_THREADS
	output:
		"{SIM_DATA_BASE_PATH}/sim_{simNum}/{bandProb}/merged/Aligned.sortedByCoord.out.bam"
	shell:
		"""
		mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/merged/
		{STAR_BINARY} --runThreadN 8 --genomeDir {input.genomeDir} --readFilesIn {input.inputFilePair1} {input.inputFilePair2} --outSAMstrandField intronMotif --outFileNamePrefix {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/merged/ --outSAMtype BAM SortedByCoordinate
		"""

rule run_STAR_Band_COMBINED:
        input:
                genomeDir = STAR_GENOME_DIR,
                inputFilePair1 = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/Comb_R1_sim.fq",
                inputFilePair2 = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/Comb_R2_sim.fq"
        threads:
                N_THREADS
        output:
                SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/Comb.Aligned.sortedByCoord.out.bam"
        shell:
                """
                {STAR_BINARY} --runThreadN 8 --genomeDir {input.genomeDir} --readFilesIn {input.inputFilePair1} {input.inputFilePair2} --outSAMstrandField intronMotif --outFileNamePrefix {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}/Comb. --outSAMtype BAM SortedByCoordinate
                """
