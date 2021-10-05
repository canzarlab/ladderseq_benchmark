### Snakemake file for ###
### 1. Running stringtie on alignment files from STAR ###

### Author: Shounak ###
### 11.07.2018 ###


rule run_stringtie_individual:
  input:
    alignmentFile = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{individualBand}/Aligned.sortedByCoord.out.bam",
    bandProb = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt"
  threads:
    N_THREADS
  output:
    SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{individualBand}/stringtie_out_individual.gtf",
    SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{individualBand}/stringtie_out_individual_filteredOut.gtf"
  run:
    shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.individualBand}")
    if STRINGTIE_LOW_THRESHOLD=="TRUE":
      shell("{STRINGTIE2_BINARY} {input.alignmentFile} -c {INDIVIDUAL_STRINGTIE_THRESHOLD} -o {output[0]}")
    else:
      shell("{STRINGTIE2_BINARY} {input.alignmentFile} -o {output[0]}")
    shell("sed -i '/^#/ d' {output[0]}")
    shell("sed -i 's/STRG/STRG{wildcards.individualBand}/g' {output[0]}")
    if FILTER_TYPE_INDIVIDUAL!="NONE":
      shell("Rscript GffLadderFilter.R {output[0]} {LADDER_RESOLUTION} {FILTER_TYPE_INDIVIDUAL} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.individualBand}/filteredAssembly_stringtie {input.bandProb} {LADDER_HIGH_FILTER_THRESHOLD} INDIVIDUAL {POLY_A} {wildcards.individualBand}")
      shell("mv {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.individualBand}/filteredAssembly_stringtie.gtf {output[0]}")
      shell("mv {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.individualBand}/filteredAssembly_stringtie_filteredOut.gtf {output[1]}")


rule run_stringtie_combined:
  input:
    alignmentFile = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/Comb.Aligned.sortedByCoord.out.bam",
    bandProb = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt"
  threads:
    N_THREADS
  output:
    SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/stringtie_out_combined.gtf",
    SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/stringtie_out_combined_filteredOut.gtf"
  run:
    shell("mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}")
    shell("{STRINGTIE2_BINARY} {input.alignmentFile} -o {output[0]}")
    shell("sed -i '/^#/ d' {output[0]}")
    shell("sed -i 's/STRG/STRG{wildcards.combBand}/g' {output[0]}")
    if FILTER_TYPE_COMBINED!="NONE":
      shell("Rscript GffLadderFilter.R {output[0]} {LADDER_RESOLUTION} {FILTER_TYPE_COMBINED} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}/filteredAssembly_stringtie {input.bandProb} {LADDER_HIGH_FILTER_THRESHOLD} COMBINED {POLY_A} {wildcards.combBand}")
      shell("mv {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}/filteredAssembly_stringtie.gtf {output[0]}")
      shell("mv {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}/filteredAssembly_stringtie_filteredOut.gtf {output[1]}")


rule run_stringtie_original:
  input:
    alignmentFile = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/merged/Aligned.sortedByCoord.out.bam"
  threads:
    N_THREADS
  output:
    SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/stringtie_out_original.gtf"
  shell:
    """
      mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/originalStringtie
      {STRINGTIE2_BINARY} {input.alignmentFile} -o {output}
      sed -i '/^#/ d' {output}
    """

rule run_splitPredictions_original:
	input:
		completeGroundTruthGff = COMPLETE_GROUND_TRUTH_FILE_PATH + "/{complexity}_iso_ground_truth.gtf",
		assemblyOutput = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/stringtie_out_original.gtf"
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/originalStringtie/finalAssembly_{complexity}.gtf"
	shell:
		"""
			bash PredictionSplitterCommands.sh {wildcards.complexity} {COMPLETE_GROUND_TRUTH_FILE_PATH} {input.assemblyOutput} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/originalStringtie/ {EXONREFINE_BINARY} {GROUPGENES_BINARY} {GTFFILTER_BINARY}
		"""


def determineBands(wildcards):
    bands = wildcards.combBand
    print(bands)
    if bands=="comb12":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band1/stringtie_out_individual.gtf"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band2/stringtie_out_individual.gtf")]
    elif bands=="comb23":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band2/stringtie_out_individual.gtf"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band3/stringtie_out_individual.gtf")]
    elif bands=="comb34":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band3/stringtie_out_individual.gtf"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band4/stringtie_out_individual.gtf")]
    elif bands=="comb45":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band4/stringtie_out_individual.gtf"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band5/stringtie_out_individual.gtf")]
    elif bands=="comb56":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band5/stringtie_out_individual.gtf"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band6/stringtie_out_individual.gtf")]
    elif bands=="comb67":
        returnList = [str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band6/stringtie_out_individual.gtf"),str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/band7/stringtie_out_individual.gtf")]

    returnList.append(str(SIM_DATA_BASE_PATH+"/sim_"+wildcards.simNum+"/"+wildcards.bandProb+"/"+wildcards.combBand+"/stringtie_out_combined.gtf"))

    return returnList


rule combineAssemblies:
    input:
        determineBands,
    output:
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/{combBand}/pair.combined.gtf"
    shell:
        """
            mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}
            {GFFCOMPARE_BINARY} -e 0 -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/{wildcards.combBand}/pair -C {input}
        """


rule run_gffCompare_merge:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/pair.combined.gtf", band=COMBINED_BAND) if PAIRWISE=='TRUE' else expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/cidaneOutput_individual.gtf", band=INDIVIDUAL_BAND)
    output:
        temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffMerged.combined.gtf"),
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffMerged/gffMergedTranscripts.gtf"
    shell:
        """
        mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/gffMerged
        {GFFCOMPARE_BINARY} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/gffMerged {input}
        grep -v '^#' {output[0]} > {output[1]}
        """


rule run_gffCompare_mergeRejectedTranscriptsToOneFile:
    input:
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/stringtie_out_combined_filteredOut.gtf", band=COMBINED_BAND),
        expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/stringtie_out_individual_filteredOut.gtf", band=INDIVIDUAL_BAND)
    output:
        temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/rejectedTranscripts.combined.gtf"),
        SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/rejectedTranscripts_combined.gtf",

    shell:
        """
        {GFFCOMPARE_BINARY} -D -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/rejectedTranscripts {input}
        grep -v '^#' {output[0]} > {output[1]}
        """
