### Snakemake file for ###
### 1. Single e ###

### Author: Shounak ###
### 23.07.2018 ###

rule run_stringtieMerge_gffMerge:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffMerged/gffMergedTranscripts.gtf"
	output:
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/concatStringtie.gtf"),
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/concatStringtieMerged.gtf"
	shell:
		"""
		mkdir -p {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly
		{STRINGTIE2_BINARY} --merge -o {output[0]} {input}
		grep -v '^#' {output[0]} > {output[1]}
		"""


rule run_gffCompare_compareGffandStringtieMerge:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffMerged/gffMergedTranscripts.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/concatStringtieMerged.gtf"
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/compareGff_str_merge.concatStringtieMerged.gtf.refmap",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/compareGff_str_merge.concatStringtieMerged.gtf.tmap"
	shell:
		"""
		{GFFCOMPARE_BINARY} -r {input[0]} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/compareGff_str_merge {input[1]}
		"""


rule run_generate_sets:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffMerged/gffMergedTranscripts.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/compareGff_str_merge.concatStringtieMerged.gtf.tmap"
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setG.txt", ## single exon transcripts in gffMerged assembly
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setS.txt", ## single exon transcripts in strMerge which are equal to single exon transcripts in gffMerged
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setU.txt", ## set of unknown single exons in stringtie merge
	shell:
		"""
		sed 's/\s/\t/g' {input[0]} | awk '$3 == "exon"'| datamash -s  -g 10 count 3 | awk '$2 == 1 {{print $1}}' | sed 's/\"//g' | sed 's/\;//g' > {output[0]}
		awk '($3 == "=" && $6 == 1) {{print $2}}' {input[1]} | sort > {output[1]}
		awk '($3 == "u" && $6 == 1) {{print $5}}' {input[1]} | sort > {output[2]}
		"""

rule run_setOperations_SetDifference:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setG.txt",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setS.txt",
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setA.txt",
	shell:
		"""
		comm -23 {input[0]} {input[1]} > {output[0]}
		"""


rule run_filterSingleExons_gffFilter:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/gffMerged/gffMergedTranscripts.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setA.txt",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/concatStringtieMerged.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/setU.txt"
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_withSubfrags.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_C.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_X.gtf"
	shell:
		"""
			{GTFFILTER_BINARY} --list {input[1]} -m blackT {input[0]} {output[0]}
			{GFFCOMPARE_BINARY} -e 0 -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/discardAllSubfrags_C -C {output[0]}
			mv {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/discardAllSubfrags_C.combined.gtf {output[1]}
			{GFFCOMPARE_BINARY} -e 0 -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/discardAllSubfrags_X -X {output[0]}
			mv {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/discardAllSubfrags_X.combined.gtf {output[2]}
			awk '{{for (i=1; i<=NF; ++i) {{ if ($i ~ "transcript_id") print $(i+1)}} }}' {output[1]} | sed 's/\"//g' | sed 's/\;//g' | sort | uniq > {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/transcriptList_C.txt
			awk '{{for (i=1; i<=NF; ++i) {{ if ($i ~ "transcript_id") print $(i+1)}} }}' {output[2]} | sed 's/\"//g' | sed 's/\;//g' | sort | uniq > {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/transcriptList_X.txt
			comm -23 {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/transcriptList_C.txt {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/transcriptList_X.txt > {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/presentInCNotX.txt
		"""

rule run_subExonFilter:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_C.gtf"
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.gtf"
	shell:
		"""
		bash FindUniqSubExons.sh {input[0]} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/ {SUBEXON_LENGTH_THRESHOLD}
		comm -12 {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/transcripts_uniqueSubexonLength_threshold.txt {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/presentInCNotX.txt > {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/uniqueNotInX.txt
		{GTFFILTER_BINARY} --list {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/uniqueNotInX.txt -m blackT {input[0]} {output[0]}
		"""



# rule make_listOfTranscriptsToBeKept:
# 	input:
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.gtf",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/rejectedTranscripts_combined.gtf"
# 	threads:
# 		N_THREADS
# 	output:
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_rejectedTrans_tagged.gtf",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/TranscriptsToBePutBack.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/TranscriptsToBePutBack_sorted.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_concat.combined.gtf",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/transId_oId_mapping.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_transcripts_uniqueSubexonLength_threshold.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_transcripts_uniqueSubexonLength_threshold_sorted.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_sortedRejected_newTransId.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_transcripts_finalList_newTransId.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_transcripts_finalList_exons.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_transcripts_finalList_together.txt",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/TranscriptsToBePutBack.gtf",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.restored.combined.gtf",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_tagged.gtf",
# 		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_restored.gtf",
# 	shell:
# 		"""
# 			sed 's/TCONS_/TCONS_rejected_/g' {input[1]} > {output[0]}
# 			python MergeRejectedTranscripts.py {input[0]} {output[0]} {output[1]}
# 			sort {output[1]} > {output[2]}
# 			{GFFCOMPARE_BINARY} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/restored_concat {input[0]} {output[0]}
# 			grep rejected {output[3]} | awk '{{print $10,$14}}' | sed 's/"//g' | sed 's/;//g' > {output[4]}
# 			bash FindUniqSubExons.sh {output[3]} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/restored_ 1000000000
# 			sort {output[5]} > {output[6]}
# 			awk '{{print $1}}' {output[4]} | sort > {output[7]}
# 			grep -w -f {output[7]} {output[6]} > {output[8]}
# 			grep -w -f {output[8]} {output[4]} | awk '{{print $2}}' > {output[9]}
# 			cat {output[2]} {output[9]} | sort | uniq > {output[10]}
# 			{GTFFILTER_BINARY} --list {output[2]} -m whiteT {output[0]} {output[11]}
# 			sed -i 's/TCONS_/TCONS_restored/g' {output[11]}
# 			{GFFCOMPARE_BINARY} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.restored {input[0]} {output[11]}
# 			Rscript TagRestoredTranscripts.R {output[12]} {output[13]}
# 			grep -v '^#' {output[13]} > {output[14]}
# 		"""

rule make_listOfTranscriptsToBeKept:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.gtf",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/rejectedTranscripts_combined.gtf"
	threads:
		N_THREADS
	output:
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restored_rejectedTrans_tagged.gtf"),
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/TranscriptsToBePutBack.txt"),
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/TranscriptsToBePutBack_sorted.txt"),
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/TranscriptsToBePutBack.gtf"),
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.restored.combined.gtf"),
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_tagged.gtf"),
		temp(SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_restored.gtf")
	shell:
		"""
			sed 's/TCONS_/TCONS_rejected_/g' {input[1]} > {output[0]}
			python MergeRejectedTranscripts.py {input[0]} {output[0]} {output[1]}
			sort {output[1]} > {output[2]}
			{GTFFILTER_BINARY} --list {output[2]} -m whiteT {output[0]} {output[3]}
			sed -i 's/TCONS_/TCONS_restored/g' {output[3]}
			{GFFCOMPARE_BINARY} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.restored {input[0]} {output[3]}
			Rscript TagRestoredTranscripts.R {output[4]} {output[5]}
			grep -v '^#' {output[5]} > {output[6]}
		"""



rule run_gffRead_SingleExon:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_restored.gtf" if RESTORE == 'TRUE' else SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.gtf",
		GENOME_FASTA_FILE
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/firstRoundAssembly_singleExons.fa"
	shell:
		"""
		{GFFREAD_BINARY}  {input[0]} -g {input[1]} -w {output[0]}
		"""



rule run_kallistoMod_SingleExon_index:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/firstRoundAssembly_singleExons.fa"
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/modKallistoIndex_singleExons.idx"
	shell:
		"""
		{KALLISTO_LS_BINARY} index -i {output[0]} {input[0]}
		"""


rule run_kallisto_SingleExon_modified:
	input:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/modKallistoIndex_singleExons.idx",
		bandProbFile = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/estimatedBeta_high_text.txt",
		ladderReads = expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/{band}/R{pair}_sim.fq",band=INDIVIDUAL_BAND,pair=[1,2]),
		mergedReads = expand(SIM_DATA_BASE_PATH+"/sim_{{simNum}}/{{bandProb}}/merged/R{pair}_sim.fq",pair=[1,2])
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/abundance.tsv"
	shell:
		"""
		{KALLISTO_LS_BINARY} quant -t {threads} -o {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/ -i {input[0]} -p {input.bandProbFile} {input.ladderReads} {input.mergedReads}
		"""






rule run_filterModKallistoGff_SingleExon:
	input:
		unfilteredGff = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly_restored.gtf" if RESTORE == 'TRUE' else SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedInitialAssembly.gtf",
		modKallistoQuantResult = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/abundance.tsv"
	threads:
	        N_THREADS
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/restoredTranscripts_tobeKept.txt",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/initialAssemblyTranscripts_tobeKept.txt",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/transcripts_tobeKept.txt",
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedAssembly.gtf"
	shell:
		"""
		grep "TCONS_restored" {input[1]} | awk '$5 >= {RESTORED_TPM_FILTER} {{print $1}}' > {output[0]}
		grep -v "TCONS_restored" {input[1]} | awk '$5 >= {MOD_KALLISTO_FILTER} {{print $1}}' > {output[1]}
		cat {output[0]} {output[1]} > {output[2]}
		{GTFFILTER_BINARY} --list {output[2]} -m whiteT {input[0]} {output[3]}
		"""



rule run_splitPredictions_singleExonsRemoved:
	input:
		completeGroundTruthGff = COMPLETE_GROUND_TRUTH_FILE_PATH + "/{complexity}_iso_ground_truth.gtf",
		assemblyOutput = SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/singleExonsRemovedAssembly.gtf"
	output:
		SIM_DATA_BASE_PATH+"/sim_{simNum}/{bandProb}/singleExonsRemovedAssembly/finalAssembly_{complexity}.gtf"
	shell:
		"""
			bash PredictionSplitterCommands.sh {wildcards.complexity} {COMPLETE_GROUND_TRUTH_FILE_PATH} {input.assemblyOutput} {SIM_DATA_BASE_PATH}/sim_{wildcards.simNum}/{wildcards.bandProb}/singleExonsRemovedAssembly/ {EXONREFINE_BINARY} {GROUPGENES_BINARY} {GTFFILTER_BINARY}
		"""
