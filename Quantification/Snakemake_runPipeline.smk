from collections import defaultdict
import random

include: "config.py"
include: "Snakefile_runKallisto.smk"
include: "Snakefile_analyse.smk"
include: "Snakefile_plot.smk"
include: "Snakefile_simulate.smk"
include: "Snakefile_prepareSimulation.smk"

rule all:
    input:
        # ## Plotting across betas all included
        SIM_DATA_BASE_PATH+"/geneComp/allIncluded_geneCompTpm.svg",
        SIM_DATA_BASE_PATH+"/geneComp/allIncluded_geneCompMARDS.svg"
