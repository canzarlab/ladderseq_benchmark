
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



# Prepare the main annotation to run RSEM on it
rule rsem_prepare:
    input:
        ANOTATION_FASTA_FILE
    output:
        ANNOTATION_RSEM+'.grp'
    run:
        shell('mkdir -p {CONFIG_DATA}/RSEMRef')
        shell('rsem-prepare-reference {ANOTATION_FASTA_FILE} {ANNOTATION_RSEM}')


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
        shell('bowtie2-build --seed 42 --offrate 1 {input} {GENOME_BWT}')


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
        shell('bowtie2-build --seed 42 --offrate 1 {input} {ANNOTATION_BWT}')



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
        'samtools view -Sb - > {output}'


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
