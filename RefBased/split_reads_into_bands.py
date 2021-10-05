
# coding: utf-8

# In[1]:


## Input :
## 1. gtf file name
## 2. isoform results file (from RSEM simulation)
## 3. simulated reads files (from RSEM or POLYESTER)
## 4. band probablility file (LOW or HIGH resolution)
## 5. output directory
## 6. simulation type
## 7. resolution
## 8. type of ladder used


# In[20]:


from Bio import SeqIO
import HTSeq
import sys
import random
import numpy as np
import csv


# In[21]:


# Gtf annotation file to compute transcript length
#gtf_filename = '/algbio1/shounak/genomes/human/release_91/Homo_sapiens.GRCh38.91.gtf'
gtf_filename = sys.argv[1]


# This file has the mapping transcript idx and id
# For RSEM data
#isofrom_results_file = '/algbio1/sophia/simulated_data/DTU_paired_cutoff_GRCh38.91/fullmodel/simulation_files/SRR493366_sample11.txt'
# For Dereks data
#isofrom_results_file = '/algbio1/shounak/raw/simulated/human/derekWrongSeparation/DerekImperfect/ground_truth.csv'
isofrom_results_file = sys.argv[2]

# paired-end input files
#in_read_files = ['/algbio1/shounak/raw/simulated/sophiaSim/SRR493366_sample11/RSEM/merged_R1_sim.fq','/algbio1/shounak/raw/simulated/sophiaSim/SRR493366_sample11/RSEM/merged_R2_sim.fq']
#in_read_files = ['/algbio1/derek/SnakeMake/data/old_reads/R1.fasta','/algbio1/derek/SnakeMake/data/old_reads/R2.fasta']
in_read_files = [sys.argv[3],sys.argv[4]]

# Band Probability File
#bandProb = np.loadtxt(open("/algbio1/shounak/raw/simulated/sophiaSim/SRR493366_sample11/RSEM/highRes/bandProb_cko2_high_csv.txt", "rb"), delimiter=",")
bandProb = np.loadtxt(open(sys.argv[5], "rb"), delimiter=",")
#print bandProb


# Output directory
#out_dir = '/algbio1/shounak/raw/simulated/sophiaSim/SRR493366_sample11/RSEM/highRes/'
out_dir = sys.argv[6]

# Simulation Type (RSEM or POLYESTER)
#simulationType = "RSEM"
simulationType = sys.argv[7]

# Resolution of band probabilities (HIGH or LOW)
#resolution = "HIGH"
resolution = sys.argv[8]

# Type of Ladder used (INITIAL or PROMEGA)
#ladder = "INITIAL"
ladder = sys.argv[9]


# In[22]:


def get_transcript_length_dict(gtf_filename):
    gff_file = HTSeq.GFF_Reader(gtf_filename, end_included=True)
    #1. Go through the GFF file and get all exons for each transcript:
    transcript_exon_dict = {}
    for feature in gff_file:
        if feature.type == "exon":
            transcript_id = feature.attr['transcript_id']
            if transcript_id not in transcript_exon_dict:
                transcript_exon_dict[transcript_id] = list()
            transcript_exon_dict[transcript_id].append(feature)

    #2. Go through the transcripts and calculate the length of each
    # by adding up the lengths of its exons
    transcript_len_dict = {}
    for transcript_id in transcript_exon_dict.keys():
        transcript_len_dict[transcript_id] = 200 # PolyA tail
        for exon in transcript_exon_dict[transcript_id]:
            transcript_len_dict[transcript_id] += exon.iv.length + 1

    return transcript_len_dict


# In[23]:


def get_transcript_idx_id_dict(isofrom_results_file):
    # map between transcript index and its id
    transcript_idx_id_dict = {}
    idx = 1
    header = True
    with open(isofrom_results_file, "r") as in_lines:
        for line in in_lines:
            if header:
                header = False
                continue
            transcript_id = line.split('\t')[0]
            transcript_idx_id_dict[idx] = transcript_id
            idx += 1
    return transcript_idx_id_dict


# In[24]:


def get_band_ideal(length):

    #print(ladder,resolution)
    if(ladder=='INITIAL'):
        if(resolution=="LOW"):
            ## Initial ladder low resolution cuts
            band_list = [0,1000,1500,2000,3000,4000,6000]

            if length <= band_list[1]:
                val = 1
            elif length > band_list[1] and length <= band_list[2]:
                val = 2
            elif length > band_list[2] and length <= band_list[3]:
                val = 3
            elif length > band_list[3] and length <= band_list[4]:
                val = 4
            elif length > band_list[4] and length <= band_list[5]:
                val = 5
            elif length > band_list[5] and length <= band_list[6]:
                val = 6
            else:
                val = 7
        elif(resolution=="HIGH"):
            ## Initial ladder High resolution cuts
            #print(length)
            #band_list = (0,800,900,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4500,5000,5500,6000,6500,7000)
            band_list = [0,910,1083,1268,1440,1590,1757,1920,2070,2273,2462,2701,2932,3230,3554,3893,4525,5324,6378]
            if length <= band_list[1]:
                val = 1
            elif length > band_list[1] and length <= band_list[2]:
                val = 2
            elif length > band_list[2] and length <= band_list[3]:
                val = 3
            elif length > band_list[3] and length <= band_list[4]:
                val = 4
            elif length > band_list[4] and length <= band_list[5]:
                val = 5
            elif length > band_list[5] and length <= band_list[6]:
                val = 6
            elif length > band_list[6] and length <= band_list[7]:
                val = 7
            elif length > band_list[7] and length <= band_list[8]:
                val = 8
            elif length > band_list[8] and length <= band_list[9]:
                val = 9
            elif length > band_list[9] and length <= band_list[10]:
                val = 10
            elif length > band_list[10] and length <= band_list[11]:
                val = 11
            elif length > band_list[11] and length <= band_list[12]:
                val = 12
            elif length > band_list[12] and length <= band_list[13]:
                val = 13
            elif length > band_list[13] and length <= band_list[14]:
                val = 14
            elif length > band_list[14] and length <= band_list[15]:
                val = 15
            elif length > band_list[15] and length <= band_list[16]:
                val = 16
            elif length > band_list[16] and length <= band_list[17]:
                val = 17
            elif length > band_list[17] and length <= band_list[18]:
                val = 18
            else:
                val = 19
    elif(ladder=='PROMEGA'):
        ## Promega ladder (low resolution) cuts
        band_list = [0,281,623,955,1383,1908,2604,3638,4981,6583]

        if length <= band_list[1]:
            val = 1
        elif length > band_list[1] and length <= band_list[2]:
            val = 2
        elif length > band_list[2] and length <= band_list[3]:
            val = 3
        elif length > band_list[3] and length <= band_list[4]:
            val = 4
        elif length > band_list[4] and length <= band_list[5]:
            val = 5
        elif length > band_list[5] and length <= band_list[6]:
            val = 6
        elif length > band_list[6] and length <= band_list[7]:
            val = 7
        elif length > band_list[5] and length <= band_list[8]:
            val = 8
        elif length > band_list[5] and length <= band_list[9]:
            val = 9
        else:
            val = 10

    return val


# In[25]:


def get_band(randomNum,band_val_ideal):
    prob = 0;
    index = 0;
    for i in bandProb[band_val_ideal-1]:
        prob +=i
        if(randomNum<=prob):
            break
        else:
            index += 1
            continue
    return index+1


# In[35]:


transcript_len_dict = get_transcript_length_dict(gtf_filename)

if(simulationType=="RSEM"):
    transcript_idx_id_dict = get_transcript_idx_id_dict(isofrom_results_file)
    fq1 = SeqIO.parse(in_read_files[0], "fastq")
    fq2 = SeqIO.parse(in_read_files[1], "fastq")

if(simulationType=="POLYESTER"):
    fq1 = SeqIO.parse(in_read_files[0], "fasta")
    fq2 = SeqIO.parse(in_read_files[1], "fasta")

d = {}
# Range (1,8) for initial band list
# Range (1,11) for Promega band list
for i in range(1,8):
    d["band_" + str(i)+'_1'] = []
    d["band_" + str(i)+'_2'] = []

count = {}
numReads = 0
cc = 0
for rec1 in fq1:
    cc = cc + 1
    rec2 = next(fq2)

    # Extracting the transcriptId to which a particular read maps to
    if(simulationType=="POLYESTER"):
        transcript_id = (rec1.id.split('/')[1])
        transcript_id = (transcript_id.split('.')[0])
        if(len(transcript_id)!=15):
            print("Problem")
            continue
    elif(simulationType=="RSEM"):
        idx = int(rec1.id.split('_')[2])
        if idx == 0:
            ## These reads are generated from background noise
            ## assigning them to each of the band uniformly
            transcript_id = "NOISE"
        else:
            transcript_id = transcript_idx_id_dict[idx]

    if (transcript_id!="NOISE"):
        if(transcript_id not in transcript_len_dict):
            print("Transcript not present: "+transcript_id)
            continue

        transcript_length = transcript_len_dict[transcript_id]
        band_val_ideal = get_band_ideal(transcript_length)

        #### Assigning the read to a particular band based on the prior probabilities
        ## Get Random number between 0 and 1
        while True:
            randomNum = random.uniform(0, 1)
            band_val = get_band(randomNum,band_val_ideal)
            #print(band_val)
            if(band_val<=7 or band_val>=1):
                break

        if(band_val>7 or band_val<1):
            print("Problem: "+str(band_val))

    elif (transcript_id=="NOISE"):
        band_val = random.randint(1,7)
        #print(band_val)



    ## add the band id at the header of each read
    rec1.id = 'Band'+str(band_val)+':'+rec1.id
    rec2.id = 'Band'+str(band_val)+':'+rec2.id
    rec1.description = ''
    rec2.description = ''


    d["band_" + str(band_val)+'_1'].append(rec1)
    d["band_" + str(band_val)+'_2'].append(rec2)


    overflow = False
    for i in range(1,8):
        if len(d["band_" + str(i)+'_1']) > 200000:
            overflow = True
            break


    if overflow:
        print(overflow)
        print(cc)
        # Range (1,8) for initial band list
        # Range (1,11) for Promega band list
        for i in range(1,8):
            key =  "band_" + str(i)
            numReads = numReads + len(d[key+'_1'])
            if(simulationType=="RSEM"):
                with open(out_dir+'band'+str(i)+'_R1_sim.fq', "a+") as handle1:
                    SeqIO.write(d[key+'_1'], handle1, "fastq")
                handle1.close()
                with open(out_dir+'band'+str(i)+'_R2_sim.fq', "a+") as handle2:
                    SeqIO.write(d[key+'_2'], handle2, "fastq")
                handle2.close()
            elif(simulationType=="POLYESTER"):
                with open(out_dir+'band'+str(i)+'_R1_sim.fasta', "a+") as handle1:
                    SeqIO.write(d[key+'_1'], handle1, "fasta-2line")
                handle1.close()
                with open(out_dir+'band'+str(i)+'_R2_sim.fasta', "a+") as handle2:
                    SeqIO.write(d[key+'_2'], handle2, "fasta-2line")
                handle2.close()



            #print 'band ',i,' reads: ',len(d[key+'_1'])


        ## Clearing the memory
        d.clear()
        # Range (1,8) for initial band list
        # Range (1,11) for Promega band list
        for i in range(1,8):
            d["band_" + str(i)+'_1'] = []
            d["band_" + str(i)+'_2'] = []


        print("Written ",numReads," reads")

if not overflow:
    # Range (1,8) for initial band list
    # Range (1,11) for Promega band list
    for i in range(1,8):
        key =  "band_" + str(i)

        if(simulationType=="RSEM"):
            with open(out_dir+'band'+str(i)+'_R1_sim.fq', "a+") as handle1:
                SeqIO.write(d[key+'_1'], handle1, "fastq")
            handle1.close()
            with open(out_dir+'band'+str(i)+'_R2_sim.fq', "a+") as handle2:
                SeqIO.write(d[key+'_2'], handle2, "fastq")
            handle2.close()
        elif(simulationType=="POLYESTER"):
            with open(out_dir+'band'+str(i)+'_R1_sim.fasta', "a+") as handle1:
                SeqIO.write(d[key+'_1'], handle1, "fasta-2line")
            handle1.close()
            with open(out_dir+'band'+str(i)+'_R2_sim.fasta', "a+") as handle2:
                SeqIO.write(d[key+'_2'], handle2, "fasta-2line")
            handle2.close()

        numReads += len(d[key+'_1'])

print(cc)
print("Written ",numReads," reads")
