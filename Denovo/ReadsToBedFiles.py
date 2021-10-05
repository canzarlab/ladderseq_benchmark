
# coding: utf-8

# In[1]:


# Generated bed files containing transcriptsOfOriginmore  and read positions
# Generates two bed files:
# 1. readPosition bed file : transcriptOfOrigin | startPosition (0 based) | endPosition
# 2. transcriptPosition bed file : transcriptId | startPosition (0 based) | endPosition


# In[2]:


from Bio import SeqIO
import HTSeq
import sys
import random
import numpy as np
import csv


# In[3]:


# Gtf annotation file to compute transcript length
#gtf_filename = '/algbio1/shounak/genomes/human/release_80/Homo_sapiens.GRCh38.80.gtf'
gtf_filename = sys.argv[1]


# This file has the mapping transcript idx and id
# For RSEM data
#isofrom_results_file = '/algbio1/shounak/raw/simulated/human/NA12716_7_new/NA12716_7_new.isoforms.results'
isofrom_results_file = sys.argv[2]

# paired-end input files
#in_read_files = ['/algbio1/shounak/raw/simulated/human/NA12716_7_new/new_NA12716_7_100M_1.fq','/algbio1/shounak/raw/simulated/human/NA12716_7_new/new_NA12716_7_100M_2.fq']
in_read_files = [sys.argv[3],sys.argv[4]]


# Output directory
#out_dir = '/algbio1/shounak/raw/simulated/human/NA12716_7_new/coverageGroundTruth/100Mil/'
out_dir = sys.argv[5]

# Simulation Type (RSEM or POLYESTER)
#simulationType = "RSEM"
simulationType = sys.argv[6]

# Read length
#readLength = 75
readLength = int(sys.argv[7])


# Output type
#outputType = "NotOnly_transcript"
outputType = sys.argv[8]


# In[4]:


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


# In[5]:


def get_transcript_length_dict(gtf_filename):
    gff_file = HTSeq.GFF_Reader(gtf_filename, end_included=True)
    #1. Go through the GFF file and get all exons for each transcript:
    transcript_exon_dict = {}
    for feature in gff_file:
        if feature.type == "exon":
#             print(feature.iv.strand)
#             print(feature.iv,"  ",feature.iv.length)
            transcript_id = feature.attr['transcript_id']
            if transcript_id not in transcript_exon_dict:
                transcript_exon_dict[transcript_id] = list()
            transcript_exon_dict[transcript_id].append(feature)

    #2. Go through the transcripts and calculate the length of each
    # by adding up the lengths of its exons also add the transcript start end and its strand
#     print(len(transcript_exon_dict.keys()))

    transcript_len_dict = {}

    for transcript_id in transcript_exon_dict.keys():
#         print(transcript_id)
        transcript_len_dict[transcript_id] = ()
        firstStep = True # using an arbitrary number
        # The length in the gff is 1 based but for us its always 0 based
        length = -1 # starting from -1 since the lengths in the gtf files are start from one and we want it 0 based
        numExons = 0
#         print(len(transcript_exon_dict[transcript_id]))
        for exon in transcript_exon_dict[transcript_id]:
            numExons += 1
            if(firstStep):
                firstExon = exon
                lastExon = exon
                firstStep = False
            elif(exon.iv.start < firstExon.iv.start):
                firstExon = exon
            elif(exon.iv.end > lastExon.iv.end):
                lastExon = exon
            length += exon.iv.length
#         print(transcript_id," : ",firstExon,":",lastExon,": Number: ",numExons," length: ",length)
        strand = firstExon.iv.strand
        if(numExons > 1):
            start = firstExon.iv.length - 1
            end = (length)-(lastExon.iv.length) + 1
        else:
            start = 2
            end = firstExon.iv.length - 2
        #print(start,"  ",end)
        transcript_len_dict[transcript_id] = (start,end,length,strand)

    return transcript_len_dict


# In[6]:


## the read files written here are 0 based

transcript_len_dict = get_transcript_length_dict(gtf_filename)

#print(len(transcript_len_dict))

count = 0
if (outputType != "Only_transcript"):
    if(simulationType=="RSEM"):
        transcript_idx_id_dict = get_transcript_idx_id_dict(isofrom_results_file)
        fq1 = SeqIO.parse(in_read_files[0], "fastq")
        fq2 = SeqIO.parse(in_read_files[1], "fastq")

    if(simulationType=="POLYESTER"):
        fq1 = SeqIO.parse(in_read_files[0], "fasta")
        fq2 = SeqIO.parse(in_read_files[1], "fasta")

    print("Starting with the reads.")
    ## Writing the read Position file for the first pair
    outFile = out_dir+"readPosFile.bed"
    file = open(outFile,'w')
    ## test file
#     outFileTest = out_dir+"test.txt"
#     filetest = open(outFileTest,'w')
    ## test file end

    for rec in fq1:
        # Extracting the transcriptId to which a particular read maps to
        if(simulationType=="POLYESTER"):
            transcript_id = (rec.id.split('/')[1])
            transcript_id = (transcript_id.split('.')[0])
            if(len(transcript_id)!=15):
                print("Problem")
                continue
            header = rec.description.split(' ')
            description = header[len(header)-1]
            #print(description)
            if("mate1Start" in description):
                continue
            pair1StartPos = int(((description.split(';')[1]).split(':')[1]).split('-')[0])
            pair1EndPos = int(((description.split(';')[1]).split(':')[1]).split('-')[1])
            pair2StartPos = int(((description.split(';')[2]).split(':')[1]).split('-')[0])
            pair2EndPos = int(((description.split(';')[2]).split(':')[1]).split('-')[1])
    #         print(pair1StartPos," ",pair1EndPos)
    #         print(pair2StartPos," ",pair2EndPos)
        elif(simulationType=="RSEM"):
            idx = int(rec.id.split('_')[2])
            if idx == 0:
                continue
                #Writing the transcript position file
#             print(rec.id)
            transcript_id = transcript_idx_id_dict[idx]
            readStrand = int(rec.id.split('_')[1])

#             print(transcript_id)
#             line = rec.id+"  "+str(transcript_id)+"  "+str(readStrand)+"\n"
#             filetest.write(line)
            transcriptStrand = transcript_len_dict[transcript_id][3]
            transcriptLength = transcript_len_dict[transcript_id][2]
            insertLength = rec.id.split('_')[4]
            insertLength = int(insertLength.split('/')[0])
            ## Calculating the positions of the reads based on their strandedness.
            ## Not sure if this needs to be done
            if((transcriptStrand=="+" and readStrand==0) or (transcriptStrand=="-" and readStrand==1)):
                pair1StartPos = int(rec.id.split('_')[3])
                pair1EndPos = pair1StartPos + readLength - 1
                pair2EndPos = pair1StartPos + insertLength -1
                pair2StartPos = pair2EndPos - readLength + 1

            elif((transcriptStrand=="+" and readStrand==1) or (transcriptStrand=="-" and readStrand==0)):
#                 print(rec.id)
                pair2EndPos = transcriptLength - int(rec.id.split('_')[3])
                pair2StartPos = pair2EndPos - readLength + 1
                pair1StartPos = pair2EndPos - insertLength + 1
                pair1EndPos = pair1StartPos + readLength - 1
#                 pair1EndPos = pair2EndPos - insertLength + 1
#                 pair1StartPos = pair1EndPos - readLength + 1


            if(pair1StartPos<0 or pair2StartPos<0):
                print(rec.id)

        linePair1 = str(transcript_id)+'\t'+str(pair1StartPos)+'\t'+str(pair1EndPos)+'\n'
        linePair2 = str(transcript_id)+'\t'+str(pair2StartPos)+'\t'+str(pair2EndPos)+'\n'
#         if(transcript_id == "ENST00000361714"):
#             print(rec.id)
#             print(transcriptLength)
#             print(linePair1)
#             print(linePair2)
#             print(insertLength)
#             count = count +1
        file.write(linePair1)
        file.write(linePair2)
#     filetest.close()
    file.close()
    print("Finished reads.")
#print  (count)


#Writing the transcript position file
outFileTranscriptPos = out_dir+"transcriptPosFile.bed"
fileTranscriptPos = open(outFileTranscriptPos,'w')
for key, value in transcript_len_dict.items():
    line = str(key)+'\t'+str(value[0])+'\t'+str(value[1])+'\n'
#     line = str(key)+'\t'+str(value[2])+'\n' # print this line for just the transcript lengths
#     print(line)
    fileTranscriptPos.write(line)
fileTranscriptPos.close()
