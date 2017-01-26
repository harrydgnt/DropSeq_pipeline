#!/usr/bin/python

import os
import sys 
from Bio import SeqIO

dev_mode = True
if dev_mode == True:
    args.qsub = False
    args.clean_up = False
    args.verbose = True

def filter_short_reads(input_fastq):
    
    count = 0
    temp_string = ''
    long_sequences = []
    num_short_read = 0

    with open(input_fastq,'r') as lines:
        sequence_length = 0
        for line in lines:
            if count % 4 == 0: # if first line - read name 
                temp_string = temp_string + line 
            elif count % 4 == 1: # if actual sequence 
                temp_string = temp_string + line 
                sequence_length = len(line)
            elif count % 4 == 2: # @ 
                temp_string = temp_string + line 
            elif count % 4 == 3: # base quality - terminate 
                temp_string = temp_string + line
                if sequence_length >= 30: 
                    long_sequences.append(temp_string)
                else:
                    num_short_read += 1
                temp_string = '' 
            count += 1
            if count % 400000 == 0: 
                print "processed %i reads" % (count/4)


    return long_sequences

def filter_only():
    sample_files='/u/home/h/harryyan/scratch/dropTA/data/filtered_fastq/filtered_fastq.samples' # list of files in fastq
    output_dir = '/u/home/h/harryyan/scratch/dropTA/data/filtered_fastq/' # the filtered files will be written her
    with open(sample_files, 'r') as samples:
        for sample in samples:
            print sample
            sample_name = sample.split('.fastq')[0]
            sample_sequences_after_len_filter = filter_short_reads(sample.rstrip())

            output_file = output_dir + sample_name + '_filtered.fastq'
            with open(output_file, 'w') as output:
                for sequence in sample_sequences_after_len_filter:
                    output.write(sequence)

def test():
    filter_only()

test()

REF_DIR = '' 
WORK_DIR = '' 


def align(input_fastq, sample_name):
    # print sample_name
    # TEMP_DIR = WORK_DIR + "/" + sample_name + "_TEMP/"
    # OUTPUT_DIR = WORK_DIR + "/" + sample_name + "_OUT/"
    ALIGN_DIR = WORK_DIR + "/" + "alignment/"
    if not os.path.exists(ALIGN_DIR): 
        os.makedirs(ALIGN_DIR)
        os.makedirs(ALIGN_DIR + "human/")
        os.makedirs(ALIGN_DIR + "mouse/")
    species_list.append('human')
    species_list.append('mouse')

    tophat_option = " --b2-D 20 --b2-R 3 --b2-N 1 --b2-L 20 --b2-i S,1,0.50 " # TODO - make it an option
    for species in species_list:
        # for mouse
        if species == 'mouse':
            bowtie_reference = REF_DIR + species + "/" + "Bowtie2Index/genome" # TODO - make this reference format
            output_dir = ALIGN_DIR + "mouse/" + sample_name
        elif species == 'human':
            bowtie_reference = REF_DIR + species + "/" + "Bowtie2Index/genome" # TODO - make this reference format
            output_dir = ALIGN_DIR + "human/" + sample_name

        align = "tophat " + bowtie_reference + input_fastq + \
                         "-o" + output_dir + " " + tophat_options
        align = align + "\n" + "mv " + output_dir + "accepted_hits.bam " + ALIGN_DIR + "/" + sample_name + "_" + species + ".bam"
        mapped_bam = ALIGN_DIR + "/" + sample_name + "_" + species + ".bam"
        if args.clean_up == True: 
            align = align + "\n" + "rm -rf %s" % (output_dir)

        if args.verbose == True:
            print "Mapping %s to %s genome." % (sample_name, species)

        if args.qsub == True:
            # TODO - write command
            # TODO - qsub that command
            continue
        else:
            os.system(align)

        if args.verbose == True:
            print "Finished mapping %s to %s genome." % (sample_name, species)


"""
 
WRDIR = "/u/home/s/ssabri/project-ernst/dropseq/2016-09-08/Drop_7_C/process"
REFDIR = "/u/home/s/ssabri/project-ernst/ref_data"
 
samples = [line.strip() for line in open(WRDIR + "/samples.txt")]
 
for idx, s in enumerate(samples):
 
    print "Processing: %s (%s/%s)" % (s, str(idx+1), str(len(samples)))
 
    TEMPDIR = WRDIR + "/" + s + "_temp/"
    OUTDIR  = WRDIR + "/" + s + "_output/"
    LOGDIR  = WRDIR + "/qsub_oe"
    INPUT = TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA_rmShort.fastq"
    bed = ""
 
    for GENOMEDIR in [REFDIR + "/Mus_musculus/UCSC/mm9", REFDIR + "/Homo_sapiens/UCSC/hg19"]:
        if not os.path.exists(TEMPDIR): os.makedirs(TEMPDIR)
        if not os.path.exists(OUTDIR): os.makedirs(OUTDIR)
        if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)
 
        id = GENOMEDIR.split("/")[9]
 
    if id == "hg19": 
        bed = "/u/home/s/ssabri/project-ernst/ref_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.bed"
 
    if id == "mm9":
        bed = "/u/home/s/ssabri/project-ernst/ref_data/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/mm9.bed"
 
        ANNOTATION = WRDIR + "/" + id + "_longest_isoform.bed"
        OUTPUT_ALIGNED = OUTDIR + s + "_" + id + "_AlignedSorted"
        OUTPUT_BED = OUTDIR + s + "_" + id + "_AlignedSorted.bed"
        INTERSECT_OUTPUT = OUTDIR + s + "_" + id + "_Final.bed"
    DIST_OUT = OUTDIR + s + "_" + id + "_read_distribution.txt"
 
        align = "bowtie2 --time --threads 5 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" + \
                         " -x " + GENOMEDIR + "/Sequence/Bowtie2Index/genome" + \
                         " -U " + INPUT + \
                         " -S " + TEMPDIR + s + "_" + id + "_" + "aligned.sam"
 
        bamtobed = "bedtools bamtobed -i " + OUTPUT_ALIGNED + ".bam > " + OUTPUT_BED
 
    dist = "/u/home/s/ssabri/bin/RSeQC-2.6.3/scripts/read_distribution.py -i " + OUTPUT_ALIGNED + ".bam -r " + bed + " > " + DIST_OUT
"""
# /u/home/h/harryyan/scratch/tools/RSeQC-2.6.4/scripts/read_distribution.py -i [sorted_bam] -r [bed] > 


"""
        map = "bedtools intersect -bed -u -wa -a " + OUTPUT_BED + " -b " + ANNOTATION + " > " + OUTDIR + id + "_temp.bed"
         
    join = "bedtools intersect -bed -wa -wb -a " + ANNOTATION + " -b " + OUTDIR + id + "_temp.bed > " + INTERSECT_OUTPUT
     
    intersect = "bedtools intersect -wa -wb -a " + OUTPUT_BED + " -b " + ANNOTATION + " > " + INTERSECT_OUTPUT
 
        os.system("echo '#!/bin/bash' > align_" + s + "_" + id)
        os.system("echo 'source ~/.bash_profile' >> align_" + s + "_" + id)
        os.system("echo '#$ -cwd' >> align_" + s + "_" + id)
        os.system("echo '#$ -o " + LOGDIR + "' >> align_" + s + "_" + id)
        os.system("echo '#$ -e " + LOGDIR + "' >> align_" + s + "_" + id)
        os.system("echo '#$ -m n' >> align_" + s + "_" + id)
        os.system("echo '#$ -l h_data=8G,h_rt=12:00:00' >> align_" + s + "_" + id)
        os.system("echo '#$ -pe shared 5' >> align_" + s + "_" + id)
        os.system("echo 'module load samtools' >> align_" + s + "_" + id)
    os.system("echo '" + align + "' >> align_" + s + "_" + id)
        os.system("echo 'samtools view -bS " + TEMPDIR + s + "_" + id + "_" + "aligned.sam | samtools sort - " + OUTPUT_ALIGNED + "' >> align_" + s + "_" + id)
        os.system("echo '" + dist + "' >> align_" + s + "_" + id)
    os.system("echo '" + bamtobed + "' >> align_" + s + "_" + id)
    os.system("echo '" + intersect + "' >> align_" + s + "_" + id) # correct
        os.system("qsub align_" + s + "_" +id)
        os.system("rm align_" + s + "_" +id)
    """
