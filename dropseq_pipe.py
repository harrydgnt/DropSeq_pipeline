import sys 
import os
import argparse

# WRITTEN BY HARRY YANG - based on pipelines from Shan Sabri 

parser = argparse.ArgumentParser()	

"""
part 0 - convert fq to sam 

import os
 
FQDIR      = "/u/home/s/ssabri/project-ernst/dropseq/2016-09-08/Drop_7_C/demux/"
WRDIR      = "/u/home/s/ssabri/project-ernst/dropseq/2016-09-08/Drop_7_C/process/"
LOGDIR     = WRDIR + "/qsub_oe/"
FASTQTOSAM = "java -jar /u/home/s/ssabri/bin/picard-tools-1.138/picard.jar FastqToSam"
 
if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)
 
samples = [line.strip() for line in open(WRDIR + "/samples.txt")]
samples = ["D34", "D9", "D15"]
 
for idx, s in enumerate(samples):
    os.system("printf 'Processing: %s (%s/%s)\n' " + s + " " + str(idx+1) + " " + str(len(samples)))
    os.system("echo '#!/bin/csh' > runPicardFastqToSam_" + s)
    os.system("echo 'source ~/.bash_profile' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -cwd' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -o " + LOGDIR + "' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -e " + LOGDIR + "' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -m n' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -l h_data=16G,h_rt=8:00:00' >> runPicardFastqToSam_" + s)
    os.system("echo '" + FASTQTOSAM + " FASTQ="  + FQDIR + "/" + s + "_1.fq.gz" +
                                      " FASTQ2=" + FQDIR + "/" + s + "_2.fq.gz" +
                                      " OUTPUT=" + WRDIR + "/" + s + "_RAW_UNALIGNED.bam"     +
                                      " SAMPLE_NAME=" + s + "' >> runPicardFastqToSam_" + s)
    #os.system("qsub runPicardFastqToSam_" + s)
    #os.system("rm runPicardFastqToSam_" + s)

"""


"""
part 1 - qc filter 

#!/usr/bin/python
 
__author__  = "Shan Sabri"
 
import os
 
WRDIR = "/u/home/s/ssabri/project-ernst/dropseq/2016-09-08/Drop_7_C/process"
DST   = "/u/home/s/ssabri/project-ernst/dropseq/v4_03042016/DST"
 
samples = [line.strip() for line in open(WRDIR + "/samples.txt")]
 
for idx, s in enumerate(samples):
    print "Processing: %s (%s/%s)" % (s, str(idx+1), str(len(samples)))
 
    TEMPDIR = WRDIR + "/" + s + "_temp/"
    OUTDIR  = WRDIR + "/" + s + "_output/"
    LOGDIR  = WRDIR + "/qsub_oe"
 
    if not os.path.exists(TEMPDIR): os.makedirs(TEMPDIR)
    if not os.path.exists(OUTDIR): os.makedirs(OUTDIR)
    if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)
 
    tag_cells = DST + "/TagBamWithReadSequenceExtended " \
                      "INPUT=" + s + "_RAW_UNALIGNED.bam " \
                      "OUTPUT=" + TEMPDIR + s + "_tagged_cell.bam " \
                      "SUMMARY=" + OUTDIR + "/tagged_cell_summary.txt " \
                      "BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1"
 
    tag_molecules = DST + "/TagBamWithReadSequenceExtended " \
                          "INPUT=" + TEMPDIR + s + "_tagged_cell.bam " \
                          "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi.bam " \
                          "SUMMARY=" + OUTDIR + "/tagged_cell_umi_summary.txt " \
                          "BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1"
 
    filter_bam = DST + "/FilterBAM TAG_REJECT=XQ " \
                       "INPUT=" + TEMPDIR + s + "_tagged_cell_umi.bam " \
                       "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered.bam" \
 
    trim_starting_sequence = DST + "/TrimStartingSequence " \
                                   "INPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered.bam " \
                                   "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart.bam " \
                                   "OUTPUT_SUMMARY=" + OUTDIR + "/adapter_trimming_report.txt " \
                                   "SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"
 
    trim_poly_a = DST + "/PolyATrimmer " \
                        "INPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart.bam " \
                        "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA.bam " \
                        "OUTPUT_SUMMARY=" + OUTDIR + "/polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"
 
    sam_to_fastq = "java -Xmx500m -jar /u/home/s/ssabri/bin/picard-tools-1.138/picard.jar SamToFastq " \
                        "INPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA.bam " \
                        "FASTQ=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA.fastq " \
 
    os.system("echo '#!/bin/bash' > process_" + s)
    os.system("echo 'source ~/.bash_profile' >> process_" + s)
    os.system("echo '#$ -cwd' >> process_" + s)
    os.system("echo '#$ -o " + LOGDIR + "' >> process_" + s)
    os.system("echo '#$ -e " + LOGDIR + "' >> process_" + s)
    os.system("echo '#$ -m n' >> process_" + s)
    os.system("echo '#$ -l h_data=16G,h_rt=12:00:00' >> process_" + s)
    os.system("echo '" + tag_cells + "' >> process_" + s)
    os.system("echo '" + tag_molecules + "' >> process_" + s)
    os.system("echo '" + filter_bam + "' >> process_" + s)
    os.system("echo '" + trim_starting_sequence + "' >> process_" + s)
    os.system("echo '" + trim_poly_a + "' >> process_" + s)
    os.system("echo '" + sam_to_fastq + "' >> process_" + s)
    os.system("qsub process_" + s)
    os.system("rm process_" + s)

"""

"""
TODO - read length filter? filter reads with 30bp or less 
"""

"""
PART 2 - alignment
#!/usr/bin/python
 
__author__  = "Shan Sabri"
 
import os
 
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

"""
convert barcodes (*_1.fq.gz) files to bam (here)

#!/usr/bin/python
 
__author__  = "Shan Sabri"
 
import os
 
FQDIR      = "/u/home/s/ssabri/project-ernst/dropseq/2016-09-08/Drop_7_C/demux"
WRDIR      = "/u/home/s/ssabri/project-ernst/dropseq/2016-09-08/Drop_7_C/process"
LOGDIR     = WRDIR + "/qsub_oe/"
FASTQTOSAM = "java -jar /u/home/s/ssabri/bin/picard-tools-1.138/picard.jar FastqToSam"
 
if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)
 
samples = [line.strip() for line in open(WRDIR + "/samples.txt")]
samples = ["D24", "D15", "D12"]
 
for idx, s in enumerate(samples):
    os.system("printf 'Processing: %s (%s/%s)\n' " + s + " " + str(idx+1) + " " + str(len(samples)))
    os.system("echo '#!/bin/csh' > runPicardFastqToSam_" + s)
    os.system("echo 'source ~/.bash_profile' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -cwd' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -o " + LOGDIR + "' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -e " + LOGDIR + "' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -m n' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -l h_data=16G,h_rt=4:00:00' >> runPicardFastqToSam_" + s)
    os.system("echo '" + FASTQTOSAM + " FASTQ="  + FQDIR + "/" + s + "_1.fq.gz" +
                                      " OUTPUT=" + WRDIR + "/" + s + "_bcs.bam"     +
                                      " SAMPLE_NAME=" + s + "' >> runPicardFastqToSam_" + s)
    os.system("qsub runPicardFastqToSam_" + s)
    os.system("rm runPicardFastqToSam_" + s)

"""
"""
#!/bin/bash 
 
source ~/.bash_profile
 
# Sort barcode bams by spot ID
for f in *_bcs.bam; do echo "Sorting: ${f%%.*}"; samtools sort -n $f ${f%%.*}_nameSorted; done
# -m 10240000000
 
# Sort BED files by spot ID
for f in *_output/*_Final.bed; do b=${f##*/}; echo "Sorting: ${b%%.*}"; sort -k4 -u $f > ${b%%.*}_nameSorted.bed; done
 
# Convert sorted barcode BAM to SAM
for f in *_nameSorted.bam; do echo "Converting: ${f%%.*}"; samtools view $f > ${f%%.*}.sam; done
 
# Preform the join for each species! 
while read s; do
    echo "Joining ${s} on mm9"; 
    LANG=en_EN join -a1 -1 4 -2 1 -o 1.4 1.1 1.2 1.3 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.10 \
        <(LANG=en_EN sort -k4 ${s}_mm9_Final_nameSorted.bed) \
        <(LANG=en_EN sort -k1 ${s}_bcs_nameSorted.sam) \
       | tr ' ' '\t' \
       | awk -F $'\t' 'BEGIN {OFS = FS} {$2="MOUSE_"$2;$9="GENE:"$9;$10="ACC:"$10;$6="MOUSE_"$6;$13="BC:"substr($12,0,12);$14="UMI:"substr($12,13,8);$15="FEAT:CODING";print}' \
       > ${s}_mm9_Joined.bed; 
done < samples.txt
 
while read s; do
        echo "Joining ${s} on hg19";
        LANG=en_EN join -a1 -1 4 -2 1 -o 1.4 1.1 1.2 1.3 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.10 \
                <(LANG=en_EN sort -k4 ${s}_hg19_Final_nameSorted.bed) \
                <(LANG=en_EN sort -k1 ${s}_bcs_nameSorted.sam) \
       | tr ' ' '\t' \
       | awk -F $'\t' 'BEGIN {OFS = FS} {$2="HUMAN_"$2;$9="GENE:"$9;$10="ACC:"$10;$6="HUMAN_"$6;$13="BC:"substr($12,0,12);$14="UMI:"substr($12,13,8);$15="FEAT:CODING";print}' \
       > ${s}_hg19_Joined.bed;
done < samples.txt
 
# combine the joins
while read samples; do echo $samples; cat ${samples}_mm9_Joined.bed ${samples}_hg19_Joined.bed > ${samples}_FinalJoined.bed; done <samples.txt
 
 
# garbage collector


"""





"""
filtering out non stranded information and uniq spot/read ID using the following command 
'for f in *.bed.gz ; do echo $f ; gzcat $f | pv | parallel --pipe --keep-order --block 500M -j4 \
-q perl -ne 'use Digest::MD5 qw(md5_base64); print unless $seen{md5_base64($_)}++' | awk '$5 == $11 {print $0}' | 
/Users/shansabri/Downloads/pigz-2.3.4/pigz -p 4 > ${f%.*.*}.uniq.stranded.bed.gz ; done` 
"""





"""
Collapse the uniq.stranded.bed.gz with dropseq_processing.py
#!/usr/bin/python

from __future__ import division
from numpy import cumsum
from collections import Counter, defaultdict
from itertools import chain, combinations, product
import pandas as pd
import datetime, os, operator, gzip, warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib.pyplot as plt

__author__  = 'Shan Sabri'
__date__    = '2016/12/09'

wrdir  = '/Users/shansabri/Dropbox/PlathLab/Analyses/Justin_Langerman/2017-01-04/DropTA2/beds'

extend_knee = 5000
take_top_n_cells = 2000
barnyard_threshold = 70
barnyard_read_cutoff = 5000

input_format = 'uniq.stranded.bed.gz'

CHARS     = 'ACGT'
bc_flag   = 'BC:'
gene_flag = 'GENE:'
umi_flag  = 'UMI:'

VERBOSE = False
with_mixed = False
with_uncollapsed_dges = False

start = datetime.datetime.now()


def init_dir(files):
    print '{}\tCreating necessary directories'.format(datetime.datetime.now()-start)
    if not os.path.exists(os.path.join(wrdir, "output")):
        os.makedirs(os.path.join(wrdir, "output"))
    outdir = os.path.join(wrdir, "output")
    for f in files:
        if not os.path.exists(os.path.join(outdir, f)):
            os.makedirs(os.path.join(outdir, f))
    return outdir


def file_len(file):
    with gzip.open(file) as f:
        for i, l in enumerate(f, start=1):
            pass
    return i

def extract_barcodes(file):
    with gzip.open(file) as f:
        for idx, l in enumerate(f):
            l = l.strip().split('\t')
            obs_bc = [x for x in l if bc_flag in x][0].split(':')[1]
            yield obs_bc


def count_cell_barcodes(file):
    print datetime.datetime.now()-start, '\t- Creating barcode counter dictionary'
    barcode_counts = Counter(extract_barcodes(file))
    print('{}\t-- Identified {} barcodes corresponding to {} reads'.format(datetime.datetime.now()-start,
                                                                           len(barcode_counts), sum(barcode_counts.values())))
    return barcode_counts


def count_umis(file, barcodes):
    print '{}\t- Creating cell:[gene:transcipt] dictionary'.format(datetime.datetime.now() - start)
    cells = {}
    for idx, l in enumerate(gzip.open(file)):
        l = l.strip().split("\t")
        gene    = [x for x in l if gene_flag in x]
        obs_umi = [x for x in l if umi_flag in x]
        obs_bc  = [x for x in l if bc_flag in x]
        gene    = gene[0].split(":")[1].upper()
        obs_bc  = obs_bc[0].split(":")[1]
        obs_umi = obs_umi[0].split(":")[1]
        if obs_bc not in barcodes: continue
        if obs_bc not in cells: cells[obs_bc] = {}
        if gene not in cells[obs_bc]: cells[obs_bc][gene] = Counter()
        cells[obs_bc][gene].update([obs_umi])
    return cells


def hamming_circle(s, n):
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(CHARS) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == CHARS[r]:
                    cousin[p] = CHARS[-1]
                else:
                    cousin[p] = CHARS[r]
            yield ''.join(cousin)


def hamming_ball(s, n):
    return chain.from_iterable(hamming_circle(s, i) for i in range(n + 1))


def find_mergable(kmers, i, edits):
    if i in kmers:
        return i
    for s in hamming_ball(i, edits):
        if s in kmers:
            return s


def find_cell_barcodes_to_collapse(x, edits=1, cutoff=5):
    print '{}\t- Identifying cell barcodes to collapse'.format(datetime.datetime.now() - start)
    seen = set()
    collapsed = {}
    for k, v in x.most_common():
        if v < cutoff: continue
        true = find_mergable(seen, k, edits)
        if true:
            if VERBOSE:
                print '{}\t-- Collapsing {} ({} occurrences) with {} ({} occurrences)'.format(datetime.datetime.now() - start,
                                                                                              k, v, true, x[true])
            collapsed[k] = true
        else:
            seen.add(k)
    return collapsed


def collapse_cell_barcodes(f, dict, out):
    print '{}\t- Collapsing barcodes'.format(datetime.datetime.now() - start)
    with gzip.open(f) as i, gzip.open(out, 'w') as o:
        for idx, l in enumerate(i, start=1):
            l = l.strip().split('\t')
            obs  = [x for x in l if bc_flag in x]
            obs  = obs[0].split(':')[1]
            true = dict.get(obs)
            if true:
                if VERBOSE:
                    print '{}\t-- Replaceing less common {} with more common {} at line {}'.format(
                        datetime.datetime.now() - start, obs, true, idx)
                l[l.index(bc_flag + obs)] = bc_flag + true
            o.write('\t'.join(l) + '\n')


def collapse_umis(x, edits=1):
    print '{}\t- Figuring which UMIs to collapse'.format(datetime.datetime.now() - start)
    collapsed = {}
    for cell, gene_dict in x.iteritems():
        collapsed[cell] = {}
        for gene, umis in gene_dict.iteritems():
            collapsed[cell][gene] = len(umis)
            if len(umis) == 1: continue
            seen = set()
            for k, v in umis.most_common():
                true = find_mergable(seen, k, edits)
                if true:
                    if VERBOSE:
                        print '{}\t\t\tCan collapse {} ({} occurrences) with {} ({} occurrences) -- {}[{}]'.format(datetime.datetime.now()-start,
                                                                                                                   k, umis[k], true, umis[true], cell, gene)
                    umis[true] = umis[k] + umis[true]
                    del umis[k]
                else:
                    seen.add(k)
            collapsed[cell][gene] = len(umis)
    return collapsed


def plot_knee(x, outfile, xlim=(0,5000), h_line=0, close=True):
    print '{}\t- Plotting knee'.format(datetime.datetime.now()-start)
    total = sum(x.values())
    cum_perc = map(lambda x: round(x*100.0/total, 10), cumsum(sorted(x.values(), reverse=True)))

    if VERBOSE:
        for i in range(100, 1100, 100):
            print '{}\t-- {} cells capture {}% of total reads'.format(datetime.datetime.now()-start, i, cum_perc[i])

    # most_common = x.most_common()
    # with open(outfile, 'w') as o:
    #     o.write('Barcode\tReads\tCumulativePercent\n')
    #     for percent, kv in zip(cumPerc, most_common[:n]):
    #         key, value = kv
    #         o.write('%s\t%s\t%s\n' % (key, value, percent))

    plt.plot(cum_perc)
    plt.title(str(h_line) + " cells captured\n" + ' total.reads:' + str(total) + ', total.bcs:' + str(len(cum_perc)))
    plt.xlim(xlim)
    plt.axvline(x=h_line, ymin=0, ymax=1, color='gray', linestyle='--')
    plt.ylabel('Cumulative Fraction of Reads')
    plt.xlabel('Cellular Barcodes Sorted by Number of Reads [decending]')
    plt.savefig(outfile)
    if close: plt.close()


def get_top_cell_barcodes(x, n=1000):
    top = dict(x.most_common(n))
    return top.keys()


def split_by_species(f, species, out):
    print '{}\t- Splitting file by {} cells'.format(datetime.datetime.now()-start, species)
    with gzip.open(f) as i, gzip.open(out, 'w') as o:
        for l in i:
            l = l.strip().split('\t')
            chr = l[1]
            if species in chr:
                o.write('\t'.join(l)+'\n')


def split_by_mixed_cells(f, mixed_cells, out):
    print '{}\t- Splitting file by MIXED cells'.format(datetime.datetime.now() - start)
    with gzip.open(f) as i, gzip.open(out, 'w') as o:
        for idx, l in enumerate(i):
            l = l.strip().split('\t')
            bc = [x for x in l if bc_flag in x]
            bc = bc[0].split(':')[1]
            if bc in mixed_cells:
                o.write('\t'.join(l) + '\n')


def get_barnyard_info(file, top_barcodes):
    barnyard = {}
    with gzip.open(file) as f:
        for l in f:
            l = l.strip().split('\t')
            bc = [x for x in l if bc_flag in x]
            bc = bc[0].split(':')[1]
            if bc in top_barcodes:
                if bc in barnyard:
                    barnyard[bc] += 1
                else:
                    barnyard[bc] = 1
            else:
                barnyard[bc] = 0
    return barnyard


def identity_x_mapped_reads(f1, f2):
    print '{}\t- Identifing species X-mapped reads'.format(datetime.datetime.now() - start)
    reads, x_mapped_reads = set(), set()
    with gzip.open(f1) as m:
        for ml in m:
            ml = ml.strip().split('\t')
            reads.add(ml[0])
    with gzip.open(f2) as h:
        for hl in h:
            hl = hl.strip().split('\t')
            if hl[0] in reads:
                x_mapped_reads.add(hl[0])
    print '{}\t-- Found {} X-mapped reads'.format(datetime.datetime.now() - start, len(x_mapped_reads))
    return x_mapped_reads


def remove_x_mapped_reads(file, x_mapped_reads, out):
    with gzip.open(file) as f, gzip.open(out, 'w') as o:
        for l in f:
            l = l.strip().split('\t')
            if l[0] in x_mapped_reads: continue
            o.write('\t'.join(l) + '\n')


def plot_barnyard(barnyard, top_barcodes, outfile, threshold=75, exclude_cells_below=10000):
    data = {'x':[], 'y':[], 'label':[]}
    for label, coord in barnyard.iteritems():
        if label not in top_barcodes: continue
        if len(coord) != 2: continue
        if coord[0] < exclude_cells_below and coord[1] < exclude_cells_below: continue
        data['x'].append(coord[0])
        data['y'].append(coord[1])
        data['label'].append(label)

    plt.figure(figsize=(10,8))
    plt.title('Barnyard Plot', fontsize=20)
    plt.xlabel('MOUSE Reads', fontsize=15)
    plt.ylabel('HUMAN Reads', fontsize=15)
    plt.xlim(0, max(max(data['x']), max(data['y'])))
    plt.ylim(0, max(max(data['x']), max(data['y'])))
    plt.axvline(x=exclude_cells_below, ymin=0, ymax=exclude_cells_below/max(max(data['x']), max(data['y'])),
                color='black', linestyle='-')
    plt.axhline(y=exclude_cells_below, xmin=0, xmax=exclude_cells_below/max(max(data['x']), max(data['y'])),
                color='black', linestyle='-')

    mouse_barcodes, human_barcodes, mixed_barcodes = [], [], []
    mouse_x, human_x, mixed_x = [], [], []
    mouse_y, human_y, mixed_y = [], [], []

    for label, x, y  in zip(data['label'], data['x'], data['y']):
        x, y = int(x), int(y)
        percent_mouse = int(round((x/(x+y))*100))
        percent_human = int(round((y/(x+y))*100))
        if percent_mouse >= threshold:
            mouse_barcodes.append(label)
            mouse_x.append(x)
            mouse_y.append(y)
        elif percent_human >= threshold:
            human_barcodes.append(label)
            human_x.append(x)
            human_y.append(y)
        else:
            mixed_barcodes.append(label)
            mixed_x.append(x)
            mixed_y.append(y)

    m = plt.scatter(x=mouse_x, y=mouse_y, label=mouse_barcodes, color="blue")
    h = plt.scatter(x=human_x, y=human_y, label=human_barcodes, color="red")
    x = plt.scatter(x=mixed_x, y=mixed_y, label=mixed_barcodes, color="gray")
    plt.legend((m, h, x), ('Mouse: ' + str(len(mouse_barcodes)),
                           'Human: ' + str(len(human_barcodes)),
                           'Mixed: ' + str(len(mixed_barcodes))),
               scatterpoints=1, loc=1, fontsize=8)
    plt.savefig(outfile)
    plt.close()

    return mouse_barcodes, human_barcodes, mixed_barcodes


def write_dge(dge, out):
    print "{}\t- Writing DGE".format(datetime.datetime.now()-start)
    df = pd.DataFrame(dge).fillna(0)
    df.to_csv(out, sep="\t", compression='gzip')


def write_lib_info(x, outfile):
    print "{}\t- Writing summary".format(datetime.datetime.now() - start)
    with open(outfile, "wb") as o:
        o.write('\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ("Total.reads", "Total.reads.collapsed",
                                                        "Total.reads.collapsed.mouse", "Total.reads.collapsed.human",
                                                        "Total.uncollapsed.UMIs", "Total.collapsed.UMIs",
                                                        "Total.reads.collapsed.mouse.cleaned", "Total.reads.collapsed.human.cleaned"))
        for k, v in x.iteritems():
            o.write('%s\t%s\n' % (k, "\t".join(str(i) for i in v)))

def clean_up(files):
    print "{}\t- Cleaning up".format(datetime.datetime.now() - start)
    for f in files:
        os.remove(f)


def pipeline():
    files = [os.path.splitext(f)[0].split(".")[0] for f in os.listdir(wrdir) if f.endswith(input_format)]
    outdir = init_dir(files)

    for f in files:
        print '{}\tProcessing {}'.format(datetime.datetime.now()-start, f)

        lib_info = defaultdict(list)

        # STEP 0: Count the number of uncollapsed transcripts
        total_uncollapsed_transcripts = file_len(os.path.join(wrdir, '.'.join((f, input_format))))
        # print "total_uncollapsed_transcripts: " + total_uncollapsed_transcripts

        # STEP 1: Count the raw occurences of each barcode
        barcode_counts = count_cell_barcodes(os.path.join(wrdir, '.'.join((f, input_format))))
        for k, v in barcode_counts.iteritems(): lib_info[k].append(v)

        # STEP 2: Plot knee
        plot_knee(barcode_counts, outfile=os.path.join(outdir, f, '.'.join((f, "knee", "uncollapsed", "pdf"))),
                  xlim=(0, extend_knee), h_line=take_top_n_cells ,close=False)

        # STEP 3: Identify which cell barcodes to collapse
        barcodes_to_collapse = find_cell_barcodes_to_collapse(barcode_counts, edits=2, cutoff=50)

        # STEP 4: Collapse on the identified barcodes
        collapse_cell_barcodes(f=os.path.join(wrdir, '.'.join((f, input_format))),
                               dict=barcodes_to_collapse,
                               out=os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))))

        # STEP 5: sort | uniq on the collapsed 20-mer
        uniq_20mers = set()
        with gzip.open(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format)))) as collapsed_bed:
                for idx, l in enumerate(collapsed_bed):
                    l = l.strip().split('\t')
                    uniq_20mers.add(l[11])
        with open(os.path.join(outdir, f, '.'.join((f, "uniq.fraction.txt"))), 'w') as o:
            o.write('Uncollapsed_Transcripts\tCollapsed_Transcripts\tCollapse_Rate\n')
            o.write('{}\t{}\t{}\n'.format(str(total_uncollapsed_transcripts), str(len(uniq_20mers)),
                                          str(len(uniq_20mers)/total_uncollapsed_transcripts)))

        collapsed_barcode_counts = count_cell_barcodes(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))))
        for k, v in collapsed_barcode_counts.iteritems(): lib_info[k].append(v)

        # STEP 5: Plot collapsed knee
        plot_knee(collapsed_barcode_counts, outfile=os.path.join(outdir, f, '.'.join((f, "knee", "collapsed", "pdf"))),
                  xlim=(0,extend_knee), h_line=take_top_n_cells, close=True)

        # STEP 6: Get most abundant barcodes from collapsed cell barcode counts dictionary
        top_collapsed_cell_barcodes = get_top_cell_barcodes(x=collapsed_barcode_counts, n=take_top_n_cells)

        # STEP 7: Split collapsed bed file by species and count number of reads within the top collapsed cell barcodes
        species = ["MOUSE", "HUMAN"]
        master_barnyard = {}
        for s in species:
            split_by_species(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))), species=s,
                             out=os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))))
            species_barcode_counts = count_cell_barcodes(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))))
            for k, v in species_barcode_counts.iteritems(): lib_info[k].append(v)
            species_barnard_info = get_barnyard_info(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))), top_barcodes=top_collapsed_cell_barcodes)
            for k, v in sorted(species_barnard_info.iteritems(), key=operator.itemgetter(1), reverse=True):
                master_barnyard.setdefault(k, []).append(v)

        # STEP 8: Plot the barnyard
        plot_barnyard(master_barnyard, top_barcodes=top_collapsed_cell_barcodes,
                      outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "pdf"))),
                      threshold=barnyard_threshold, exclude_cells_below=0)

        mouse, human, mixed = plot_barnyard(master_barnyard, top_barcodes=top_collapsed_cell_barcodes,
                                            outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "cutoff", "pdf"))),
                                            threshold=barnyard_threshold, exclude_cells_below=barnyard_read_cutoff)

        # STEP 9: Write out bed files with MIXED cells
        if with_mixed:
            split_by_mixed_cells(os.path.join(outdir, f, '.'.join((f, "collapsed", input_format))), mixed_cells=mixed,
                                 out=os.path.join(outdir, f, '.'.join((f, "collapsed", "MIXED", input_format))))


        # STEP 9: Collapse UMIs and write out DGEs in a species-specific manner
        species = {'mouse': mouse, 'human': human}
        if with_mixed:
            species['mixed'] = mixed

        for s, v in species.iteritems():
            print '{}\t- Computing DGE for {} - {}'.format(datetime.datetime.now() - start, f, s)
            umi_counts = count_umis(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))), v)

            umi_counts_uncollapsed = {}
            for cell, gene_dict in umi_counts.iteritems():
                umi_counts_uncollapsed[cell] = {}
                for gene, umis in gene_dict.iteritems():
                    umi_counts_uncollapsed[cell][gene] = sum(umis.values())
            for k, v in umi_counts_uncollapsed.iteritems():  lib_info[k].append(sum(v.values()))

            umi_counts_collapsed = collapse_umis(umi_counts, edits=1)
            for k, v in umi_counts_collapsed.iteritems():  lib_info[k].append(sum(v.values()))

            write_dge(umi_counts_collapsed, os.path.join(outdir, f, '.'.join((f, "collapsed", s, "dge", "tsv.gz"))))
            if with_uncollapsed_dges:
                write_dge(umi_counts, os.path.join(outdir, f, '.'.join((f, "uncollapsed", s, "dge.showUMIs", "tsv.gz"))))
                write_dge(umi_counts_uncollapsed, os.path.join(outdir, f, '.'.join((f, "uncollapsed", s, "dge", "tsv.gz"))))


        # STEP 10: Remove cross-mapped reads from species collapsed files
        x_mapped_reads = identity_x_mapped_reads(f1=os.path.join(outdir, f, '.'.join((f, "collapsed.MOUSE", input_format))),
                                                 f2=os.path.join(outdir, f, '.'.join((f, "collapsed.HUMAN", input_format))))
        cleaned_master_barnyard = {}
        for s in ["MOUSE", "HUMAN"]:
            remove_x_mapped_reads(os.path.join(outdir, f, '.'.join((f, "collapsed", s, input_format))), x_mapped_reads,
                                  out=os.path.join(outdir, f, '.'.join((f, "collapsed", s, "cleaned", input_format))))
            cleaned_species_barcode_counts = count_cell_barcodes(os.path.join(outdir, f, '.'.join((f, "collapsed", s, "cleaned", input_format))))
            for k, v in cleaned_species_barcode_counts.iteritems(): lib_info[k].append(v)

            cleaned_species_barnard_info = get_barnyard_info(os.path.join(outdir, f, '.'.join((f, "collapsed", s, "cleaned", input_format))),
                                                             top_barcodes=top_collapsed_cell_barcodes)
            for k, v in sorted(cleaned_species_barnard_info.iteritems(), key=operator.itemgetter(1), reverse=True):
                cleaned_master_barnyard.setdefault(k, []).append(v)

        # STEP 11: Plot the barnyard
        plot_barnyard(cleaned_master_barnyard, top_barcodes=top_collapsed_cell_barcodes,
                      outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "cleaned", "pdf"))),
                      threshold=barnyard_threshold, exclude_cells_below=0)

        plot_barnyard(cleaned_master_barnyard, top_barcodes=top_collapsed_cell_barcodes,
                      outfile=os.path.join(outdir, f, '.'.join((f, "barnyard", "collapsed", "cleaned", "cutoff", "pdf"))),
                      threshold=barnyard_threshold, exclude_cells_below=barnyard_read_cutoff)

        # STEP 13: Write out library info
        write_lib_info(lib_info, os.path.join(outdir, f, '.'.join((f, "summary", "tsv")))) # not sure how to have 0 place holders for rows with no entry?

        # STEP 14: garbage collect
        # delete_these_files = [os.path.join(outdir, f, '.'.join((f, "knee.uncollapsed.pdf"))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed", input_format)))]
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.MOUSE", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.HUMAN", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.MIXED", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.MOUSE.cleaned", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "collapsed.HUMAN.cleaned", input_format))),
        #                       os.path.join(outdir, f, '.'.join((f, "uncollapsed.mouse.dge.tsv.gz"))),
        #                       os.path.join(outdir, f, '.'.join((f, "uncollapsed.human.dge.tsv.gz"))),
        #                       os.path.join(outdir, f, '.'.join((f, "uncollapsed.mixed.dge.tsv.gz")))]
        # clean_up(delete_these_files)


if __name__ == '__main__':

    print 'START: ', datetime.datetime.now()
    pipeline()
    print 'FINISH: ', datetime.datetime.now()
"""





		