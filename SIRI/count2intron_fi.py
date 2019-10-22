#!/bin/python
import getopt
import copy
import re
import os
import sys
import time
import datetime
import pysam
from collections import defaultdict
from multiprocessing import Pool

def makedirs(_dir):
    try:
        os.stat(_dir)
    except:
        os.makedirs(_dir)

def get_intron2anno(gtf, lib):
    gtf_path, gtf = os.path.split(gtf.split(',')[1])
    if gtf_path == '':
        gtf_path = '.'
    intron_anno = dict()
    fp = open("%s/Intron_Annotated_%s" % (gtf_path, gtf))
    for info in fp:
        a = info.strip().split("\t")
        ann = re.sub('.*annotated_IR "|\".*', '', a[8])
        if lib == 'unstrand':
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], 'unstrand')
        else:
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], a[6])
        if (key in intron_anno):
            if (ann == "true"):
                intron_anno[key] = ann
        else:
            intron_anno[key] = ann
    fp.close()
    intron_clean = dict()
    fp = open("%s/Intron_attri_%s" % (gtf_path, gtf))
    for info in fp:
        a = info.strip().split("\t")
        C1 = re.sub('.*clean "|\".*', '', a[8])
        cc = re.sub('.*clean_simple "|\".*', '', a[8])
        if lib == 'unstrand':
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], 'unstrand')
        else:
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], a[6])
        if (key in intron_clean):
            if (C1 == "false"):
                intron_clean[key][0] = C1
        else:
            intron_clean[key] = [C1, cc]
    fp.close()
    return intron_anno, intron_clean

def bam2count(output, gtf, bam_files, read, lib, length, anchor, thread, bin_size=1000):
    bam_files = bam_files.split(',')
    sample_number = len(bam_files)
    run_cmd_list = []
    for sample_index in xrange(0, sample_number):
        bam = bam_files[sample_index]
        if bam.endswith('.bam') or bam.endswith('.sam'):
            cmd = 'python {}/parse_bam.py --gtf {} --length {} --anchor {} --bam {} -o {}/count_{} --lib {} --read {} --Total {}/Total_{}.txt'.format(
                os.path.dirname(os.path.realpath(__file__)), gtf, length, anchor, bam, output, sample_index, lib, read, output, sample_index)
            run_cmd_list.append(cmd)

    if len(run_cmd_list) > 0:
        pool = Pool(processes=int(thread))
        pool.map(os.system, run_cmd_list)
        pool.close()
        pool.join()

    intron_header = dict()
    for sample_index in xrange(0, sample_number):
        fp = open("%s/count_%s_intron.txt" % (output, sample_index))
        for info in fp:
            sp = info.strip().split("\t")
            if (intron_header.has_key(sp[0])):
                t = sp[3].split(",")
                if (t[0] == "false"):
                    intron_header[sp[0]][2] = re.sub("true,", "false,", intron_header[sp[0]][2])
                if (t[1] == "false"):
                    intron_header[sp[0]][2] = re.sub(",true", ",false", intron_header[sp[0]][2])
            else:
                intron_header[sp[0]] = sp[1:7]
        fp.close()

    fw = open('{}/intron_header'.format(output), 'w')
    fp = open('{}/count_0_intron.txt'.format(output))
    for line in fp:
        sp = line.strip().split('\t')
        fw.write('{}\t{}\n'.format(sp[0], '\t'.join(intron_header[sp[0]])))
    fw.close()
    fp.close()

    cmd = "cut -f 1 %s/count_0_exon.txt > %s/exon_header" % (output, output)
    os.system(cmd)

    cmd = "cut -f 3,4 %s/count_0_exon.txt > %s/exon_tail" % (output, output)
    os.system(cmd)

    cmd1 = "paste "
    cmd2 = "paste %s/intron_header " % output
    cmd3 = "paste %s/exon_header " % output

    for sample_index in xrange(0, sample_number):
        cmd = "cut -f 8,9,10,11,12,13 %s/count_%s_intron.txt > %s/intron_%s.val" % (output, sample_index, output, sample_index)
        os.system(cmd)
        cmd = "cut -f 2 %s/count_%s_exon.txt > %s/exon_%s.val" % (output, sample_index, output, sample_index)
        os.system(cmd)
        cmd1 += "%s/Total_%s.txt " % (output, sample_index)
        cmd2 += "%s/intron_%s.val " % (output, sample_index)
        cmd3 += "%s/exon_%s.val " % (output, sample_index)

    cmd1 += " > %s/Total.txt" % (output)
    cmd2 += " > %s/count_intron.txt" % (output)
    cmd3 += " %s/exon_tail  > %s/count_exon.txt" % (output, output)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

    return sample_number

def count2fi(gtf, lib, length, anchor, sample_number, output, update):
    intron_anno, intron_clean = get_intron2anno(gtf, lib)
    intron_count_file = open('{}/count_intron.txt'.format(output))
    result_path = output + '/../' + 'results/'
    makedirs(result_path)
    fw = open("%s/intron.PI.txt" % result_path, "w")
    fw.write(
        "Intron_id\tGene_id\tStrand\tChr\tStart\tEnd\tAnnotated\tAttributes\tInclusion_counts\tSkipping_counts\tInclusion_counts_with_intron_body\tInclusion_length\tSkipping_length\tIntron_body_length\tPI_Junction\tPI_JunctionIntron\n")
   
    for intron_count_line in intron_count_file:
        intron_count_list = intron_count_line.strip().split("\t")
        if update != "false":
            attri = intron_count_list[3].split(",")
            if (attri[0] == "false"):
                intron_clean[intron_count_list[0]][0] = "false"
            if (attri[1] == "false"):
                intron_clean[intron_count_list[0]][1] = "false"

        skp = [0] * sample_number
        inc = [0] * sample_number
        inc_intron = [0] * sample_number
        PI_J = [0] * sample_number
        PI_J_Intron = [0] * sample_number
        sk_l = length - 2 * anchor + 1
        intron_length = int(intron_count_list[6]) - int(intron_count_list[5]) + 1
        if intron_length == 0:
            continue
        in_l = length - 4 * anchor + 2 + min(intron_length, length)
       	#in_l = 2 * (length - 2 * anchor + 1)
        intron_body_length = intron_length + length - 2 * anchor + 1

        for i in range(0, sample_number):
            inc[i] = str(int(intron_count_list[i * 6 + 7]) + int(intron_count_list[i * 6 + 9]))
            inc_intron[i] = str(int(intron_count_list[i * 6 + 12]) + int(intron_count_list[i * 6 + 9]) + int(intron_count_list[i * 6 + 7]))
            skp[i] = intron_count_list[i * 6 + 11]
            if update != "false":
                if (intron_count_list[i * 6 + 8] != intron_count_list[i * 6 + 11] or intron_count_list[i * 6 + 10] !=
                    intron_count_list[i * 6 + 11]):
                    intron_clean[intron_count_list[0]][1] = "false"
            if ((inc[i] + skp[i]) > 0) & ((float(inc[i]) / in_l + float(skp[i]) / sk_l) > 0):
                PI_J[i] = str((float(inc[i]) / in_l) / (float(inc[i]) / in_l + float(skp[i]) / sk_l))
            else:
                PI_J[i] = 'NA'
            if int(inc_intron[i]) + int(skp[i]) > 0:
                PI_J_Intron[i] = str(float(inc_intron[i]) / intron_body_length / (float(inc_intron[i]) / intron_body_length + float(skp[i]) / sk_l))
            else:
                PI_J_Intron[i] = 'NA'

        attri = "U"
        if (intron_clean[intron_count_list[0]][0] == "false" and intron_clean[intron_count_list[0]][1] == "false"):
            attri = "EI"
        if (intron_clean[intron_count_list[0]][0] == "false" and intron_clean[intron_count_list[0]][1] == "true"):
            attri = "E"
        if (intron_clean[intron_count_list[0]][0] == "true" and intron_clean[intron_count_list[0]][1] == "false"):
            attri = "I"

        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            intron_count_list[0], intron_count_list[1], intron_count_list[2], intron_count_list[4],
            intron_count_list[5], intron_count_list[6], intron_anno[intron_count_list[0]], attri, ",".join(inc),
            ",".join(skp), ",".join(inc_intron), in_l,
            sk_l, intron_body_length, 
            ",".join(PI_J), ",".join(PI_J_Intron)))

    fw.close()
    intron_count_file.close()

def parse_args():
    options, args = getopt.getopt(sys.argv[1:], 'o:',
                                  ['gtf=', 'anchor=', 'lib=', 'read=', 'length=', 'bam_files=', 'output=', 'thread=', 'update='])
    gtf = ''
    read = ''
    bam_files = ''
    lib = 'unstrand'
    length = 0
    anchor = 0
    output = ''
    thread = 1
    update = 'false'

    for opt, arg in options:
        if opt in ('-o', '--output'):
            output = arg
        elif opt in ('--gtf'):
            gtf = arg
        elif opt in ('--bam_files'):
            bam_files = arg
        elif opt in ('--lib'):
            lib = arg
        elif opt in ('--read'):
            read = arg
        elif opt in ('--length'):
            length = int(arg)
        elif opt in ('--anchor'):
            anchor = int(arg)
        elif opt in ('--thread'):
            thread = arg
        elif opt in ('--update'):
            update = arg

    if not gtf or not read or not bam_files or not output or not length or not anchor:
        print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "A python program to count the reads for retained intron events for varities of junction from a series of bam file."
        print "Usage :", sys.argv[0], " --gtf: the intron gtf file and gtf file seperated by comma;"
        print "Usage :", sys.argv[0], " --length:the length of reads;"
        print "Usage :", sys.argv[0], " --anchor:the anchor length of the read;"
        print "Usage :", sys.argv[0], " --bam_files: the bam file,multiple bam file seperated by commas;"
        print "Usage :", sys.argv[0], " --lib: the library type;"
        print "Usage :", sys.argv[0], " --read: The sequencing strategy of producing reads with choices of P/S;"
        print "Usage :", sys.argv[0], ' --output: intron_id, gene_id,strand,chr,start,end,5SS inclusion counts,5SS skipping counts, 3SS includion counts,3SS skipping counts,skipping counts,intron counts;'
        print "Usage :", sys.argv[0], " --thread: the thread to get bam reads counts."
        print 'Usage :', sys.argv[0], " --update: update intron attributes based on data"
        print datetime.datetime.now()
        sys.exit()

    print 'parsing bam files...'
    sample_number = bam2count(output, gtf, bam_files, read, lib, length, anchor, thread, bin_size=1000)
    print 'converting counts to PI'
    count2fi(gtf, lib, length, anchor, sample_number, output, update)

if __name__ == "__main__":
    parse_args()
