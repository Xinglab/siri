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

def read_file(filename):
    with open(filename) as fp:
        for line in fp:
            yield line.strip().split('\t')


def get_read_dictionary(gtf_file, lib, bin_size=1000):
    dict_count = defaultdict(lambda: [])  # the dict to store the intron count, first is the inclusion at the left side,
    # second is the skipping
    # at the left side, third is inclusion at right side, fourth is the skipping at right side, fifth is skipping counts,
    # sixth is the intron counts, at the end of the three column are geneID, geneStrand, and clean introns based on bam files
    dict_intron_left = defaultdict(lambda: [])  # the dict to store left edge position of intron
    dict_intron_right = defaultdict(lambda: [])  # the dict to store right edge position of intron
    dict_intron_left_right = defaultdict(lambda: [])  # the dict to store the two edge position of intron
    dict_intron_window = defaultdict(lambda: [])  # the dict to store the introns in a certain window
    for sp in read_file(gtf_file):
        if lib == 'unstrand':
            key = (sp[0], int(sp[3]), int(sp[4]), 'unstrand')
        else:
            key = (sp[0], int(sp[3]), int(sp[4]), sp[6])
        gene_id = sp[8].split('"')[1]
        if key in dict_count:
            gene_id_list = dict_count[key][6].split(',')
            if gene_id in gene_id_list:
                continue
            dict_count[key][6] += ',' + gene_id
            dict_count[key][7] += ',' + sp[6]
        else:
            dict_count[key] = [0] * 6 + [gene_id, sp[6], 'true', 'true', 'true']
            dict_intron_left[(sp[0], int(sp[3]))].append(key)
            dict_intron_right[(sp[0], int(sp[4]))].append(key)
            dict_intron_left_right[(sp[0], int(sp[3]), int(sp[4]))].append(key)
            index_start = int(sp[3]) / bin_size
            index_end = int(sp[4]) / bin_size
            for i in xrange(index_start, index_end + 1):
                dict_intron_window[(sp[0], i)].append(key)
    return dict_count, dict_intron_left, dict_intron_right, dict_intron_left_right, dict_intron_window


def get_exon_gene_dictionary(gtf_file, bin_size=1000):
    dict_gene_count = {}  # dict to store gene count
    dict_gene_count_windows = defaultdict(lambda: [])  # dict to store gene count by window
    dict_gene2exon = defaultdict(lambda: [])  # dict map from gene to exon
    dict_exon2gene = {}  # dict map from exon to gene
    for sp in read_file(gtf_file):
        if len(sp) < 9 or sp[2] != 'exon':
            continue
        gene_id = sp[8].split('"')[1]
        exon = (int(sp[3]), int(sp[4]))
        dict_gene2exon[gene_id].append(exon)
        dict_gene_count[gene_id] = [0]
        dict_gene_count[gene_id].append(sp[6])
        dict_gene_count[gene_id].append(sp[0])
    for gene in dict_gene2exon:  # merge the overlapping exons
        dict_gene2exon[gene] = list(set(dict_gene2exon[gene]))
        exonStart = []
        exonEnd = {}
        for exon_pos in dict_gene2exon[gene]:
            exonStart.append(exon_pos[0])
            if exon_pos[0] in exonEnd and exon_pos[1] > exonEnd[exon_pos[0]]:
                exonEnd[exon_pos[0]] = exon_pos[1]
            if exon_pos[0] not in exonEnd:
                exonEnd[exon_pos[0]] = exon_pos[1]
        exonStart = list(set(exonStart))
        if len(exonStart) == 1:
            dict_gene2exon[gene] = [(exonStart[0], exonEnd[exonStart[0]])]
        else:
            exonStart = sorted(exonStart)
            exonStart1 = []
            exonEnd1 = {}
            s_1, e_1 = exonStart[0], exonEnd[exonStart[0]]
            for i in range(1, len(exonStart)):
                s_2 = exonStart[i]
                s_min = min(s_2, s_1)
                e_2 = exonEnd[s_2]
                e_max = max(e_1, e_2)
                if e_max - s_min - 1 <= e_1 - s_1 + e_2 - s_2:
                    s_1 = s_min
                    e_1 = e_max
                else:
                    exonStart1.append(s_1)
                    exonEnd1[s_1] = e_1
                    s_1 = s_2
                    e_1 = e_2
                if i == len(exonStart) - 1:
                    exonStart1.append(s_1)
                    exonEnd1[s_1] = e_1
            dict_gene2exon[gene] = []
            for s_pos in exonStart1:
                dict_gene2exon[gene].append((s_pos, exonEnd1[s_pos]))
    for gene in dict_gene2exon:
        for exon_pos in dict_gene2exon[gene]:
            dict_exon2gene[(dict_gene_count[gene][2], exon_pos[0], exon_pos[1])] = gene
            index_start = exon_pos[0] / bin_size
            index_end = exon_pos[1] / bin_size
            for i in range(index_start, index_end + 1):
                dict_gene_count_windows[(dict_gene_count[gene][2], i)].append(
                    (dict_gene_count[gene][2], exon_pos[0], exon_pos[1]))
    return dict_gene_count, dict_gene_count_windows, dict_gene2exon, dict_exon2gene


def parse_bam_file(Total, output, gtf, bam_file, read, lib, length, anchor, bin_size=1000):
    dict_count, dict_intron_left, dict_intron_right, dict_intron_left_right, dict_intron_window = get_read_dictionary(
        gtf.split(',')[0], lib, bin_size)
    dict_gene_count, dict_gene_count_windows, dict_gene2exon, dict_exon2gene = get_exon_gene_dictionary(
        gtf.split(',')[1], bin_size)
    samfile = pysam.AlignmentFile(bam_file, 'rb')
    sam_list = samfile.fetch()
    total_reads = 0
    for reads_index, iters in enumerate(sam_list):
        if reads_index > 0 and reads_index % 100000 == 0:
            print '...parsing reads {}'.format(reads_index)
        if iters.get_tag('NH') != 1:
            continue
        if 'D' in iters.cigarstring or 'H' in iters.cigarstring or 'S' in iters.cigarstring or 'I' in iters.cigarstring:
            continue
        if read == 'P':
            if (iters.flag / 2) % 2 == 0:
                continue
        elif read == 'S':
            if not re.search(str(iters.flag), '0,16'):
                continue
        total_reads += 1
        if lib != 'unstrand':
            ss = (iters.flag) / 16 % 2
            if read == 'S' and lib == 'first':
                if ss == 0:
                    strand = "-"
                else:
                    strand = "+"
            if read == 'S' and lib == 'second':
                if ss == 0:
                    strand = "+"
                else:
                    strand = "-"
            if read == 'P' and lib == 'first':
                f = int(iters.flag) / 64 % 2
                s = int(iters.flag) / 128 % 2
                if ss == 0 and f == 1:
                    strand = "-"
                if ss == 0 and s == 1:
                    strand = "+"
                if ss == 1 and f == 1:
                    strand = "+"
                if ss == 1 and s == 1:
                    strand = "-"
            if read == 'P' and lib == 'second':
                f = int(iters.flag) / 64 % 2
                s = int(iters.flag) / 128 % 2
                if ss == 0 and f == 1:
                    strand = "+"
                if ss == 0 and s == 1:
                    strand = "-"
                if ss == 1 and f == 1:
                    strand = "-"
                if ss == 1 and s == 1:
                    strand = "+"
        match_split = iters.cigarstring.split('M')
        junction_split = iters.cigarstring.split('N')
        junction_split_length = len(junction_split)
        index = (iters.get_reference_positions()[0] + 1) / bin_size
        chrom = samfile.get_reference_name(iters.reference_id)
        reads_start_position = iters.get_reference_positions()[0]
        if (chrom, index) in dict_gene_count_windows:
            for key in dict_gene_count_windows[(chrom, index)]:
                if key[1] <= reads_start_position + 1 and key[2] >= reads_start_position + 1:
                    if lib == 'unstrand':
                        dict_gene_count[dict_exon2gene[key]][0] += 1
                    elif strand in dict_gene_count[dict_exon2gene[key]][1]:
                        dict_gene_count[dict_exon2gene[key]][0] += 1
        if len(match_split) == 2:
            for pos in xrange(reads_start_position + 1 + anchor,
                              reads_start_position + 1 + int(match_split[0]) - anchor + 1):
                if (chrom, pos) in dict_intron_left:
                    for key in dict_intron_left[(chrom, pos)]:
                        if lib == 'unstrand':
                            dict_count[key][0] += 1
                            continue
                        elif strand in dict_count[key][7]:
                            dict_count[key][0] += 1
                            continue
            for pos in xrange(reads_start_position + 1 + anchor - 1,
                              reads_start_position + 1 + int(match_split[0]) - anchor):
                if (chrom, pos) in dict_intron_right:
                    for key in dict_intron_right[(chrom, pos)]:
                        if lib == 'unstrand':
                            dict_count[key][2] += 1
                        elif strand in dict_count[key][7]:
                            dict_count[key][2] += 1
            if (chrom, index) in dict_intron_window:
                for key in dict_intron_window[(chrom, index)]:
                    if key[1] <= reads_start_position + 1 and key[2] >= reads_start_position + 1 + length - 1:
                        if lib == 'unstrand':
                            dict_count[key][5] += 1
                        elif strand in dict_count[key][7]:
                            dict_count[key][5] += 1
        start_position = reads_start_position + 1
        for i in xrange(0, junction_split_length - 1):
            n1 = int(junction_split[i].split('M')[0])
            n2 = int(junction_split[i].split('M')[1])
            n3 = int(junction_split[i + 1].split('M')[0])
            if n1 >= anchor and n3 >= anchor:
                index_1 = (start_position + n1 - 1) / bin_size
                index_2 = (start_position + n1 + n2) / bin_size
                if (chrom, index_1) in dict_intron_window:
                    for key in dict_intron_window[(chrom, index_1)]:
                        if key[1] < start_position + n1 and key[2] > start_position + n1:
                            if lib == 'unstrand':
                                dict_count[key][8] = 'false'
                                dict_count[key][10] = 'false'
                            elif strand in dict_count[key][7]:
                                dict_count[key][8] = 'false'
                                dict_count[key][10] = 'false'
                if (chrom, index_2) in dict_intron_window:
                    for key in dict_intron_window[(chrom, index_2)]:
                        if key[1] < start_position + n1 + n2 - 1 and key[2] > start_position + n1 + n2 - 1:
                            if lib == 'unstrand':
                                dict_count[key][8] = 'false'
                                dict_count[key][10] = 'false'
                            elif strand in dict_count[key][7]:
                                dict_count[key][8] = 'false'
                                dict_count[key][10] = 'false'
                if (chrom, start_position + n1) in dict_intron_left:
                    for key in dict_intron_left[(chrom, start_position + n1)]:
                        if lib == 'unstrand':
                            dict_count[key][1] += 1
                        elif strand in dict_count[key][7]:
                            dict_count[key][1] += 1
                if (chrom, start_position + n1 + n2 - 1) in dict_intron_right:
                    for key in dict_intron_right[(chrom, start_position + n1 + n2 - 1)]:
                        if lib == 'unstrand':
                            dict_count[key][3] += 1
                        elif strand in dict_count[key][7]:
                            dict_count[key][3] += 1
                if (chrom, start_position + n1, start_position + n1 + n2 - 1) in dict_intron_left_right:
                    for key in dict_intron_left_right[(chrom, start_position + n1, start_position + n1 + n2 - 1)]:
                        if lib == 'unstrand':
                            dict_count[key][4] += 1
                        elif strand in dict_count[key][7]:
                            dict_count[key][4] += 1
            if n1 - anchor + 1 > anchor:
                _start = start_position
                for pos in xrange(_start + anchor, _start + n1 - anchor + 1):
                    if (chrom, pos) in dict_intron_left:
                        for key in dict_intron_left[(chrom, pos)]:
                            if lib == 'unstrand':
                                dict_count[key][0] += 1
                            elif strand in dict_count[key][7]:
                                dict_count[key][0] += 1
                for pos in xrange(_start + anchor - 1, _start + n1 - anchor):
                    if (chrom, pos) in dict_intron_right:
                        for key in dict_intron_right[(chrom, pos)]:
                            if lib == 'unstrand':
                                dict_count[key][2] += 1
                            elif strand in dict_count[key][7]:
                                dict_count[key][2] += 1
            if n3 - anchor + 1 > anchor:
                _start = start_position + n1 + n2
                for pos in xrange(_start + anchor, _start + n3 - anchor + 1):
                    if (chrom, pos) in dict_intron_left:
                        for key in dict_intron_left[(chrom, pos)]:
                            if lib == 'unstrand':
                                dict_count[key][0] += 1
                            elif strand in dict_count[key][7]:
                                dict_count[key][0] += 1
                for pos in xrange(_start + anchor - 1, _start + n3 - anchor):
                    if (chrom, pos) in dict_intron_right:
                        for key in dict_intron_right[(chrom, pos)]:
                            if lib == 'unstrand':
                                dict_count[key][2] += 1
                            elif strand in dict_count[key][7]:
                                dict_count[key][2] += 1
            start_position += n1 + n2
    samfile.close()

    # write to file
    fw = open('{0}_intron.txt'.format(output), 'w')
    for gtf_sp in read_file(gtf.split(',')[0]):
        if lib == 'unstrand':
            key = (gtf_sp[0], int(gtf_sp[3]), int(gtf_sp[4]), 'unstrand')
        else:
            key = (gtf_sp[0], int(gtf_sp[3]), int(gtf_sp[4]), gtf_sp[6])
        if dict_count[key][9] == 'true':
            key_str = '{0}_{1}_{2}_{3}'.format(key[0], key[1], key[2], key[3])
            if '+' in dict_count[key][7]:
                fw.write('{0}\t{1}\t{2}\t{3},{4}\t{5}\t{6}\t{7}\t{8}\n'.format(key_str, dict_count[key][6],dict_count[key][7], dict_count[key][8],dict_count[key][10], gtf_sp[0],gtf_sp[3], gtf_sp[4],'\t'.join([str(x) for x in dict_count[key][0:6]])))
            else:
                fw.write('{0}\t{1}\t{2}\t{3},{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(key_str, dict_count[key][6],dict_count[key][7], dict_count[key][8],dict_count[key][10], gtf_sp[0],gtf_sp[3], gtf_sp[4], dict_count[key][2],dict_count[key][3], dict_count[key][0],dict_count[key][1], dict_count[key][4],dict_count[key][5]))
            dict_count[key][9] = 'false'
    fw.close()
    fw = open(Total, 'w')
    fw.write(str(total_reads))
    fw.close()
    fw = open('{0}_exon.txt'.format(output), 'w')
    for gene in dict_gene_count:
        fw.write('{0}\t{1}\n'.format(gene, '\t'.join([str(x) for x in dict_gene_count[gene]])))
    fw.close()


def parse_args():
    options, args = getopt.getopt(sys.argv[1:], 'o:',
                                  ['gtf=', 'anchor=', 'lib=', 'read=', 'length=', 'bam=', 'output=', 'Total='])
    gtf = ''
    read = ''
    bam = ''
    lib = ''
    length = 0
    anchor = 0
    output = ''
    Total = ''

    for opt, arg in options:
        if opt in ('-o', '--output'):
            output = arg
        elif opt in ('--gtf'):
            gtf = arg
        elif opt in ('--bam'):
            bam = arg
        elif opt in ('--lib'):
            lib = arg
        elif opt in ('--read'):
            read = arg
        elif opt in ('--length'):
            length = int(arg)
        elif opt in ('--anchor'):
            anchor = int(arg)
        elif opt in ('--Total'):
            Total = arg

    if not gtf or not read or not bam or not output or not length or not anchor or not Total:
        print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "A python program to count the reads for retained intron events for varities of junction from a series of bam file."
        print "Usage :", sys.argv[0], " --gtf: the intron gtf file and gtf file seperated by comma;"
        print "Usage :", sys.argv[0], " --length:the length of reads;"
        print "Usage :", sys.argv[0], " --anchor:the anchor length of the read;"
        print "Usage :", sys.argv[0], " --bam: the bam file,multiple bam file seperated by commas;"
        print "Usage :", sys.argv[0], " --lib: the library type;"
        print "Usage :", sys.argv[0], " --read: The sequencing strategy of producing reads with choices of P/S;"
        print "Usage :", sys.argv[0], ' --output: intron_id, gene_id,strand,chr,start,end,5SS inclusion counts,5SS skipping counts, 3SS includion counts,3SS skipping counts,skipping counts,intron counts;'
        print "Usage :", sys.argv[0], " --Total: the file store the total uniquely mapped reads."
        print datetime.datetime.now()
        sys.exit()
    if not os.path.exists(bam + '.bai'):
        if not bam.endswith('.sam'):
            cmd = 'samtools index %s' % bam
            os.system(cmd)
    parse_bam_file(Total, output, gtf, bam, read, lib, length, anchor, bin_size=1000)

if __name__ == "__main__":
    parse_args()
