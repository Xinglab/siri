#!/usr/bin/python
import getopt, re, os, sys, logging, time, datetime, copy
from collections import defaultdict

options, args = getopt.getopt(sys.argv[1:], '', ['gtf=', 'path=', 'strand='])
gtf, path, strand = '', '', 'unstrand'
for opt, arg in options:
    if opt in '--gtf':
        gtf = arg
    elif opt in '--path':
        path = arg
    elif opt in '--strand':
        strand = arg
if not gtf or not path:
    print "Not enough parameters!"
    print "Program : ", sys.argv[0]
    print "          A python program to get the clean intron  from given intron gtf file."
    print "Usage :", sys.argv[0], " --gtf: The gtf file;"
    print "Usage :", sys.argv[0], " --path: The directory to the input and output gtf file;"
    print "Usage :", sys.argv[0], " --strand: The library type."
    print datetime.datetime.now()
    sys.exit()


def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List


def intron_gtf(gtf_file_list):
    exon_s, exon_e = defaultdict(lambda: []), defaultdict(lambda: [])
    parse_line_start = 0
    for i in xrange(0, len(gtf_file_list)):
        if gtf_file_list[i].startswith('#') or gtf_file_list[i].split('\t')[2] != 'exon':
            continue
        parse_line_start = i
        break
    for i in xrange(parse_line_start, len(gtf_file_list)):
        gtf_list = gtf_file_list[i].strip().split('\t')
        if gtf_list[2] != 'exon':
            continue
        transcript_id = re.sub('.*transcript_id "|\".*', '', gtf_list[8])
        exon_s[transcript_id].append(int(gtf_list[3]))
        exon_e[transcript_id].append(int(gtf_list[4]))
    trans = {}
    fw = open(("%s/Intron_%s" % (path, gtf)), "w")
    for i in xrange(parse_line_start, len(gtf_file_list)):
        gtf_list = gtf_file_list[i].strip().split('\t')
        transcript_id = re.sub('.*transcript_id "|\".*', '', gtf_list[8])
        gene_id = re.sub('.*gene_id "|\".*', '', gtf_list[8])
        gene_name = re.sub('.*gene_name "|\".*', '', gtf_list[8])
        if transcript_id in trans:
            trans[transcript_id] = 'false'
        elif len(exon_s[transcript_id]) > 1:
            exon_s[transcript_id].sort()
            exon_e[transcript_id].sort()
            for i in xrange(1, len(exon_s[transcript_id])):
                if gtf_list[6] == '+':
                    fw.write(
                        "%s\t%s\tintron\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; transcript_id \"%s\"; intron_number \"%s\"; total_intron_number \"%s\"; exon_1 \"%s_%s\"; exon_2 \"%s_%s\";\n" % (
                            gtf_list[0], gtf_list[1], str(int(exon_e[transcript_id][i - 1]) + 1),
                            str(int(exon_s[transcript_id][i]) - 1),
                            gtf_list[6], gene_id, gene_name, transcript_id, i, len(exon_s[transcript_id]) - 1,
                            exon_s[transcript_id][i - 1], exon_e[transcript_id][i - 1], exon_s[transcript_id][i],
                            exon_e[transcript_id][i]))
                else:
                    fw.write(
                        "%s\t%s\tintron\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; transcript_id \"%s\"; intron_number \"%s\"; total_intron_number \"%s\"; exon_1 \"%s_%s\"; exon_2 \"%s_%s\";\n" % (
                            gtf_list[0], gtf_list[1], str(int(exon_e[transcript_id][i - 1]) + 1),
                            str(int(exon_s[transcript_id][i]) - 1),
                            gtf_list[6], gene_id, gene_name, transcript_id, len(exon_s[transcript_id]) - i,
                            len(exon_s[transcript_id]) - 1, exon_s[transcript_id][i - 1], exon_e[transcript_id][i - 1],
                            exon_s[transcript_id][i], exon_e[transcript_id][i]))
        trans[transcript_id] = 'true'
    fw.close()

def overlap(min1, max1, min2, max2):
    return max(0, min(max1 + 1, max2 + 1) - max(min1, min2))

def transcript_intron(gtf_file_list):
    dict_intron2transcript = defaultdict(lambda:[])
    for gtf_line in gtf_file_list:
        gtf_list = gtf_line.strip().split('\t')
        transcript_id = re.sub('.*transcript_id "|\".*', '', gtf_list[8])
        if strand == 'unstrand':
            key = '{}_{}_{}_{}'.format(gtf_list[0], gtf_list[3], gtf_list[4], '*')
        else:
            key = '{}_{}_{}_{}'.format(gtf_list[0], gtf_list[3], gtf_list[4], gtf_list[6])
        dict_intron2transcript[key].append(transcript_id)
    fw = open("%s/Intron_transcript.txt" % (path), "w")
    for intron_id in dict_intron2transcript:
        fw.write('%s\t%s\n' % (intron_id, "\t".join(list(set(dict_intron2transcript[intron_id])))))
    fw.close()

def annotated_intron(gtf_file_list):
    fw = open("%s/Intron_Annotated_%s" % (path, gtf), "w")
    pos = defaultdict(lambda: [])
    bin = 1000
    intron = defaultdict(lambda: [])
    for gtf_line in gtf_file_list:
        gtf_list = gtf_line.strip().split('\t')
        gene_id = re.sub('.*gene_id "|\".*', '', gtf_list[8])
        exon1 = int(re.sub('.*exon_1 "|\".*', '', gtf_list[8]).split("_")[0])
        exon2 = int(re.sub('.*exon_2 "|\".*', '', gtf_list[8]).split("_")[1])
        key = "%s_%s_%s_%s_%s" % (gtf_list[0], gtf_list[3], gtf_list[4], gene_id, gtf_list[6])
        intron[key] = ["false", gtf_list[0], gtf_list[3], gtf_list[4], gene_id, gtf_list[6]]
        index_s = exon1 / bin
        index_e = exon2 / bin
        for i in xrange(index_s, index_e + 1):
            pos[(gtf_list[0], gtf_list[6], i)].append(key)
    exon_file = open("%s/Exon_%s" % (path, gtf))
    for line in exon_file:
        sp = line.strip().split('\t')
        gene_id = re.sub('.*gene_id "|\".*', '', sp[8])
        index_s = int(sp[3]) / bin
        index_e = int(sp[4]) / bin
        for i in range(index_s, index_e + 1):
            if (sp[0], sp[6], i) in pos:
                for j in pos[sp[0], sp[6], i]:
                    if intron[j][4] == gene_id:
                        if int(sp[3]) <= int(intron[j][2]) and int(sp[4]) >= int(intron[j][3]):
                            intron[j][0] = "true"
    exon_file.close()
    intron_file = open("%s/Intron_%s" % (path, gtf))
    for line in intron_file:
        sp = line.strip().split('\t')
        gene_id = re.sub('.*gene_id "|\".*', '', sp[8])
        key = "%s_%s_%s_%s_%s" % (sp[0], sp[3], sp[4], gene_id, sp[6])
        fw.write('%s annotated_IR "%s";\n' % (line.strip(), intron[key][0]))
    fw.close()


def attribute_intron(gtf_file_list):
    pos = defaultdict(lambda: [])
    bin = 1000
    intron = defaultdict(lambda: [])
    for line in gtf_file_list:
        gtf_list = line.strip().split('\t')
        if strand == 'unstrand':
            key = "%s_%s_%s_%s" % (gtf_list[0], gtf_list[3], gtf_list[4], '*')
        else:
            key = "%s_%s_%s_%s" % (gtf_list[0], gtf_list[3], gtf_list[4], gtf_list[6])
        gene_id = re.sub('.*gene_id "|\".*', '', gtf_list[8])
        if key in intron:
            gene_id_list = intron[key][5].split(',')
            id_duplicate = 'true'
            for _gene_id in gene_id_list:
                if _gene_id == gene_id:
                    id_duplicate = 'false'
            if id_duplicate == 'true':
                intron[key][5] += ','
                intron[key][5] += gene_id
        else:
            intron[key] = ['true', 'true', gtf_list[0], int(gtf_list[3]), int(gtf_list[4]), gene_id]
        index_s = int(gtf_list[3]) / bin
        index_e = int(gtf_list[4]) / bin
        if strand == 'unstrand':
            for i in xrange(index_s, index_e + 1):
                pos[(gtf_list[0], '*', i)].append(key)
        else:
            for i in xrange(index_s, index_e + 1):
                pos[(gtf_list[0], gtf_list[6], i)].append(key)

    exon_file = open("%s/Exon_%s" % (path, gtf))
    for line in exon_file:
        sp = line.strip().split("\t")
        index_s = int(sp[3]) / bin
        index_e = int(sp[4]) / bin
        exon_strand = sp[6]
        if strand == "unstrand":
            exon_strand = "*"
        for i in range(index_s, index_e + 1):
            if (sp[0], exon_strand, i) in pos:
                for j in pos[sp[0], exon_strand, i]:
                    if int(sp[3]) < int(intron[j][3]) and int(sp[4]) > int(intron[j][4]):
                        continue
                    if overlap(int(sp[3]), int(sp[4]), intron[j][3], intron[j][4]) > 0:
                        intron[j][0] = 'false'
    exon_file.close()

    for line in gtf_file_list:
        sp = line.strip().split('\t')
        index_s = int(sp[3]) / bin
        index_e = int(sp[4]) / bin
        intron_strand = sp[6]
        if (strand == "unstrand"):
            intron_strand = "*"
        for i in range(index_s, index_e + 1):
            if (sp[0], intron_strand, i) in pos:
                for j in pos[sp[0], intron_strand, i]:
                    if (int(sp[3]) == intron[j][3] and int(sp[4]) == intron[j][4]):
                        continue
                    if overlap(int(sp[3]), int(sp[4]), intron[j][3], intron[j][4]) > 0:
                        intron[j][1] = "false"

    fw = open(("%s/Intron_attri_%s" % (path, gtf)), "w")
    for line in gtf_file_list:
        sp = line.strip().split('\t')
        if strand == 'unstrand':
            key = "%s_%s_%s_%s" % (sp[0], sp[3], sp[4], '*')
        else:
            key = "%s_%s_%s_%s" % (sp[0], sp[3], sp[4], sp[6])
        fw.write('%s clean "%s"; clean_simple "%s";\n' % (line.strip(), intron[key][0], intron[key][1]))
    fw.close()


if __name__ == '__main__':
    gtf_file_list = read_file('%s/%s' % (path, gtf))
    intron_gtf(gtf_file_list)
    gtf_file_list = read_file("%s/Intron_%s" % (path, gtf))
    transcript_intron(gtf_file_list)
    annotated_intron(gtf_file_list)
    attribute_intron(gtf_file_list)
