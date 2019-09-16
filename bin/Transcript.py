#!/usr/bin/python
import getopt, re, os, sys, logging, time, datetime, copy

options, args = getopt.getopt(sys.argv[1:], '', ['FPKM=', 'intron=', 'output='])
FPKM = ''
intron = ''
output = ''
for opt, arg in options:
    if opt in ('--FPKM'):
        FPKM = arg
    elif opt in ('--intron'):
        intron = arg
    elif opt in ('--output'):
        output = arg
if not FPKM or not intron or not output:
    print "not enough parameters"
    print "Program : ", sys.argv[0]
    print "a python program to get the FPKM value from give intron annotation file and FPKM of transcript"
    print "usage :", sys.argv[0], " --FPKM: transcript expression file"
    print "usage :", sys.argv[0], " --intron: intron file with overlapping transcript"
    print "usage :", sys.argv[0], " --output: the output file"
    print datetime.datetime.now()
    sys.exit()
exp = dict()
sample_number = 0
fr = open(FPKM)
for info in fr:
    sp = info.strip().split("\t")
    if sample_number == 0:
        sample_number = len(sp) - 1
    exp[sp[0]] = sp[1:]
fr.close()
if sample_number == 0:
    print 'Empty FPKM files'
    sys.exit()
exp_intron = dict()
fr = open(intron)
for info in fr:
    a = info.strip().split("\t")
    t = [0] * sample_number
    for i in a[1:]:
        if (not exp.has_key(i)):
            continue
        for j in range(0, sample_number):
            t[j] += float(exp[i][j])

    exp_intron[a[0]] = t
fr.close()

fw = open(output, "w")
fr = open(intron)
for info in fr:
    a = info.strip().split("\t")
    fw.write("%s\t%s\n" % (a[0], "\t".join(str(x) for x in exp_intron[a[0]])))
fr.close()
fw.close()
