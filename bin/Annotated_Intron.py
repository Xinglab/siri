#!/bin/python
import copy, getopt,re,os,sys,logging,time,datetime;

options, args = getopt.getopt(sys.argv[1:],'', ['gtf=','path='])
gtf = ''
path = ''
for opt, arg in options:
        if opt in('--gtf'):
                gtf = arg
        elif opt in('--path'):
                path = arg
if (not gtf or not path):
        print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "          A python program to get the annotated intron retention events from the given intron gtf file."
        print "Usage :", sys.argv[0], " --gtf: The gtf file, intron and exon gtf file will be located according to the name;"
        print "Usage :", sys.argv[0], " --path: The directory for the input and output gtf file."
        print datetime.datetime.now()
        sys.exit()

fr = open ("%s/Intron_%s" % (path,gtf))
fw = open(("%s/Intron_Annotated_%s" %(path,gtf)),"w")



pos = dict()
bin = 1000
intron = dict()
for info1 in fr:
	a1 = info1.strip().split("\t")
	gene_id = re.sub('.*gene_id "|\".*','',a1[8])
	exon1 = int(re.sub('.*exon_1 "|\".*','',a1[8]).split("_")[0])
        exon2 = int(re.sub('.*exon_2 "|\".*','',a1[8]).split("_")[1])
	key = "%s_%s_%s_%s_%s" % (a1[0],a1[3],a1[4],gene_id,a1[6])
	intron[key]=["false",a1[0],a1[3],a1[4],gene_id,a1[6]]
	index_s =  exon1/ bin
        index_e =  exon2/ bin
        for i in  range(index_s , index_e+1):
                pos.setdefault((a1[0],a1[6],i),[]).append(key)
		
fr.close()

fr2 = open ("%s/Exon_%s" % (path,gtf))
for info2 in fr2:
        a2 = info2.strip().split("\t")
	gene_id = re.sub('.*gene_id "|\".*','',a2[8])
        index_s =  int(a2[3]) / bin
        index_e =  int(a2[4]) / bin
	for i in  range(index_s , index_e+1):
                if (a2[0], a2[6], i) in pos:
                                for j in pos[a2[0],a2[6], i]:
					if(intron[j][4]== gene_id):
						if (int(a2[3]) <= int(intron[j][2]) and int(a2[4]) >= int(intron[j][3])):
                                                	intron[j][0] = "true"
fr2.close()

fr = open ("%s/Intron_%s" % (path,gtf))
for info1 in fr:
	a1 = info1.strip().split("\t")
	gene_id = re.sub('.*gene_id "|\".*','',a1[8])
	exon1 = int(re.sub('.*exon_1 "|\".*','',a1[8]).split("_")[0])
        exon2 = int(re.sub('.*exon_2 "|\".*','',a1[8]).split("_")[1])
	key = "%s_%s_%s_%s_%s" % (a1[0],a1[3],a1[4],gene_id,a1[6])
	fw.write('%s annotated_IR "%s";\n' %(info1.strip(),intron[key][0]))
fr.close()
fw.close()

