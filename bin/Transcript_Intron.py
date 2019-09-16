#!/bin/python
import copy, getopt,re,os,sys,logging,time,datetime;
options, args = getopt.getopt(sys.argv[1:],'', ['gtf=','path=','strand='])
gtf = ''
path=''
strand = 'unstrand'
for opt, arg in options:
        if opt in('--gtf'):
                gtf = arg
        elif opt in('--path'):
                path = arg
	elif opt in('--strand'):
                strand = arg
if (not gtf or not path):
        print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "          A python program to get the clean intron  from given intron gtf file."
        print "Usage :", sys.argv[0], " --gtf: The gtf file;"
        print "Usage :", sys.argv[0], " --path: The directory to the input and output gtf file;"
        print "Usage :", sys.argv[0], " --strand: The library type."
	print datetime.datetime.now()
        sys.exit()

trans = dict()
fr = open("%s/%s" % (path,gtf))
for info in fr:
	a = info.strip().split("\t")
	if(len(a) < 9):
                continue 
	if(a[2]!="exon"):
                continue
	transcript_id = re.sub('.*transcript_id "|\".*','',a[8])
	t_strand = a[6]
	if(strand == "unstrand"):
                t_strand = "*"
	if(trans.has_key(transcript_id)):
		if(int(a[3]) < trans[transcript_id][0]):
			trans[transcript_id][0] = int(a[3])
		if(int(a[4]) > trans[transcript_id][1]):
                        trans[transcript_id][1] = int(a[4])
	else:	
		trans[transcript_id]= [int(a[3]), int(a[4]),a[0], t_strand]
fr.close()

pos = dict()
bin = 1000
for t in trans:
	index_s = trans[t][0]/bin
	index_e = trans[t][1]/bin
	for i in  range(index_s , index_e+1):
                        pos.setdefault((trans[t][2],trans[t][3],i),[]).append(t)


fr2 = open ("%s/Intron_%s" % (path,gtf))
intron = dict()
for info2 in fr2:
	a2 = info2.strip().split("\t")
	key = "%s_%s_%s" % (a2[0],a2[3],a2[4])
        index_s =  int(a2[3]) / bin 
        index_e =  int(a2[4]) / bin 
        intron_strand = a2[6]
        if(strand == "unstrand"):
                intron_strand = "*" 
        for i in  range(index_s , index_e+1):
                if (a2[0], intron_strand, i) in pos:
                                for j in pos[a2[0],intron_strand, i]: 
                                        if(int(a2[3]) >= trans[j][0] and int(a2[4])<= trans[j][1]):
                                                intron.setdefault(key,[]).append(j)

fr2.close()




fw = open(("%s/Intron_transcript.txt" %(path)),"w")
for i in intron:
	fw.write('%s\t%s\n' %(i, "\t".join(list(set(intron[i])))))
fw.close()















