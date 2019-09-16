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

fr1 = open ("%s/Intron_%s" % (path,gtf))
pos = dict()
bin = 1000
intron = dict()
for info1 in fr1:
	a1 = info1.strip().split("\t")
	key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
	gene_id = re.sub('.*gene_id "|\".*','',a1[8])
	if(intron.has_key(key)):
                gg = intron[key][5].split(",")
                gg_l = "true"
                for g_id in gg:
                        if(gene_id == g_id):
                                gg_l = "false"
                if(gg_l =="true"):
                        intron[key][5] += ","
                        intron[key][5] += gene_id
        else:
		intron[key]=["true","true",a1[0],int(a1[3]),int(a1[4]),gene_id]
	index_s =  int(a1[3])/ bin
        index_e =  int(a1[4])/ bin
	if(strand=="unstrand"):
		for i in  range(index_s , index_e+1):
                        pos.setdefault((a1[0],"*",i),[]).append(key)
	else:
        	for i in  range(index_s , index_e+1):
                	pos.setdefault((a1[0],a1[6],i),[]).append(key)
		
fr1.close()

fr2 = open ("%s/Exon_%s" % (path,gtf))
for info2 in fr2:
        a2 = info2.strip().split("\t")
	gene_id = re.sub('.*gene_id "|\".*','',a2[8])
        index_s =  int(a2[3]) / bin
        index_e =  int(a2[4]) / bin
	exon_strand = a2[6]
	if(strand == "unstrand"):
		exon_strand = "*"
	for i in  range(index_s , index_e+1):
		if (a2[0], exon_strand, i) in pos:
				for j in pos[a2[0],exon_strand, i]:
					#if(int(a2[3]) <= intron[j][3] and int(a2[4]) >= intron[j][3]):
					if(int(a2[3]) >= intron[j][3] and int(a2[3]) < intron[j][4]):
						intron[j][0]="false"
					if( (int(a2[4]) > intron[j][3]) and  (int(a2[4]) <= intron[j][4])):
					#if( (int(a2[3]) > intron[j][3]) and  (int(a2[3]) <= intron[j][4])):
						intron[j][0]="false"
                                                

fr2.close()

fr2 = open ("%s/Intron_%s" % (path,gtf))
for info2 in fr2:
        a2 = info2.strip().split("\t")
        gene_id = re.sub('.*gene_id "|\".*','',a2[8])
        index_s =  int(a2[3]) / bin 
        index_e =  int(a2[4]) / bin 
        intron_strand = a2[6]
        if(strand == "unstrand"):
                intron_strand = "*" 
        for i in  range(index_s , index_e+1):
                if (a2[0], intron_strand, i) in pos:
                                for j in pos[a2[0],intron_strand, i]: 
					if(int(a2[3]) == intron[j][3] and int(a2[4])== intron[j][4]):
						continue
					if(int(a2[3]) == intron[j][3] or int(a2[4])==intron[j][4]):
                                               intron[j][1]="false"
					if(int(a2[3]) < intron[j][3] and int(a2[4]) > intron[j][3]):
                                                intron[j][1]="false"
                                        if( (int(a2[3]) > intron[j][3]) and  (int(a2[3]) < intron[j][4])):
                                                intron[j][1]="false"        
					
fr2.close()



fw = open(("%s/Intron_attri_%s" %(path,gtf)),"w")
fr = open ("%s/Intron_%s" % (path,gtf))
for info1 in fr:
	a1 = info1.strip().split("\t")
	key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
	fw.write('%s clean "%s"; clean_simple "%s";\n' %(info1.strip(),intron[key][0],intron[key][1]))
fr.close()
fw.close()















