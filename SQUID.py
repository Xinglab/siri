#!/bin/python
import getopt, copy, re, os, sys, logging, time, datetime
from multiprocessing import Pool

options, args = getopt.getopt(sys.argv[1:], 'o:',
                              ['align=', 'GTF=', 'fastq=', 'check_len=', 'index_star=', 'index_kallisto=', 'l=', 's=',
                               'output=', 'update=', 'lib=', 'read=', 'length=', 'anchor=', 'Cal=', 'FPKM=',
                               'Comparison=', 'analysis=', 'c1=', 'p=', 'resume=', 'F_deltaPI=', 'F_FDR='])
## initial parameters
align = ''
GTF = ''
fastq = ''
check_len = 'false'
index_star = ''
index_kallisto = ''
l = 200
s = 100
output = '.'
update = 'false'
lib = 'unstrand'
read = 'P'
length = 100
anchor = 8
Cal = 'All'
FPKM = ''
Comparison = ''
analysis = 'U'
c1 = 0.0001
p = 1
F_deltaPI = 0.05
F_FDR = 0.1
resume = "false"
for opt, arg in options:
    if opt in ('--align'):
        align = arg
    elif opt in ('--GTF'):
        GTF = arg
    elif opt in ('--fastq'):
        fastq = arg
    elif opt in ('--check_len'):
        check_len = arg
    elif opt in ('--index_star'):
        index_star = arg
    elif opt in ('--index_kallisto'):
        index_kallisto = arg
    elif opt in ('--l'):
        l = int(arg)
    elif opt in ('--s'):
        s = int(arg)
    elif opt in ('-o', '--output'):
        output = arg
    elif opt in ('--update'):
        update = arg
    elif opt in ('--lib'):
        lib = arg
    elif opt in ('--read'):
        read = arg
    elif opt in ('--length'):
        length = int(arg)
    elif opt in ('--anchor'):
        anchor = int(arg)
    elif opt in ('--Cal'):
        Cal = arg
    elif opt in ('--FPKM'):
        FPKM = arg
    elif opt in ('--Comparison'):
        Comparison = arg
    elif opt in ('--analysis'):
        analysis = arg
    elif opt in ('--c1'):
        c1 = float(arg)
    elif opt in ('--p'):
        p = int(arg)
    elif opt in ('--F_deltaPI'):
        F_deltaPI = float(arg)
    elif opt in ('--F_FDR'):
        F_FDR = float(arg)
    elif opt in ('--resume'):
        resume = arg
run = "true"
if not align:
    if not fastq or not index_star:
        print "Please provide either alignment file or fastq and star index to do the alignment"
        run = "false"
if not FPKM and not align:
    if not fastq or not index_kallisto:
        print "please provide either FPKM file or fastq and kallisto index to generate FPKM file"
        run = "false"
if not GTF:
    print "please provide gtf annotation file"
    run = "false"
if run == "false":
    print "Not enough parameters!"
    print "Program : ", sys.argv[0]
    print "          A python program to calculate the retained intron level and differential retained introns.\n"
    print "Usage :", sys.argv[
        0], " --align: s1.bam/s1.sam[,s2.bam/s2.sam,...]. Mapping results for all of samples in bam/sam format. Different samples are sepreated by commas;"
    print "Usage :", sys.argv[0], " --GTF: The gtf file;"
    print "Usage :", sys.argv[
        0], " --fastq: s1_1.fq[:s1_2.fq][,s1_1.fq[:s2_2.fq],...]. The raw sequencing reads in fasta or fastq format that is required to call kallisto to calculate FPKM values;"
    print "Usage :", sys.argv[
        0], " --check_len: Whether to generate new fastq files to with equal read length. The default value is false;"
    print "Usage :", sys.argv[
        0], " --index_star: The path to the star index that is required to do the alignment using STAR;"
    print "Usage :", sys.argv[
        0], " --index_kallisto: The path to the kallisto index that is required to run kallisto from raw reads;"
    print "Usage :", sys.argv[
        0], " --l: Estimated average fragment length. The parameter to run kallisto with default value of 200;"
    print "Usage :", sys.argv[
        0], " --s: Estimated standard deviation of fragment length. The parameter to run kallisto with default value of 100;"
    print "Usage :", sys.argv[0], " -o/--output: The output directory. The default is current directory;"
    print "Usage :", sys.argv[
        0], " --update: Whether to update the attributes of introns using spliced reads. The default is false;"
    print "Usage :", sys.argv[
        0], " --lib: The library type with choices of unstrand/first/second. The details are explained in the parameter of library-type in tophat2. The default is unstrand;"
    print "Usage :", sys.argv[
        0], " --read: The sequencing strategy of producing reads with choices of (paired end) or S (single end). The default is P;"
    print "Usage :", sys.argv[0], " --length: The read length of sequencing reads. The default length is 100;"
    print "Usage :", sys.argv[
        0], " --anchor: The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8;"
    print "Usage :", sys.argv[
        0], " --Cal: Which  part of the program user choose to run, the choices are All/count/DSI. All means run the whole program, count means only run the PI value calculation part, DSI means only run the differential analysis of spliced introns. The default is All;"
    print "Usage :", sys.argv[
        0], " --FPKM: A file providing the FPKM value for each sample, the first column is transcript ID with the following column being the FPKM value for each sample. If it is not provided, kallisto will be called to calculate FPKM value;"
    print "Usage :", sys.argv[
        0], " --Comparison: A file providing the sample pairs to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the align files replicates seperated by commas), column 3 (sample 2 order in the align files replicates seperated by commas), column 4 (optional, if present as 'pool', the replicates are combined together in rMATS calculation). If absent, the step of calculation of differential spliced introns  will be skipped;"
    print "uasge: ", sys.argv[
        0], " --analysis: Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U;"
    print "Usage :", sys.argv[
        0], " --c1: The cutoff of splicing difference using Junction method. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001;"
    print "Usage :", sys.argv[0], " --p: The number of threads used to run rMATS. The default is 1;"
    print "Usage :", sys.argv[
        0], " --F_deltaPI: The cutoff of delta to output differential spliced introns.The default is 0.05;"
    print "Usage :", sys.argv[
        0], " --F_FDR: The cutoff of combined FDR to output differential spliced introns.The default is 0.05;"
    print "Usage :", sys.argv[0], " --resume: Whether to resume previous run. The default is false;"
    print datetime.datetime.now()
    sys.exit()


def listToString(str_list):
    return ' '.join(str_list)


if not os.path.exists(output):
    os.system("mkdir %s" % output)

# setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=output + '/log.SQUID' + str(datetime.datetime.now()) + '.txt',
                    filemode='w')

# getting Start Time
if (resume == "true"):
    logging.debug('Resume the program with [%s]\n', listToString(sys.argv))
else:
    logging.debug('Start the program with [%s]\n', listToString(sys.argv))
startTime = time.time()

# get path of the main program
path = os.path.abspath(os.path.dirname(__file__))

# get the path of the python programs

bin_path = "%s/bin" % path
samples = align.split(",")
fq = fastq.split(",")
ALIGN = align.split(",")
if (not align):
    samples = fastq.split(",")
    ALIGN = []
sample_number = len(samples)

# get the path  of count files
count_path = "%s/counts" % output
if (not os.path.exists(count_path)):
    os.system("mkdir %s" % count_path)

# get the path of the results
output_path = "%s/Result" % output
if (not os.path.exists(output_path)):
    os.system("mkdir %s" % output_path)
normF = [0] * sample_number

if Cal == "All" or Cal == "count":
    if Cal == "All":
        logging.debug("Run the whole program\n")
    else:
        logging.debug("Run the PI value calculation parts only\n")
    # make directory for gtf files
    gtf_path = "%s/gtf_files" % output
    if not os.path.exists(gtf_path):
        os.system("mkdir %s" % gtf_path)
    gtf = re.sub(".*/", "", GTF)

    attri = "%s/Intron_attri_%s" % (gtf_path, gtf)
    if (resume == "false" or not os.path.exists(attri)):
        resume = "false"
        cmd = "cp %s %s/%s" % (GTF, gtf_path, gtf)
        os.system(cmd)

        # use the awk command to generate Exon.gtf file
        logging.debug("################ Generating  the intron gtf files #######################\n")
        cmd = "less %s/%s | awk '{if($3==\"exon\"){print $0}}'> %s/Exon_%s" % (gtf_path, gtf, gtf_path, gtf)
        os.system(cmd)
        logging.debug("gtf_files\Exon_" + gtf)

        # use the Intron_gtf.py to generate Intron.gtf file
        cmd = "python %s/Intron_gtf.py --gtf %s --path %s --strand %s" % (bin_path, gtf, gtf_path, lib)
        os.system(cmd)
        logging.debug("gtf_files\Intron_" + gtf)
        logging.debug("gtf_files\Intron_transcript.txt")
        logging.debug("gtf_files\Intron_Annotated_" + gtf)
        logging.debug("gtf_files\Intron_attri_" + gtf)
        logging.debug("#########################################################################\n")

    if fq and check_len != "false":
        fq_path = "%s/fq" % output
        if not os.path.exists(fq_path):
            os.system("mkdir %s" % fq_path)
        start = 0
        if (resume == "true"):
            for ss in range(0, sample_number):
                fq1 = "%s/sample_%s.fq" % (fq_path, ss)
                fq2 = "%s/sample_%s_2.fq" % (fq_path, ss)
                if (os.path.exists(fq1) or os.path.exists(fq2)):
                    start += 1
                else:
                    break
        for ss in range(start, sample_number):
            format = "cat "
            if (re.sub(".*\\.", "", fq[ss]) == "gz"):
                format = "zcat "
            if (read == "S"):
                cmd = format + fq[
                    ss] + ' | awk \'{a=$0; getline b;if(length(b)==' + str(
                    length) + '){print a"\\n"b}}\' > ' + fq_path + '/sample_' + str(
                    ss) + '.fq'
                logging.debug(cmd)
                os.system(cmd)
                fq[ss] = "%s/sample_%s.fq" % (fq_path, ss)
            else:

                cmd = format + fq[ss].split(":")[
                    0] + ' | awk \'{a=$0; getline b;if(length(b)==' + str(
                    length) + '){print a"\\n"b}}\' > ' + fq_path + '/sample_' + str(
                    ss) + '_1.fq'
                logging.debug(cmd)
                os.system(cmd)
                cmd = format + fq[ss].split(":")[
                    1] + ' | awk \'{a=$0; getline b;if(length(b)==' + str(
                    length) + '){print a"\\n"b}}\' > ' + fq_path + '/sample_' + str(
                    ss) + '_2.fq'
                logging.debug(cmd)
                os.system(cmd)
                fq[ss] = "%s/sample_%s_1.fq:%s/sample_%s_2.fq:" % (fq_path, ss, fq_path, ss)

    # check whether the alignment file exist:
    if not align:
        align_path = "%s/align" % output
        if (not os.path.exists(align_path)):
            os.system("mkdir %s" % align_path)

        start = -1
        for ss in range(0, sample_number):
            align_file = "%s/star_%s" % (align_path, ss)
            if (os.path.exists(align_file)):
                start += 1
            else:
                break
        if (start == -1 or resume == "false"):
            start = 0
        for ss in range(start, sample_number):
            if (re.sub(".*\\.", "", fq[ss]) == "gz"):
                cmd = "STAR --genomeDir %s --readFilesIn %s --readFilesCommand zcat --runThreadN 4 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --outFilterMismatchNmax 8 --alignEndsType EndToEnd --outFileNamePrefix %s/S%s --outSAMattributes All" % (
                    index_star, re.sub(":", " ", fq[ss]), align_path, ss)
            else:
                cmd = "STAR --genomeDir %s --readFilesIn %s --runThreadN 4 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --outFilterMismatchNmax 8 --alignEndsType EndToEnd --outFileNamePrefix %s/S%s --outSAMattributes All" % (
                    index_star, re.sub(":", " ", fq[ss]), align_path, ss)

            logging.debug(cmd)
            os.system(cmd)
            ALIGN.append("%s/S%sAligned.sortedByCoord.out.bam" % (align_path, ss))
            print ALIGN[-1]

    # generate the counts files
    COUNT_file = "%s/Total.txt" % count_path
    if resume == "false" or not os.path.exists(COUNT_file):
        start = 0
        for ss in range(0, sample_number):
            exon = "%s/count_%s_exon.txt" % (count_path, ss)
            if (os.path.exists(exon)):
                start += 1
            else:
                break
        run_cmd_list = []
        for ss in range(start, sample_number):
            print ss, ALIGN[ss]
            if ALIGN[ss].endswith('.bam') or ALIGN[ss].endswith('.sam'):
                cmd = "python %s/bam_count.py --gtf %s/Intron_%s,%s/%s --length %s --anchor %s --bam %s -o %s/count_%s --lib %s --read %s --Total %s/Total_%s.txt" % (
                    bin_path, gtf_path, gtf, gtf_path, gtf, length, anchor, ALIGN[ss], count_path, ss, lib, read,
                    count_path, ss)
                run_cmd_list.append(cmd)
                logging.debug(cmd)
        if len(run_cmd_list) > 0:
            pool = Pool(processes=int(p))
            pool.map(os.system, run_cmd_list)
            pool.close()
            pool.join()
        logging.debug("Generate the count file")

        intron_header = dict()
        for ss in range(0, sample_number):
            FR = open("%s/count_%s_intron.txt" % (count_path, ss))
            for info in FR:
                a = info.strip().split("\t")
                if (intron_header.has_key(a[0])):
                    t = a[3].split(",")
                    if (t[0] == "false"):
                        intron_header[a[0]][2] = re.sub("true,", "false,", intron_header[a[0]][2])
                    if (t[1] == "false"):
                        intron_header[a[0]][2] = re.sub(",true", ",false", intron_header[a[0]][2])
                else:
                    intron_header[a[0]] = a[1:7]
            FR.close()

        fw = open("%s/intron_header" % count_path, "w")
        FR = open("%s/count_0_intron.txt" % count_path)

        for info in FR:
            a = info.strip().split("\t")
            fw.write("%s\t%s\n" % (a[0], "\t".join(intron_header[a[0]])))
        fw.close()
        FR.close()

        cmd = "cut -f 1 %s/count_0_exon.txt > %s/exon_header" % (count_path, count_path)
        logging.debug(cmd)
        os.system(cmd)

        cmd = "cut -f 3,4 %s/count_0_exon.txt > %s/exon_tail" % (count_path, count_path)
        logging.debug(cmd)
        os.system(cmd)

        cmd1 = "paste "
        cmd2 = "paste %s/intron_header " % count_path
        cmd3 = "paste %s/exon_header " % count_path
        for ss in range(0, sample_number):
            cmd = "cut -f 8,9,10,11,12,13 %s/count_%s_intron.txt > %s/intron_%s.val" % (count_path, ss, count_path, ss)
            logging.debug(cmd)
            os.system(cmd)
            cmd = "cut -f 2 %s/count_%s_exon.txt > %s/exon_%s.val" % (count_path, ss, count_path, ss)
            logging.debug(cmd)
            os.system(cmd)
            cmd1 += "%s/Total_%s.txt " % (count_path, ss)
            cmd2 += "%s/intron_%s.val " % (count_path, ss)
            cmd3 += "%s/exon_%s.val " % (count_path, ss)

        cmd1 += " > %s/Total.txt" % (count_path)
        cmd2 += " > %s/count_intron.txt" % (count_path)
        cmd3 += " %s/exon_tail  > %s/count_exon.txt" % (count_path, count_path)
        logging.debug(cmd1)
        os.system(cmd1)
        logging.debug(cmd2)
        os.system(cmd2)
        logging.debug(cmd3)
        os.system(cmd3)
        logging.debug("Finish generating the count files\n")
        logging.debug("#########################################################################\n")

    fr = open("%s/Total.txt" % count_path)
    info = fr.readline()
    fr.close()
    normF = map(int, info.strip().split("\t"))

    # generate the PI_Density counts
    if (not os.path.exists(FPKM) or FPKM == ''):
        FPKM_path = "%s/FPKM" % output
        trans_FPKM = dict()
        FPKM = "%s/transcript_exp.txt" % FPKM_path
        if (not os.path.exists(FPKM_path)):
            os.system("mkdir %s" % FPKM_path)
        if (resume == "false" or not os.path.exists(FPKM)):
            if (fastq == '' or index_kallisto == ''):
                l_type = "fr-unstranded"
                if (lib == "first"):
                    l_type = "fr-firststrand"
                if (lib == "second"):
                    l_type = "fr-secondstrand"
                start = -1
                for ss in range(0, sample_number):
                    cufflinks_file = "%s/cufflinks_%s" % (FPKM_path, ss)
                    if (os.path.exists(cufflinks_file)):
                        start += 1
                    else:
                        break

                if (start == -1 or resume == "false"):
                    start = 0
                for ss in range(start, sample_number):
                    cmd = "cufflinks --GTF %s/%s -p %s --library-type %s --multi-read-correct -o %s/cufflinks_%s %s" % (
                        gtf_path, gtf, p, l_type, FPKM_path, ss, ALIGN[ss])
                    logging.debug(cmd)
                    os.system(cmd)
                for ss in range(0, sample_number):
                    fr = open("%s/cufflinks_%s/isoforms.fpkm_tracking" % (FPKM_path, ss))
                    info = fr.readline()
                    for info in fr:
                        a = info.strip().split("\t")
                        if (trans_FPKM.has_key(a[0])):
                            trans_FPKM[a[0]][ss] = a[9]
                        else:
                            trans_FPKM[a[0]] = [0] * sample_number
                            trans_FPKM[a[0]][ss] = a[9]
                    fr.close()

                fw = open(FPKM, "w")
                for rp in trans_FPKM:
                    fw.write("%s\t%s\n" % (rp, "\t".join(str(x) for x in trans_FPKM[rp])))
                fw.close()
            else:
                start = -1
                for ss in range(0, sample_number):
                    kallisto_file = "%s/kallisto_%s" % (FPKM_path, ss)
                    if (os.path.exists(kallisto_file)):
                        start += 1
                    else:
                        break

                if (start == -1 or resume == "false"):
                    start = 0
                for ss in range(start, sample_number):
                    if (read == "P"):
                        cmd = "kallisto quant --index=%s --output-dir=%s/kallisto_%s --threads=%s --plaintext %s" % (
                            index_kallisto, FPKM_path, ss, p, re.sub(":", " ", fq[ss]))
                    else:
                        cmd = "kallisto quant --index=%s --output-dir=%s/kallisto_%s --single -l %s -s %s --threads=%s --plaintext %s" % (
                            index_kallisto, FPKM_path, ss, l, s, p, re.sub(":", " ", fq[ss]))
                    logging.debug(cmd)
                    os.system(cmd)
                for ss in range(0, sample_number):
                    fr = open("%s/kallisto_%s/abundance.tsv" % (FPKM_path, ss))
                    info = fr.readline()
                    for info in fr:
                        a = info.strip().split("\t")
                        if (a[0] in trans_FPKM):
                            if (read == "P"):
                                trans_FPKM[a[0]][ss] = float(a[3]) * 2 / float(a[2]) * 1000 * 1000000 / normF[ss]
                            else:
                                trans_FPKM[a[0]][ss] = float(a[3]) / float(a[2]) * 1000 * 1000000 / normF[ss]
                        else:
                            trans_FPKM[a[0]] = [0] * sample_number
                            if (read == "P"):
                                trans_FPKM[a[0]][ss] = float(a[3]) * 2 / float(a[2]) * 1000 * 1000000 / normF[ss]
                            else:
                                trans_FPKM[a[0]][ss] = float(a[3]) / float(a[2]) * 1000 * 1000000 / normF[ss]
                    fr.close()
                fw = open(FPKM, "w")
                for rp in trans_FPKM:
                    fw.write("%s\t%s\n" % (rp, "\t".join(str(x) for x in trans_FPKM[rp])))
                fw.close()

    if (os.path.exists(FPKM)):
        Trans = dict()
        fr = open(FPKM)
        for info in fr:
            a = info.strip().split("\t")
            Trans[a[0]] = a[1:]
        fr.close()
        Intron = dict()

        fr = open("%s/Intron_transcript.txt" % gtf_path)
        for info in fr:
            a = info.strip().split("\t")
            Intron[a[0]] = a[1:]
        fr.close()

        intron_obs = dict()
        intron_exp = dict()
        fr = open("%s/counts/count_intron.txt" % output)
        for info in fr:
            a = info.strip().split("\t")
            ri_obs = []
            ri_exp = []
            Intron_l = int(a[6]) - int(a[5]) + 1
            gene_id = a[1].split(",")
            for i in range(0, sample_number):
                ri_obs.append(int(a[i * 6 + 12]) + int(a[i * 6 + 9]) + int(a[i * 6 + 7]))
                exp_expression = 0

                for e in Intron[a[0]]:
                    if (e in Trans):
                        exp_expression += float(Trans[e][i])
                ri_exp.append(int(exp_expression * Intron_l * normF[i] / (1000 * 1000000)))
            intron_obs[a[0]] = ri_obs
            intron_exp[a[0]] = ri_exp
        fr.close()
        fw1 = open("%s/counts/count_all_Density.txt" % output, "w")
        fr = open("%s/counts/count_intron.txt" % output)
        for info in fr:
            a = info.strip().split("\t")
            Intron_l = int(a[6]) - int(a[5]) + 1
            fw1.write("%s\t%s\t%s\t%s\n" % (
                a[0], Intron_l, "\t".join(str(x) for x in intron_obs[a[0]]),
                "\t".join(str(x) for x in intron_exp[a[0]])))
        fw1.close()
        fr.close()
        logging.debug("Generate the Density counts")

    # generate the attributes of the intron
    intron_anno = dict()
    fr = open("%s/Intron_Annotated_%s" % (gtf_path, gtf))
    for info in fr:
        a = info.strip().split("\t")
        ann = re.sub('.*annotated_IR "|\".*', '', a[8])
        if lib == 'unstrand':
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], '*')
        else:
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], a[6])
        if (key in intron_anno):
            if (ann == "true"):
                intron_anno[key] = ann
        else:
            intron_anno[key] = ann
    fr.close()
    intron_clean = dict()
    fr = open("%s/Intron_attri_%s" % (gtf_path, gtf))
    for info in fr:
        a = info.strip().split("\t")
        C1 = re.sub('.*clean "|\".*', '', a[8])
        cc = re.sub('.*clean_simple "|\".*', '', a[8])
        if lib == 'unstrand':
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], '*')
        else:
            key = "%s_%s_%s_%s" % (a[0], a[3], a[4], a[6])
        if (key in intron_clean):
            if (C1 == "false"):
                intron_clean[key][0] = C1
        else:
            intron_clean[key] = [C1, cc]
    fr.close()
    logging.debug("Finish generating the attributes of each introns\n")

    # Output PI values
    intron_count_file = open("%s/counts/count_intron.txt" % output)
    density_count_file = open("%s/counts/count_all_Density.txt" % output)
    fw = open("%s/intron_PI.txt" % output_path, "w")
    fw.write(
        "Intron_id\tGene_id\tStrand\tChr\tStart\tEnd\tAnnotated\tAttributes\tInclusion_counts\tSkipping_counts\tInclusion_length\tSkipping_length\tPI_Junction\tObserved_counts\tExpected_counts\tPI_Density\n")
    intron_count_line = intron_count_file.readline()
    density_count_line = density_count_file.readline()

    while density_count_line:
        intron_count_list = intron_count_line.strip().split("\t")
        if update != "false":
            attri = intron_count_list[3].split(",")
            if (attri[0] == "false"):
                intron_clean[intron_count_list[0]][0] = "false"
            if (attri[1] == "false"):
                intron_clean[intron_count_list[0]][1] = "false"
        skp = [0] * sample_number
        inc = [0] * sample_number
        PI_J = [0] * sample_number
        sk_l = length - 2 * anchor + 1
        in_l = 2 * (length - 2 * anchor + 1)
        for i in range(0, sample_number):
            inc[i] = str(int(intron_count_list[i * 6 + 7]) + int(intron_count_list[i * 6 + 9]))
            skp[i] = intron_count_list[i * 6 + 11]
            if (update != "false"):
                if (intron_count_list[i * 6 + 8] != intron_count_list[i * 6 + 11] or intron_count_list[i * 6 + 10] !=
                    intron_count_list[i * 6 + 11]):
                    intron_clean[intron_count_list[0]][1] = "false"
            if ((inc[i] + skp[i]) > 0) & ((float(inc[i]) / in_l + float(skp[i]) / sk_l) > 0):
                PI_J[i] = str((float(inc[i]) / in_l) / (float(inc[i]) / in_l + float(skp[i]) / sk_l))
            else:
                PI_J[i] = 'NA'

        density_count_list = density_count_line.strip().split("\t")
        obs = [0] * sample_number
        exp = [0] * sample_number
        PI_D = [0] * sample_number
        for i in range(0, sample_number):
            obs[i] = density_count_list[i + 2]
            exp[i] = density_count_list[i + sample_number + 2]
            if (int(density_count_list[i + sample_number + 2]) < int(density_count_list[i + 2])):
                PI_D[i] = "1"
            elif (int(density_count_list[i + sample_number + 2]) > 0):
                PI_D[i] = str(float(density_count_list[i + 2]) / int(density_count_list[i + sample_number + 2]))
            else:
                PI_D[i] = 'NA'
        intron_count_line = intron_count_file.readline()
        density_count_line = density_count_file.readline()
        ATTRI = "U"
        if (intron_clean[intron_count_list[0]][0] == "false" and intron_clean[intron_count_list[0]][1] == "false"):
            ATTRI = "EI"
        if (intron_clean[intron_count_list[0]][0] == "false" and intron_clean[intron_count_list[0]][1] == "true"):
            ATTRI = "E"
        if (intron_clean[intron_count_list[0]][0] == "true" and intron_clean[intron_count_list[0]][1] == "false"):
            ATTRI = "I"

        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            intron_count_list[0], intron_count_list[1], intron_count_list[2], intron_count_list[4],
            intron_count_list[5], intron_count_list[6], intron_anno[intron_count_list[0]], ATTRI, ",".join(inc),
            ",".join(skp), in_l,
            sk_l,
            ",".join(PI_J), ",".join(obs), ",".join(exp), ",".join(PI_D)))
    fw.close()
    intron_count_file.close()
    density_count_file.close()
    logging.debug("Finished output PI results\n")

logging.debug("#########################################################################\n");

if (Cal == "count"):
    currentTime = time.time()
    runningTime = currentTime - startTime  ## in seconds
    logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime / 3600, (runningTime % 3600) / 60, runningTime % 60))
    logging.debug("Program finished")
    sys.exit()
if not os.path.exists(Comparison) or not Comparison:
    logging.debug("The comparion file was not provided")
    currentTime = time.time()
    runningTime = currentTime - startTime  ## in seconds
    logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime / 3600, (runningTime % 3600) / 60, runningTime % 60))
    logging.debug("Program finished")
    sys.exit()

logging.debug("Calculation of differential spliced introns")
if (not os.path.exists("%s/rMATS_files" % output)):
    os.system("mkdir %s/rMATS_files" % output)
if (not os.path.exists("%s/DEXSeq_files" % output)):
    os.system("mkdir %s/DEXSeq_files" % output)
comparison_param_file = open(Comparison)

for param in comparison_param_file:
    comparison_param = param.strip().split()
    comparison_name = comparison_param[0]
    rep1 = comparison_param[1].split(',')
    rep2 = comparison_param[2].split(',')

    # run rMATS
    fw = open("%s/rMATS_files/rMATS_%s_Junction.txt" % (output, comparison_name), "w")
    fw.write("ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n")
    intron_PI_file = open("%s/intron_PI.txt" % (output_path))
    intron_PI_file.readline()
    for intron_PI_line in intron_PI_file:
        intron_PI_list = intron_PI_line.split("\t")
        SUM = 0
        skp1, inc1 = [0] * len(rep1), [0] * len(rep1)
        skp2, inc2 = [0] * len(rep2), [0] * len(rep2)
        in_level = intron_PI_list[8].split(",")
        sk_level = intron_PI_list[9].split(",")
        for i in range(0, len(rep1)):
            inc1[i] = in_level[int(rep1[i]) - 1]
            skp1[i] = sk_level[int(rep1[i]) - 1]
            SUM += int(inc1[i]) + int(skp1[i])
        for i in range(0, len(rep2)):
            inc2[i] = in_level[int(rep2[i]) - 1]
            skp2[i] = sk_level[int(rep2[i]) - 1]
            SUM += int(inc2[i]) + int(skp2[i])
        if (SUM > 0):
            if (len(comparison_param) > 3 and comparison_param[3] == "pool"):
                fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    intron_PI_list[0], sum(map(lambda x: int(x), inc1)), sum(map(lambda x: int(x), skp1)),
                    sum(map(lambda x: int(x), inc2)), sum(map(lambda x: int(x), skp2)), intron_PI_list[10],
                    intron_PI_list[11]))
            else:
                fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    intron_PI_list[0], ",".join(inc1), ",".join(skp1), ",".join(inc2), ",".join(skp2),
                    intron_PI_list[10], intron_PI_list[11]))
    intron_PI_file.close()
    fw.close()

    cmd = "bash %s/MATS/rMATS.sh -d %s/rMATS_files/rMATS_%s_Junction.txt -o %s/rMATS_files/rMATS_%s_Junction -p %s  -t %s -c %s" % (
        bin_path, output, comparison_name, output, comparison_name, p, analysis, c1)
    logging.debug(cmd)
    os.system(cmd)
    logging.debug("Done running the rMATS for " + comparison_name + " using Junction methods")

    # run DEXSeq
    output_DEXSeq = "%s/DEXSeq_files/DEXSeq_%s" % (output, comparison_name)
    if not os.path.exists(output_DEXSeq):
        os.system("mkdir %s" % output_DEXSeq)

    # generate the count files for DEXSeq
    density_count_file = open("%s/counts/count_all_Density.txt" % output)
    fw = open("%s/DEXSeq_counts.txt" % output_DEXSeq, "w")
    fw.write(
        "ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\tPI_SAMPLE1\tPI_SAMPLE2\tPI_Diff\n")
    for count_info in density_count_file:
        density_count_list = count_info.strip().split("\t")
        exp1 = [0] * len(rep1)
        obs1 = [0] * len(rep1)
        exp2 = [0] * len(rep2)
        obs2 = [0] * len(rep2)
        pi1 = [0] * len(rep1)
        pi2 = [0] * len(rep2)
        if (re.search("NA", count_info)):
            "print error of PI_Density"
        for i in range(0, len(rep1)):
            obs1[i] = density_count_list[int(rep1[i]) + 1]
            exp1[i] = density_count_list[int(rep1[i]) + 1 + sample_number]
            if (int(exp1[i]) < int(obs1[i])):
                pi1[i] = "1"
            elif (int(exp1[i]) == 0):
                pi1[i] = "NA"
            else:
                pi1[i] = str(float(obs1[i]) / (float(exp1[i])))
        for i in range(0, len(rep2)):
            obs2[i] = density_count_list[int(rep2[i]) + 1]
            exp2[i] = density_count_list[int(rep2[i]) + 1 + sample_number]
            if (int(exp2[i]) < int(obs2[i])):
                pi2[i] = "1"
            elif (int(exp2[i]) == 0):
                pi2[i] = "NA"
            else:
                pi2[i] = str(float(obs2[i]) / (float(exp2[i])))

        diff = "NA"
        s1, s2 = [], []
        for i in range(0, len(pi1)):
            if (pi1[i] != "NA"):
                s1.append(float(pi1[i]))
        for i in range(0, len(pi2)):
            if (pi2[i] != "NA"):
                s2.append(float(pi2[i]))
        if (len(s1) > 0 and len(s2) > 0):
            diff = sum(s1) / len(s1) - sum(s2) / len(s2)
        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            density_count_list[0], ",".join(obs1), ",".join(exp1), ",".join(obs2), ",".join(exp2),
            density_count_list[1], density_count_list[1], ",".join(pi1),
            ",".join(pi2), diff))
    density_count_file.close()
    fw.close()

    if len(rep1) + len(rep2) == 2:
        print "DEXSeq does not work on comparison without replicates"
        deseq_counts_file = open("%s/DEXSeq_counts.txt" % (output_DEXSeq))
        intron_PI_file = open("%s/intron_PI.txt" % (output_path))
        rMATS_file = open("%s/rMATS_files/rMATS_%s_Junction/rMATS_Result.txt" % (output, comparison_name))

        fw = open("%s/Result/Diff_%s_intron_PI.txt" % (output, comparison_name), "w")
        fw.write(
            "Intron_id\tGene_id\tStrand\tChr\tStart\tEnd\tAnnotated\tAttributes\tInclusion_counts_SAMPLE1\tSkipping_counts_SAMPLE1\tInclusion_counts_SAMPLE2\tSkipping_counts_SAMPLE2\tInclusion_length\tSkipping_length\tPValue_rMATS\tFDR_rMATS\tPI_Junction_SAMPLE1\tPI_Junction_SAMPLE2\tDiff_PI_Junction\tObserved_counts_SAMPLE1\tExpected_counts_SAMPLE1\tObserved_counts_SAMPLE2\tExpected_counts_SAMPLE2\tPI_Density_SAMPLE1\tPI_Density_SAMPLE2\tDiff_PI_Density\n")

        # skip header line
        deseq_counts_file.readline()
        intron_PI_file.readline()
        rMATS_file.readline()

        deseq_counts_line = deseq_counts_file.readline()
        intron_PI_line = intron_PI_file.readline()
        rMATS_line = rMATS_file.readline()
        inc1 = ['0'] * len(rep1)
        skp1 = ['0'] * len(rep1)
        inc2 = ['0'] * len(rep2)
        skp2 = ['0'] * len(rep2)
        PI1 = ["NA"] * len(rep1)
        PI2 = ["NA"] * len(rep2)
        if (len(comparison_param) > 3 and comparison_param[3] == "pool"):
            inc1 = ['0']
            skp1 = ['0']
            inc2 = ['0']
            skp2 = ['0']
            PI1 = ['NA']
            PI2 = ['NA']
        while rMATS_line:
            deseq_counts_list = deseq_counts_line.strip().split("\t")
            intron_PI_list = intron_PI_line.strip().split("\t")
            rMATS_list = rMATS_line.strip().split("\t")
            if (intron_PI_list[0] != rMATS_list[0]):
                fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\t%s\t%s\tNA\t%s\t%s\n" % (
                    "\t".join(intron_PI_list[0:8]), ",".join(inc1), ",".join(skp1), ",".join(inc2), ",".join(skp2),
                    rMATS_list[5], rMATS_list[6],
                    ",".join(PI1), ",".join(PI2), "\t".join(deseq_counts_list[1:5]),
                    "\t".join(deseq_counts_list[7:10])))
                deseq_counts_line = deseq_counts_file.readline()
                intron_PI_line = intron_PI_file.readline()
            if (intron_PI_list[0] == rMATS_list[0]):
                fw.write("%s\t%s\t%s\t%s\n" % (
                    "\t".join(intron_PI_list[0:8]), "\t".join(rMATS_list[1:]), "\t".join(deseq_counts_list[1:5]),
                    "\t".join(deseq_counts_list[7:10])))
                deseq_counts_line = deseq_counts_file.readline()
                intron_PI_line = intron_PI_file.readline()
                rMATS_line = rMATS_file.readline()
        deseq_counts_file.close()
        intron_PI_file.close()
        rMATS_file.close()
        fw.close()
        logging.debug("Output the result of differential spliced intron analysis of " + comparison_name + "\n")
    else:
        in1 = [str(int(x) + 2) for x in rep1]
        in2 = [str(int(x) + 2) for x in rep2]

        cmd = "cut -f 1,2 %s/count_intron.txt > %s/intron_id.txt" % (count_path, output_DEXSeq)
        os.system(cmd)
        logging.debug(cmd)
        cmd = "cut -f %s,%s %s/count_all_Density.txt > %s/intron_count.txt" % (
            ",".join(in1), ",".join(in2), count_path, output_DEXSeq)
        os.system(cmd)
        logging.debug(cmd)
        cmd = "paste %s/intron_id.txt %s/intron_count.txt  > %s/intron_data.txt" % (
            output_DEXSeq, output_DEXSeq, output_DEXSeq)
        os.system(cmd)
        logging.debug(cmd)

        DEX_Gene = dict()
        count_exon_file = open("%s/count_exon.txt" % count_path)
        for count_exon in count_exon_file:
            count_exon_list = count_exon.strip().split("\t")
            DEX_Gene[count_exon_list[0]] = [0] * (len(rep1) + len(rep2))
            for j in range(0, len(rep1)):
                DEX_Gene[count_exon_list[0]][j] = int(count_exon_list[int(rep1[j])])
            for j in range(0, len(rep2)):
                DEX_Gene[count_exon_list[0]][len(rep1) + j] = int(count_exon_list[int(rep2[j])])
        count_exon_file.close()

        fw = open("%s/alter_data.txt" % (output_DEXSeq), "w")
        intron_data_file = open("%s/intron_data.txt" % (output_DEXSeq))
        for intron_data in intron_data_file:
            intron_data_list = intron_data.strip().split("\t")
            gene = intron_data_list[1].split(",")
            for j in range(2, len(intron_data_list)):
                alt = 0
                for g in gene:
                    if g in DEX_Gene:
                        alt += DEX_Gene[g][j - 2]
                fw.write("%s\t" % (alt))
            fw.write("\n")
        intron_data_file.close()
        fw.close()

        fw = open("%s/SampleData.txt" % output_DEXSeq, "w")
        fw.write("\tcondition\tlibType\n")
        libType = "single-end"
        if (read == "P"):
            libType = "paired-end"

        for i in range(0, len(in1)):
            fw.write("sample1_%s\tSAMPLE1\t%s\n" % (i, libType))
        for i in range(0, len(in2)):
            fw.write("sample2_%s\tSAMPLE2\t%s\n" % (i, libType))
        fw.close()

        # run DEXSeq
        cmd = "Rscript %s/DEXSeq.R %s/intron_data.txt %s/SampleData.txt %s/alter_data.txt %s/DEXSeq.txt %s" % (
            bin_path, output_DEXSeq, output_DEXSeq, output_DEXSeq, output_DEXSeq, p)
        logging.debug(cmd)
        os.system(cmd)
        logging.debug("Done running the DEXSeq for " + comparison_name + " using Denstiy methods")

        # generate final output
        deseq_file = open("%s/DEXSeq.txt" % (output_DEXSeq))
        deseq_counts_file = open("%s/DEXSeq_counts.txt" % (output_DEXSeq))
        intron_PI_file = open("%s/intron_PI.txt" % (output_path))
        rMATS_file = open("%s/rMATS_files/rMATS_%s_Junction/rMATS_Result.txt" % (output, comparison_name))

        fw = open("%s/Diff_%s_intron_temp.txt" % (output_DEXSeq, comparison_name), "w")
        fw.write(
            "Intron_id\tGene_id\tStrand\tChr\tStart\tEnd\tAnnotated\tAttributes\tInclusion_counts_SAMPLE1\tSkipping_counts_SAMPLE1\tInclusion_counts_SAMPLE2\tSkipping_counts_SAMPLE2\tInclusion_length\tSkipping_length\tPValue_rMATS\tFDR_rMATS\tPI_Junction_SAMPLE1\tPI_Junction_SAMPLE2\tDiff_PI_Junction\tObserved_counts_SAMPLE1\tExpected_counts_SAMPLE1\tObserved_counts_SAMPLE2\tExpected_counts_SAMPLE2\tPValue_DEXSeq\tFDR_DEXSeq\tPI_Density_SAMPLE1\tPI_Density_SAMPLE2\tDiff_PI_Density\n")

        deseq_file.readline()
        deseq_counts_file.readline()
        intron_PI_file.readline()
        rMATS_file.readline()

        deseq_line = deseq_file.readline()
        deseq_counts_line = deseq_counts_file.readline()
        intron_PI_line = intron_PI_file.readline()
        rMATS_line = rMATS_file.readline()

        inc1 = ['0'] * len(rep1)
        skp1 = ['0'] * len(rep1)
        inc2 = ['0'] * len(rep2)
        skp2 = ['0'] * len(rep2)
        PI1 = ["NA"] * len(rep1)
        PI2 = ["NA"] * len(rep2)
        rMATS_list = rMATS_line.strip().split("\t")
        len1 = rMATS_list[5]
        len2 = rMATS_list[6]
        if (len(comparison_param) > 3 and comparison_param[3] == "pool"):
            inc1 = ['0']
            skp1 = ['0']
            inc2 = ['0']
            skp2 = ['0']
            PI1 = ['NA']
            PI2 = ['NA']

        while (intron_PI_line):
            deseq_list = deseq_line.strip().split("\t")
            deseq_counts_list = deseq_counts_line.strip().split("\t")
            intron_PI_list = intron_PI_line.strip().split("\t")
            rMATS_list = rMATS_line.strip().split("\t")
            if (intron_PI_list[0] != rMATS_list[0]):
                fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\t%s\t%s\tNA\t%s\t%s\t%s\n" % (
                    "\t".join(intron_PI_list[0:8]), ",".join(inc1), ",".join(skp1), ",".join(inc2), ",".join(skp2),
                    len1, len2,
                    ",".join(PI1), ",".join(PI2), "\t".join(deseq_counts_list[1:5]), "\t".join(deseq_list[6:8]),
                    "\t".join(deseq_counts_list[7:10])))
                deseq_line = deseq_file.readline()
                deseq_counts_line = deseq_counts_file.readline()
                intron_PI_line = intron_PI_file.readline()
            if (intron_PI_list[0] == rMATS_list[0]):
                fw.write("%s\t%s\t%s\t%s\t%s\n" % (
                    "\t".join(intron_PI_list[0:8]), "\t".join(rMATS_list[1:]), "\t".join(deseq_counts_list[1:5]),
                    "\t".join(deseq_list[6:8]), "\t".join(deseq_counts_list[7:10])))
                deseq_line = deseq_file.readline()
                deseq_counts_line = deseq_counts_file.readline()
                intron_PI_line = intron_PI_file.readline()
                rMATS_line = rMATS_file.readline()
        deseq_file.close()
        deseq_counts_file.close()
        intron_PI_file.close()
        rMATS_file.close()
        fw.close()
        gtf_path = "%s/gtf_files" % output
        cmd = "python %s/Transcript.py --FPKM  %s --intron %s/Intron_transcript.txt --output %s/FPKM_intron_transcript.txt" % (
            bin_path, FPKM, gtf_path, output_DEXSeq)
        logging.debug(cmd)
        os.system(cmd)
        cmd = "Rscript %s/RP_value.R %s/Diff_%s_intron_temp.txt 100 %s/rank_product_test_%s.txt %s/Result/Diff_%s_intron_PI.txt %s %s %s/Result/Increase_%s_intron_PI.txt %s/Result/Decrease_%s_intron_PI.txt" % (
            bin_path, output_DEXSeq, comparison_name, output_DEXSeq, comparison_name, output, comparison_name,
            F_deltaPI, F_FDR, output, comparison_name, output,
            comparison_name)
        logging.debug(cmd)
        os.system(cmd)
        logging.debug("Output the result of differential spliced intron analysis of " + comparison_name + "\n")
comparison_param_file.close()

currentTime = time.time()
runningTime = currentTime - startTime  ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime / 3600, (runningTime % 3600) / 60, runningTime % 60))
logging.debug("Program finished")
sys.exit(0)
