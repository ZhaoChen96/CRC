#!/usr/bin/python
#######################################################################
# File Name: tgev.py
# Author: Qinglan Li
# mail: liqinglan@whu.edu.cn
# Usage: 
# Created Time: Thu Sep 22 20:56:03 2016
# Description: This is a python script to help you download data from 
#              TCGA and analyize the expression of specific gene in 
#              different type of tissue(normal and cancer). 
# History: Fri Nov 11 2016 revised, gathering all three python scripts and adding the function of gunzip.
#######################################################################

import sys
import os
import getopt

########################## default ###########################
##setting python file path
python_file_name = sys.argv[0]
python_path_part = len(python_file_name.split('/'))
python_path_list = []
if python_path_part == 1:
    python_file_path = '.'
else:    
    for n in range(0,python_path_part-1):
        python_path_list.append(python_file_name.split('/')[n])
    python_file_path = '/'.join((python_path_list))

########################## optional ###########################
##setting work path
def get_work_path(manifest):
    path_part = len(manifest.split('/'))
    if path_part == 1:
        work_path = '.'
        clean_manifest = 'clean_manifest.txt'
        json_type = '%s/json_type.txt' % (python_file_path)
        json = 'clean.json'
    else:
        path_list = []
        for N in range(0,path_part-1):
            path_list.append(manifest.split('/')[N])
        work_path = '/'.join((path_list))
        clean_manifest = '%s/clean_manifest.txt' % (work_path)
        json_type = '%s/json_type.txt' % (python_file_path)
        json = '%s/clean.json' % (work_path)     #work_path is character indicating the work dir
    return (work_path,manifest,clean_manifest,json_type,json)

########################## optional ###########################
##module of JSON file preperation
#locating the relative file position
def get_JOSN(Tuple):
#get args from get_work_path
    work_path = Tuple[0]
    manifest = Tuple[1]
    clean_manifest = Tuple[2]
    json_type = Tuple[3]
    json = Tuple[4]
    i1 = open(manifest,'r')
    i2 = open(json_type,'r')
    o1 = open(clean_manifest,'w')
    o2 = open(json,'w')
#creating cleaned manifest
    o1.write(i1.readline())
    U = []
    for l1 in i1:
        UUID = l1.split('\t')[0]
        name = l1.split('\t')[1]
        if name.split('.')[1] == 'FPKM': # check whether this is a FPKM file we need
            U.append(UUID)
            o1.write(l1)
#calculating the tab numbers in json_type
    j = i2.readlines()
    all_len = len(j[5])
    remain_len = len(j[5].strip('\t'))
    tab_number = all_len - remain_len + 1
    tab = tab_number*'\t'
#writing new json file
    for n in range(0,6):
        o2.write(j[n])
    num = 0
    for u in U:
        num += 1
        if num < len(U):
            line = '%s"%s",\n' % (tab,u)
            o2.write(line)
        elif num == len(U):
            line = '%s"%s"\n' % (tab,u)
            o2.write(line)
    for n in range(6,len(j)):
        o2.write(j[n])

    i1.close()
    i2.close()
    o1.close()
    o2.close()
##module of creating file metadata
    os.system("curl --request POST --header 'Content-Type: application/json' --data @%s 'https://gdc-api.nci.nih.gov/files' > %s/File_message.txt" % (json,work_path))
##module of downloading data
    os.system('gdc-client download -m %s' % (clean_manifest))

########################## optional ###########################
##module of dividing data type
def gunzip_divide(Tuple):
    work_path = Tuple[0]
    file_massage = '%s/File_message.txt' % (work_path)
    #locating and gunzip file
    for d in os.listdir(work_path):
        if os.path.isdir('%s/%s' % (work_path,d)):
            for f in os.listdir('%s/%s' % (work_path,d)):
                if os.path.isfile('%s/%s/%s' % (work_path,d,f)):
                    os.system('gunzip %s/%s/%s' % (work_path,d,f))
    #choose group of cancer and normal
    os.system('mkdir %s/raw_data' % (work_path))
    os.system('mkdir %s/cancer' % (work_path))
    os.system('mkdir %s/normal' % (work_path))
    i_massage = open(file_massage,'r')
    i_massage.readline()
    for l in i_massage:
        sample_type = l.split('\t')[2]
        UUID = l.split('\t')[7]
        if os.path.isdir('%s/%s' % (work_path,UUID)):
            if 'Normal' in sample_type.split(' '):
                os.system('mv %s/%s/*.txt %s/normal' % (work_path,UUID,work_path))
                os.system('rm -rf %s/%s' % (work_path,UUID))
            else:
                os.system('mv %s/%s/*.txt %s/cancer' % (work_path,UUID,work_path))
                os.system('rm -rf %s/%s' % (work_path,UUID))
    os.system('mv %s/normal %s/raw_data' % (work_path,work_path))
    os.system('mv %s/cancer %s/raw_data' % (work_path,work_path))

########################## optional ###########################
def gunzip(path):   #the path of dir saved unzip data
    for d in os.listdir(path):
        if os.path.isdir('%s/%s' % (path,d)):
            for f in os.listdir('%s/%s' % (path,d)):
                if os.path.isfile('%s/%s/%s' % (path,d,f)):
                    os.system('gunzip %s/%s/%s' % (path,d,f))
                    os.system('mv %s/%s/*.txt ..' % (path,d))
                    os.system('rm -rf %s/%s' % (path,d))

########################## optional ###########################
##gene FPKM statistics
def gene_stat(genename,Tuple,other_dir_path):
    work_path = Tuple[0]
    default_path = '%s/raw_data' % (work_path)
    if os.path.isdir(default_path):
        dirname = default_path
    else:
        dirname = other_dir_path    #path of other data saved dir
    annotation_file = '%s/offcial_ensembl_id.txt' % (python_file_path)
    #getting gene ensembl ID using popen
    popen_gene_id = os.popen('grep %s %s' % (genename,annotation_file))
    gene_id = popen_gene_id.readline().split('\t')[1].split('\n')[0]    #ensembl name
    #normal
    normal = []
    cancer = []
    for f in os.listdir('%s/normal' % (dirname)):
        print f
        if f == 'annotations.txt':
            continue
        popen_gene_FPKM = os.popen('grep %s %s/normal/%s' % (gene_id,dirname,f))
        normal_FPKM = popen_gene_FPKM.readline().split('\t')[1].split('\n')[0] #normal FPKM value
        normal.append(normal_FPKM)
        popen_gene_FPKM.close()
    #cancer
    for f in os.listdir('%s/cancer' % (dirname)):
        print f
        if f == 'annotations.txt':
            continue
        popen_gene_FPKM = os.popen('grep %s %s/cancer/%s' % (gene_id,dirname,f))
        cancer_FPKM = popen_gene_FPKM.readline().split('\t')[1].split('\n')[0] #cancer FPKM value
        cancer.append(cancer_FPKM)
        popen_gene_FPKM.close()
    #output
    num_cancer = len(cancer)
    num_normal = len(normal)
    n = num_cancer - num_normal
    for e in range(0,n):
        normal.append('NA')
    result = open('%s/../result_%s.txt' % (dirname,genename),'w')
    result.write('\t'.join(('normal','cancer','\n')))
    for e in range(0,num_cancer):
        result.write('\t'.join((normal[e],cancer[e],'\n')))
    result.close()

######################## usage ########################
def usage():
    print ''' 
 Usage: python3 Your_dir/tgev.py [-h|-a|-d] [--all|--gunzip] args...
 
 Options:
    -h  None
        Calling help.
    -a  [data_saved_path(absolute or relative)] [offical_gene_name]
        Doing analysis. If you have saved your interseted data and you want to get the statstical file for specific gene, you can call this option. The dir of saving data should contain two already divided dir named 'normal' and 'cancer'.
    -d  [manifest_file_path]
        Downloading and dividing your interested data without analysis.
    -A  [manifest_file_path(absolute or relative)] [offical_gene_name]
        Doing whole process of this task. All the thing you need to do is just providing the path of manifest file and the gene name you interested in.
    --gunzip [data_saved_path(absolute or relative)]
        Gunzip the downloaded data. This function was design for who downloaded his interested data by himself, it will help you gunzip your data without dividing.

 Version 1.1
 Qinglan Li
 '''

######################## main ########################
if '__main__' == __name__:
    opts,args = getopt.getopt(sys.argv[1:],'hAad',["gunzip"])
    if not(opts):
        usage()
        sys.exit(0)
    for name,value in opts:
        if name in ('-h'):
            usage()
            sys.exit(0)
        elif name in ('-A'):
            Tuple = get_work_path(args[0])
            get_JOSN(Tuple)
            gunzip_divide(Tuple)
            gene_stat(args[1],Tuple,' ')
            print 'Done!'
            sys.exit(0)
        elif name in ('-d'):
            Tuple = get_work_path(args[0])
            print(Tuple[0]) 
            get_JOSN(Tuple)
            gunzip_divide(Tuple)
            print 'Done!'
            sys.exit(0)
        elif name in ('-a'):
            gene_stat(args[1],' ',args[0])
            print 'Done!'
            sys.exit(0)
        elif name in ('--gunzip'):
            gunzip(args[0])
            print 'Done!'
            sys.exit(0)














