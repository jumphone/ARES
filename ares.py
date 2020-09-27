#----------------------------------------------------------------------------
# ARES
# Annotation-free toolkit for identifying Rna Editing Sites
# Author: Feng Zhang
# Date: 2020.9
# Requirements:
# python3=3.8.2
# pysam=0.16.0.1
# numpy=1.19.1
# bedtools=v2.26.0
# blat=v.36x7
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
import re
import os
import sys
import time
import pysam
import numpy as np
import subprocess 
from multiprocessing import Process, Queue

#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
print('Started !\n')
print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def loadRef(fa_in_path):
    fref=open(fa_in_path)
    chrom={}
    chrr=''
    line=fref.read()
    line=line.split('>')
    for seq in line:
        if ' ' in seq:
            chrr=seq[0:seq.find(' ')]
        else:
            chrr=seq[0:seq.find('\n')]
        chrom[chrr]=seq[seq.find('\n'):].replace('\n','')
    #0 base
    fref.close()
    return chrom
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def bam2zz(bam_in_path=0,fa_in_path=0,zz_out_path=0):
    #-----------------------------------------------------------
    #Load refgenome
    chrom=loadRef(fa_in_path=fa_in_path)
    #------------------------------------------------------------
    #
    #donee:dolist;lst:donumlist
    def doCG(a):
        donee=[]
        lst=re.findall( '(\d+|\+|-|\*|/)', a )
        for i in a:
            if i == 'I' or i == 'D' or i== 'M' or i=='S' or i=='P' or i=='N' :
                donee.append(i)
        return donee,lst
	#donefunction

    def doneCG(CG,chrr,pos,seq,qseq):#pos is 1 base
        donee,lst=doCG(CG)
        #print(donee, lst)
        errorsite=''
        intersite=''
        quasite=''
        locsite=''
        pieceloc=''
        refseq=''
        seqseq=''
        refpos=int(pos)-1
        seqpos=0
        step=0
        while step<len(donee):
            if donee[step]=='I':
                seqpos=seqpos+int(lst[step])
            elif donee[step]=='D':
                refpos=refpos+int(lst[step])
            elif donee[step]=='N':
                refpos=refpos+int(lst[step])
            elif donee[step]=='S':
                seqpos=seqpos+int(lst[step])
            elif donee[step]=='M':
                refseq=refseq+chrom[chrr][refpos:refpos+int(lst[step])]
                seqseq=seqseq+seq[seqpos:seqpos+int(lst[step])]
                j=refpos
                jj=seqpos
                while j<refpos+int(lst[step]):
                    try:
                        if chrom[chrr][j].upper() != seq[jj].upper() and chrom[chrr][j].upper() !='N' and seq[jj].upper() != 'N':
                            errorsite=errorsite+chrom[chrr][j].upper()+seq[jj].upper()+':'+str(j+1)+';'
                            quasite=quasite+','+str(ord(qseq[jj]))
                            locsite=locsite+','+str(jj+1)
                            pieceloc=pieceloc+','+str( min( jj+1-seqpos,seqpos+int(lst[step])-jj  ) )
                    except ValueError:
                        pass
                    j=j+1
                    jj=jj+1
                intersite=intersite+str(refpos+1)+':'+str(refpos+int(lst[step]))+';'
                refpos=refpos+int(lst[step])
                seqpos=seqpos+int(lst[step])
            step=step+1
        refseq=refseq.upper()
        seqseq=seqseq.upper()
        return refseq,seqseq,errorsite,intersite,quasite,locsite,pieceloc
	#################################################################################		
    
    fi = pysam.AlignmentFile(bam_in_path, "rb") #bam
    fo=open(zz_out_path,'w') #zz
    
    fo.write('#chr\tsam_flag\tsam_mapq\tread_interval\tmismatch_list\tmismatch_quality\tmismatch_read_pos\tseq\tread_name\tfragment_end_dist\n') 
    
    for line in fi:
        this_name=str(line.query_name)
        this_chr=str(line.reference_name)
        this_CG=str(line.cigarstring)
        this_pos=str(line.pos+1)
        this_seq=str(line.seq)
        this_qseq=str(line.qual)
        this_mapq=str(line.mapq)
        this_flag=str(line.flag)
        
        if this_chr!='*' and this_CG!='*' and this_CG!='None' :
                refseq,seqseq,errorsite,intersite,quasite,locsite,pieceloc=doneCG(this_CG, this_chr, this_pos, this_seq, this_qseq)
                quasite=quasite[1:]
                locsite=locsite[1:]
                pieceloc=pieceloc[1:]
				
                if len(intersite[0:-1])==0:
                    intersite='*;'
                if len(errorsite[0:-1])==0:
                    errorsite='*;'
                if len(quasite)==0:
                    quasite='*'
                if len(locsite)==0:
                    locsite='*'
                if len(pieceloc)==0:
                    pieceloc='*'
				
                if len(bin(int(this_flag)))>=9 and bin(int(this_flag))[-7]=='1':
                    this_name=this_name+'_1'
                elif len(bin(int(this_flag)))>=10 and bin(int(this_flag))[-8]=='1':
                    this_name=this_name+'_2'
                this_out=this_chr+'\t'+this_flag+'\t'+this_mapq+'\t'+intersite[0:-1]+'\t'+errorsite[0:-1]+'\t'+quasite+'\t'+locsite+'\t'+this_seq+'\t'+this_name+'\t'+pieceloc+'\n'
                if errorsite !='*;':
                    fo.write(this_out)	
    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def zz2snv(zz_in_path=0,bed_out_path=0):

    fi=open(zz_in_path)
    fo = open(bed_out_path,'w')

    fo.write('#chr\tstart_0\tend_1\tmismatch\tquality\tfragment_end_dist\n')

    header=fi.readline().replace('#','').rstrip().split('\t')
    for line in fi:
        seq=line.rstrip().split('\t')
        this_mismatch_list=seq[header.index('mismatch_list')]
        if this_mismatch_list !='*':
            this_mismatch_list=this_mismatch_list.split(';')
            this_chr=seq[header.index('chr')]
            this_fragment_end_dist=seq[header.index('fragment_end_dist')]
            this_fragment_end_dist=this_fragment_end_dist.split(',')
            this_mismatch_quality=seq[header.index('mismatch_quality')]
            this_mismatch_quality=this_mismatch_quality.split(',')
            this_read_name=seq[header.index('read_name')]
            i=0
            while i <len(this_mismatch_list):
                i_mismatch_type=this_mismatch_list[i].split(':')[0]
                i_mismatch_pos=this_mismatch_list[i].split(':')[1]
                i_fragment_end_dist=this_fragment_end_dist[i]
                i_mismatch_quality=this_mismatch_quality[i]
                i_out=[this_chr, str(int(i_mismatch_pos)-1),i_mismatch_pos, i_mismatch_type, i_mismatch_quality, i_fragment_end_dist]
                fo.write('\t'.join(i_out)+'\n')
                i=i+1

    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def filterSnv(bed_in_path=0, bed_out_path=0, CUT=10):
    CUT=CUT
    bed_in_path=bed_in_path
    bed_out_path=bed_out_path

    TYPE={}
    fi=open(bed_in_path)
    header=fi.readline
    for line in fi:
        seq=line.rstrip().split('\t')
        this_tag=seq[0]+':'+seq[2]
        if this_tag in TYPE:
            TYPE[this_tag].append(seq[3])
        else:
            TYPE[this_tag]=[seq[3]]
    fi.close()

    ###################
    fi=open(bed_in_path)
    QUA=[] #Base Quality
    FED=[] #Fragment End Dist
    header=fi.readline()
    for line in fi:
        seq=line.rstrip().split('\t')
        QUA.append(float(seq[4]))
        FED.append(float(seq[5]))
    fi.close()

    Qcut=np.percentile(QUA, CUT)
    Fcut=np.percentile(FED, CUT)

    print('Qcut Fcut\n')
    print(Qcut, Fcut)

    #####################
    fi=open(bed_in_path)
    fo=open(bed_out_path,'w')
    header=fi.readline()
    PRINTED=set()
    for line in fi:
        seq=line.rstrip().split('\t')
        this_tag=seq[0]+':'+seq[2]
        if float(seq[4])>Qcut and float(seq[5])>Fcut and len(set(TYPE[this_tag]))==1:
            this_snv='\t'.join(seq[0:4])
            if this_snv not in PRINTED:
                PRINTED.add(this_snv)
                fo.write(this_snv+'\n')
    fi.close()
    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def sortBed(bed_in_path=0, bed_out_path=0, bedtools_path=0):
    this_step =subprocess.Popen(' '.join([bedtools_path,'sort -i', bed_in_path, '>', bed_out_path]),shell=True)
    this_step.wait()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def dupletCluster(bed_in_path=0,bed_out_path=0,cluster_distance=200,cluster_size=2):
    fi=open(bed_in_path)
    fo=open(bed_out_path,'w')
    tmp='chr0:0:AA'
    limitdistance=int(cluster_distance)
    limitnum=int(cluster_size)
    lst=[]
    cluster_id=1
    for line in fi:
            seq=line.split('\t')
            tmpseq=tmp.split(':')
            if seq[0]==tmpseq[0] and int(seq[2])-int(tmpseq[1])<=limitdistance and seq[3]==tmpseq[2]:
                lst.append(line)
            else:
                if len(lst)>=limitnum:
                        begin=float(lst[0].split('\t')[1])
                        end=float(lst[-1].split('\t')[2])
                        density=len(lst)/(end-begin)
                        for one in lst:
                            #fo.write(one[0:-1]+'\t'+str(len(lst))+'\t'+str(density)+'\tC'+str(cluster_id)+'E\n')
                            fo.write(one[0:-1]+'\t'+str(len(lst))+'\tC'+str(cluster_id)+'E\n')
                        cluster_id = cluster_id +1
                lst=[]
                lst.append(line)
            tmp=seq[0]+':'+seq[2]+':'+seq[3]
    if len(lst)>=limitnum:
            begin=float(lst[0].split('\t')[1])
            end=float(lst[-1].split('\t')[2])
            density=len(lst)/(end-begin)
            for one in lst:
                #fo.write(one[0:-1]+'\t'+str(len(lst))+'\t'+str(density)+'\tC'+str(cluster_id)+'E\n')
                fo.write(one[0:-1]+'\t'+str(len(lst))+'\tC'+str(cluster_id)+'E\n')
            cluster_id = cluster_id +1
    fi.close()
    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def getClusterBed(bed_in_path=0, bed_out_path=0):

    CLST_POS={}
    CLST_CHR={}
    CLST_TYPE={}
    CLST_NUM={}
    CLST_LIST=[]
    CLST_SET=set()

    fi=open(bed_in_path)
    for line in fi:
        seq=line.rstrip().split('\t')
        this_id=seq[5]
        CLST_SET.add(this_id)

    TMP=[]
    for this_id in CLST_SET:
        TMP.append([int(this_id.replace('C','').replace('E','')), this_id])
    TMP.sort()

    for one in TMP:
        this_id=one[1]
        CLST_LIST.append(this_id)
        CLST_POS[this_id]=[]
        CLST_CHR[this_id]=[]
        CLST_NUM[this_id]=[]
        CLST_TYPE[this_id]=[]

    fi=open(bed_in_path)
    fo=open(bed_out_path,'w')

    i=0
    for line in fi:
        seq=line.rstrip().split('\t')
        this_start=int(seq[1])
        this_end  =int(seq[2])
        this_chr  =seq[0]
        this_type =seq[3]
        this_num  =seq[4]
        this_id   =seq[5]

        CLST_POS[this_id].append(this_start)
        CLST_POS[this_id].append(this_end)
        CLST_CHR[this_id].append(this_chr)
        CLST_TYPE[this_id].append(this_type)
        CLST_NUM[this_id].append(this_num)
    
        if i % 1000==1:
            print(i)
        i=i+1


    for this_id in CLST_LIST:
        this_start=str(min(CLST_POS[this_id]))
        this_end  =str(max(CLST_POS[this_id]))
        this_chr  =CLST_CHR[this_id][0]
        this_type =CLST_TYPE[this_id][0]
        this_num  =CLST_NUM[this_id][0]
        fo.write(this_chr+'\t'+this_start+'\t'+this_end+'\t'+this_type+'\t'+this_num+'\t'+this_id+'\n')

    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def bed2flankfa(bed_in_path=0, ref_in_path=0, fa_out_path=0, AROUND=0, SEP='r'):

    #-----------------------------------------------------------
    #Load refgenome
    chrom=loadRef(fa_in_path=ref_in_path)
    #------------------------------------------------------------
    fi=open(bed_in_path)
    fo=open(fa_out_path,'w')
    for line in fi:
        seq=line.rstrip().split('\t')
        this_tag=SEP.join(seq)
        this_chr=seq[0]
        this_start=int(seq[1])
        this_end=int(seq[2])

        this_add_start=max([this_start-AROUND,0])
        this_add_end=this_end+AROUND

        this_seq_before = chrom[this_chr][this_add_start:this_start].upper()
        this_seq_after = chrom[this_chr][this_end:this_add_end].upper()

        fo.write('>'+this_tag+SEP+'BEFORE\n')
        fo.write(this_seq_before+'\n')
        fo.write('>'+this_tag+SEP+'AFTER\n')
        fo.write(this_seq_after+'\n')

    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def bed2addfa(bed_in_path=0, ref_in_path=0, fa_out_path=0, AROUND=0, SEP='q'):
    #-----------------------------------------------------------
    #Load refgenome
    chrom=loadRef(fa_in_path=ref_in_path)
    #------------------------------------------------------------
    
    fi=open(bed_in_path)
    fo=open(fa_out_path,'w')
    for line in fi:
        seq=line.rstrip().split('\t')
        this_chr=seq[0]
        this_start=int(seq[1])
        this_end=int(seq[2])
        this_tag=SEP.join(seq)


        this_add_start=max([this_start-AROUND,0])
        this_add_end=this_end+AROUND
 
        this_seq = chrom[this_chr][this_add_start:this_add_end].upper()
        fo.write('>'+this_tag+'\n')
        fo.write(this_seq+'\n')

    fo.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def blatWorker(q, ARG):
    ref_in_path=ARG[0]
    query_in_path=ARG[1]
    map_out_path=ARG[2]
    blat_path=ARG[3]

    this_step =subprocess.Popen(' '.join([blat_path,'-out=blast9','-minScore=25','-minIdentity=70', '-fastMap', ref_in_path, query_in_path, 
	                                   map_out_path+'.tmpOut', '>', map_out_path+'.tmpInfo' ]),shell=True)
    this_step.wait()
    fmap=open(map_out_path,'a')
    ftmp=open(map_out_path+'.tmpOut')

    for l_tmp in ftmp:
        if l_tmp[0]!='#':
            seq_tmp=l_tmp.rstrip().split('\t')
            fmap.write(l_tmp)

    ftmp.close()
    fmap.close()
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def blatAlign(ref_in_path=0, query_in_path=0, out_dir=0, blat_path=0, CPU=10):
    out_dir=out_dir+'/'
    try:  
        os.mkdir(out_dir)  
    except OSError as error:  
        print(error)

    count_out_path=out_dir+'/finished_count.txt'

    ##########################
    REF={}
    fi=open(ref_in_path)
    l1=fi.readline()
    l2=fi.readline()
    while l1 !='':
        this_id=l1.rstrip().split('r')[5]
        REF[this_id]=''
        l1=fi.readline()
        l2=fi.readline()
    fi.close()

    fi=open(ref_in_path)
    l1=fi.readline()
    l2=fi.readline()
    while l1 !='':
        this_id=l1.rstrip().split('r')[5]
        REF[this_id]+=l1+l2
        l1=fi.readline()
        l2=fi.readline()
    fi.close()
    #######################
    #import subprocess

    QUERY={}
    fi=open(query_in_path)
    l1=fi.readline()
    l2=fi.readline()
    while l1 !='':
        this_id=l1.rstrip().split('q')[5]
        QUERY[this_id]=''
        l1=fi.readline()
        l2=fi.readline()
    fi.close()


    fi=open(query_in_path)
    l1=fi.readline()
    l2=fi.readline()
    while l1 !='':
        this_id=l1.rstrip().split('q')[5]
        QUERY[this_id]+=l1+l2
        l1=fi.readline()
        l2=fi.readline()
    fi.close()
    #######################

    JOBS=[]
    queue = Queue()
    i=0
    for this_id in QUERY:

        this_id=this_id

        this_ref=REF[this_id]
        this_query=QUERY[this_id]

        ##################################
        this_job_id     =str(len(JOBS))
        this_ref_path   =out_dir+'/r.'+this_job_id+'.fa'
        this_query_path =out_dir+'/q.'+this_job_id+'.fa'
        this_out_path   =out_dir+'/o.'+this_job_id+'.blast9'

        open(this_ref_path,  'w').write(this_ref)
        open(this_query_path,'w').write(this_query)
        if i % 100 ==1:
            open(count_out_path,'a').write(str(i)+'\n')
        i=i+1

        ##################
        #blatWorker(ref_in_path=0, query_in_path=0, map_out_path=0, blat_path=0)
        p=Process(target=blatWorker, args=(queue, [this_ref_path, this_query_path, this_out_path, blat_path] ))
        p.start()
        #this_job =subprocess.Popen(' '.join(['python3',BLAT_WORKER, this_ref_path, this_query_path, this_out_path]),shell=True)
        JOBS.append(p)
        ##################################
        if len(JOBS)>=CPU:
            for this_job in JOBS:
                this_job.join()
            JOBS=[]
        ##################################
        ##################################

    for this_job in JOBS:
        this_job.join()
#----------------------------------------------------------------------------



#PATH='/home/Lilab/Docker/zhangfeng/database/RNAseq/U87MG_ADAR/hisat2/SRR388226.ADAR.CTRL.airOut/'

#blatAlign(ref_in_path=PATH+'/f4_SNV_sorted_duplet.bed.c.bed.20000.fa', 
#     query_in_path=PATH+'/f4_SNV_sorted_duplet.bed.fa', 
#     out_dir='/home/Lilab/Docker/zhangfeng/database/RNAseq/U87MG_ADAR/hisat2/ARES/TMP_OUT', 
#     blat_path='/usr/local/bin/blat', 
#     CPU=10)

print('''
    
bam_in_path   =  sys.argv[1]
ref_in_path   =  sys.argv[2]
OUT_DIR       =  sys.argv[3]
bedtools_path =  sys.argv[4]
blat_path     =  sys.argv[5]

	''')


bam_in_path   =  sys.argv[1]
ref_in_path   =  sys.argv[2]
OUT_DIR       =  sys.argv[3]
bedtools_path =  sys.argv[4]
blat_path     =  sys.argv[5]

CPU=10
CUT=10
CLUSTER_DISTANCE=200
SITE_ADD_AROUND=20
CLUSTER_FLANK_AROUND=20000

OUT_DIR=OUT_DIR+'/'
try:  
    os.mkdir(OUT_DIR)  
except OSError as error:  
    print(error)


bam2zz(bam_in_path=bam_in_path, fa_in_path=ref_in_path, zz_out_path=OUT_DIR+'/f1_read.zz')
zz2snv(zz_in_path=OUT_DIR+'/f1_read.zz', bed_out_path=OUT_DIR+'/f2_snv.bed')
filterSnv(bed_in_path=OUT_DIR+'/f2_snv.bed', bed_out_path=OUT_DIR+'/f3_snv_filtered.bed', CUT=CUT)
sortBed(bed_in_path=OUT_DIR+'/f3_snv_filtered.bed', bed_out_path=OUT_DIR+'/f4_snv_sorted.bed', bedtools_path=bedtools_path)
dupletCluster(bed_in_path=OUT_DIR+'/f4_snv_sorted.bed', bed_out_path=OUT_DIR+'/f5_snv_duplet_site.bed', cluster_distance=CLUSTER_DISTANCE,cluster_size=2)
getClusterBed(bed_in_path=OUT_DIR+'/f5_snv_duplet_site.bed', bed_out_path=OUT_DIR+'/f6_snv_duplet_cluster.bed')
bed2addfa(bed_in_path=OUT_DIR+'/f5_snv_duplet_site.bed', ref_in_path=ref_in_path, fa_out_path=OUT_DIR+'/f7_snv_duplet_site.a20.fa', AROUND=SITE_ADD_AROUND, SEP='q')
bed2flankfa(bed_in_path=OUT_DIR+'/f6_snv_duplet_cluster.bed', ref_in_path=ref_in_path, fa_out_path=OUT_DIR+'/f8_snv_duplet_cluster.f20000.fa', AROUND=CLUSTER_FLANK_AROUND, SEP='r')
blatAlign(ref_in_path=OUT_DIR+'/f8_snv_duplet_cluster.f20000.fa', query_in_path=OUT_DIR+'/f7_snv_duplet_site.a20.fa', out_dir=OUT_DIR+'/f9_blat_out', blat_path=blat_path, CPU=CPU)



#----------------------------------------------------------------------------
print('Finished !\n')
print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


