#!/usr/bin/env python
#fn: crispr.py
# produced by niutj, date: 20160721

import os
import copy


def getOff(start,seq, name):
    # one mismatch 
    if name == 'CTCATCTGAGCTGCAGGGACCGG':
        print('off 1')
    seq = ''.join(seq) 
    OUTFILE = open('./crispr_off.fastq','a')
    #FAFILE = open('./crispr_off.fasta','a')
    quality = 'IIIIIIIIIIIIIIIIIIIIIII'
    seq_ls = list(seq)
    old_seq = copy.deepcopy(seq_ls)
    i = start + 1
    if seq_ls == list('ACTGCAGCGTCGTAGTTTTTGAG'):
        print('hello')
    
    for e in old_seq[start+1:21]:
        seq_ls = copy.deepcopy(old_seq)
        
        if e == 'A':
            seq_ls[i] = 'T'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'G'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'C'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

        elif e == 'T':
            seq_ls[i] = 'A'
            seq_ls[21] = 'G'
            #seq2 = ''.join(seq_ls)
            #print(type(seq2))
            #print(type(seq))
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'G'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            
            seq_ls[i] = 'C'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

        elif e == 'G':
            seq_ls[i] = 'T'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'A'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'C'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
        
        elif e == 'C':
            seq_ls[i] = 'T'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'A'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

            seq_ls[i] = 'G'
            seq_ls[21] = 'G'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')
            seq_ls[21] = 'A'
            #FAFILE.write(''.join(seq_ls)+'\n')
            OUTFILE.write('@'+name+'\n'+''.join(seq_ls)+'\n'+'+'+name+'\n'+quality+'\n')

        i+=1
    OUTFILE.close()
    # two mismatch
def getOff_2(start, seq, name):
    if name == 'CTCATCTGAGCTGCAGGGACCGG':
        print('off 2')
    i = start + 1
    seq_ls = list(seq)
    old_seq = copy.deepcopy(seq_ls)
    for e in old_seq[start+1:20]:
        seq_ls = copy.deepcopy(old_seq)
        if seq_ls[i] == 'A':
            seq_ls[i] = 'T'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'G'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'C'
            getOff(i, seq_ls, name)
        elif seq_ls[i] == 'T':
            seq_ls[i] = 'A'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'G'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'C'
            getOff(i, seq_ls, name)
        elif seq_ls[i] == 'G':
            seq_ls[i] = 'T'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'A'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'C'
            getOff(i, seq_ls, name)
        elif seq_ls[i] == 'C':
            seq_ls[i] = 'T'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'G'
            getOff(i, seq_ls, name)
            seq_ls[i] = 'A'
            getOff(i, seq_ls, name)
        i += 1
def getOff_3(start, seq, name):
    if name == 'CTCATCTGAGCTGCAGGGACCGG':
        print('off 3')
    i = start + 1
    seq_ls = list(seq)
    old_seq = copy.deepcopy(seq_ls)
    for e in old_seq[start+1:19]:
        seq_ls = copy.deepcopy(old_seq)
        if seq_ls[i] == 'A':
            seq_ls[i] = 'T'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'G'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'C'
            getOff_2(i, seq_ls, name)

        elif seq_ls[i] == 'T':
            seq_ls[i] = 'A'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'G'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'C'
            getOff_2(i, seq_ls, name)

        elif seq_ls[i] == 'G':
            seq_ls[i] = 'A'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'T'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'C'
            getOff_2(i, seq_ls, name)

        elif seq_ls[i] == 'C':
            seq_ls[i] = 'A'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'G'
            getOff_2(i, seq_ls, name)
            seq_ls[i] == 'T'
            getOff_2(i, seq_ls, name)

        i += 1

def getOff_4(seq,name):
    if name == 'CTCATCTGAGCTGCAGGGACCGG':
        print('off 4')
    i = 0
    seq_ls = list(seq)
    old_seq = copy.deepcopy(seq_ls)
    for e in old_seq[0:18]:
        seq_ls = copy.deepcopy(old_seq) 
        #print(seq_ls)
        if seq_ls[i] == 'A':
            seq_ls[i] = 'T'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'G'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'C'
            getOff_3(i, seq_ls, name)
        elif seq_ls[i] == 'T':
            seq_ls[i] = 'A'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'G'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'C'
            getOff_3(i, seq_ls, name)
        elif seq_ls[i] == 'G':
            seq_ls[i] = 'T'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'A'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'C'
            getOff_3(i, seq_ls, name)
        elif seq_ls[i] == 'C':
            seq_ls[i] = 'T'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'G'
            getOff_3(i, seq_ls, name)
            seq_ls[i] == 'A'
            getOff_3(i, seq_ls, name)
        i += 1
    #print(seq_ls)

def reverse(seq):
    seq = seq[::-1]
    rseq = ''
    for i in seq:
        if i == 'A':
            rseq += 'T'
        elif i == 'T':
            rseq += 'A'
        elif i == 'G':
            rseq += 'C'
        elif i == 'C':
            rseq += 'G'
    #print(rseq)
    return rseq
def getMd(ln_ls):
    count = 0
    for i in ln_ls:
        count += 1
        if i[:2] == 'MD':
            return i 
def getMisLoc(seq, misSeq):
    i = 0
    mis_ls = []
    for j in range(len(seq)-3):
        if seq[j] == misSeq[j]:
            continue
        else:
            mis_ls.append(j)
    return mis_ls

def getScore(mis_pos):
    M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    mis_len = len(mis_pos)
    dis = 0.0
    misScore = 1
    for i in range(mis_len):
        #print(mis_pos[i])
        misScore = misScore * (1 - M[mis_pos[i]])
        j = i + 1
        while j < mis_len:
            dis += abs(mis_pos[i] - mis_pos[j])-1
            j += 1
    mean_dis = dis/mis_len
    score = misScore * (1/((19 - mean_dis)*4/19 + 1 ))/(mis_len*mis_len)
    return score*100


# rm crispr_off.fastq 
li = os.listdir('.')
if 'crispr_off.fastq' in li:
    print('crispr_off.fastq')
    os.system('rm crispr_off.fastq')

a = os.listdir('./offTarget')
for fl in a:
    os.remove('./offTarget/'+fl)
FILE = open('crispr_seq.txt','r')
OUTFILE = open('guide_seq.txt', 'w')

dic = {}
for line in FILE:
    line = line.strip('\n')
    ln_len = len(line)
    print('line length: '+str(ln_len))
    i = 0
    while i < ln_len - 22:
        seq = line[i:i+23]
        if seq[21:] == 'GG':
            dic[seq] = ''
            getOff(-1, seq,seq)
            getOff_2(-1, seq,seq)
            getOff_3(-1, seq,seq)
            getOff_4(seq,seq)
        elif seq[0:2] == 'CC':
            seq = reverse(seq)
            getOff(-1, seq,seq)
            getOff_2(-1, seq,seq)
            getOff_3(-1, seq,seq)
            getOff_4(seq,seq)
            dic[seq] = '' 
        #elif seq[21:] == 'AA':
            #print(seq)
        i += 1
    for i in dic:
        OUTFILE.write(dic[i] + '\n')
    print('dic length: '+str(len(dic)))
OUTFILE.close()
FILE.close()

os.system('bowtie2 -x /1/public_resources/hg19/hg19 -p 20 -N 0 -U ./crispr_off.fastq --no-unal --no-hd -S ./crispr_off.sam')

#os.system('bwa aln -t 20 /1/niutj/ref/bwa_hg19/hg19bwaidx ./crispr_off.fastq > crispr_off.bwa')
#os.system('bwa samse  /1/niutj/ref/bwa_hg19/hg19bwaidx ./crispr_off.bwa ./crispr_off.fastq > crispr_off.sam')

SAMFILE = open('./crispr_off.sam','r')
dic_sam = {}
lnum = 0
for ln in SAMFILE:
    ln = ln.strip('\n')
    ln_ls = ln.split('\t')
    if len(ln_ls) == 19:
        md = ln_ls[17]
    elif len(ln_ls) == 20:
        md = ln_ls[18]
    else: 
        md = getMd(ln_ls)
    if md == 'MD:Z:23':
        lnum += 1
        dic_sam[ln_ls[9]] = ln_ls
SAMFILE.close()
#print(dic_sam)
print(len(dic_sam))
print('lnum: '+str(lnum))

lnum = 0
for e in dic_sam:
    '''
    lnum += 1
    ln = ln.strip('\n')
    '''
    ln_ls = dic_sam[e]
    if len(ln_ls) == 19:
        md = ln_ls[17]
    elif len(ln_ls) == 20:
        md = ln_ls[18]
    else:
        md = getMd(ln_ls)
    mis_loc = {}
    if md == 'MD:Z:23':
        if int(ln_ls[1]) == 0:
            mis_loc = getMisLoc(ln_ls[0], ln_ls[9])
        else:
            mis_loc = getMisLoc(ln_ls[0], reverse(ln_ls[9]))
    else:
        continue
    score = 0
    if len(mis_loc) == 0:
        continue
    score = getScore(mis_loc)
    chr_pos = ''
    if int(ln_ls[1]) == 0:
        chr_pos = ln_ls[2] + ':+' + str(ln_ls[3])
    elif int(ln_ls[1]) == 16:
        chr_pos = ln_ls[2] + ':-' + str(ln_ls[3])
    dic[ln_ls[0]] = dic[ln_ls[0]] + ln_ls[9] + '\t'+str(score) + '\t' + str(mis_loc) + '\t' + chr_pos + ';'

print(len(dic))
#SAMFILE.close()
dic_2 = {}
for e in dic:
    score_sum = 0
    dic_ls = dic[e].split(';')
    key_str = e+';'
    for i in dic_ls:
        i_ls = i.split('\t')
        if len(i_ls) < 3:
            continue
        score_sum += float(i_ls[1])
        key_str = key_str+i+';'
    dic_2[key_str] = 100/(100+score_sum)
dic_2_sort = sorted(dic_2.iteritems(), key = lambda asd:asd[1], reverse = True)
RESULT = open('crispr_result.txt', 'w')
for i in dic_2_sort:
    info_ls = i[0].split(';')
    RESULT.write('Guide Sequence: '+info_ls[0]+'\n')
    RESULT.write('Guide Score: '+str(i[1])+'\n')
    if float(i[1]) == 100.0:
        RESULT.write('--------------------------------------------\n')
        continue
    RESULT.write('Sequence\tScore\tMismatches\tLocus\n')
    for i in info_ls[1:]:
        RESULT.write(i+'\n')
    RESULT.write('--------------------------------------------\n')
RESULT.close()
    
