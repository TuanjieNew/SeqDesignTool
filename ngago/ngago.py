#!/usr/bin/env python
# fn: ngago.py

import os 
import copy
li = os.listdir('.')
if 'ngago_off.fastq' in li:
    os.remove('ngago_off.fastq')
    #os.mkdir('./ngago_off')

#seq = 'CCAAAGTCCAAGGTATTAGCAGCC'
#seq = 'GTCGCTCGTCGGAGCTGCAGGGAC'
#seq = 'TGCGTCGTCGTAGTTTTTTGGGGG'

def getMisSeq(seq):
    OUTFILE = open('./ngago_off.fastq','a')
    quality = 'IIIIIIIIIIIIIIIIIIIIIIII'
    seqlen = len(seq)
    seqls = list(seq)
    seq_old = copy.deepcopy(seqls)
    #print(seqlen)
    #print(seqls)
    for i in range(seqlen):
        seqls = copy.deepcopy(seq_old)
        #print('seq_old: '+''.join(seq_old))
        if seqls[i] == 'A':
            seqls[i] = 'G'
            seq2 = ('').join(seqls)
            #print(seq2)
            
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'T'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'C'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
        elif seqls[i] == 'T':
            seqls[i] = 'G'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'A'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'C'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
        elif seqls[i] == 'G':
            seqls[i] = 'A'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'T'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'C'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
        elif seqls[i] == 'C':
            seqls[i] = 'A'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'T'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')
            seqls[i] = 'G'
            seq2 = ('').join(seqls)
            OUTFILE.write('@'+seq+'\n'+seq2+'\n'+'+'+seq+'\n'+quality+'\n')

    OUTFILE.close()
def getMd(ln_ls):
    count = 0
    for i in ln_ls:
        count += 1
        if i[:2] == 'MD':
            return i
def getMisLoc(seq, misSeq):
    i = 0
    mis_ls = []
    for j in range(len(seq)):
        if seq[j] == misSeq[j]:
            continue
        else:
            mis_ls.append(j)
        #i += 1
    return mis_ls

def getScore(mis_loc):
    M = [ 0.19736842,  0.19407895,  0.25328947,  0.26315789,  0.19078947,0.22697368,  0.25328947,  0.        ,  0.09539474,  0.14802632,0.11842105,  0.25      ,  0.24671053,  0.26644737,  0.16776316,0.21052632,  0.23684211,  0.20394737,  0.21381579,  0.24342105,0.22368421,  0.20394737,  0.1875    ,  0.17434211]
    score = 0
    for i in mis_loc:
        score += M[i]
    return score

q_seq = 'GTCGCTCGTCGGAGCTGCAGGGACCGGCGCGAGCGAGTGCTGGACTGTTTGTGCAGGGCTCCGAGGGGACCCATGTGGCTCAGGGTGGCTAAGGGGGCAATGCTGCGTCGTCGTAGTTTTTTGGGGG'
GFILE = open('ngago_gseq.txt','w')
# get mismatch sequence
dic = {}
for i in range(len(q_seq)):
    if i < len(q_seq) - 23:
        seq = q_seq[i:i+24]
        GFILE.write(seq+'\n')
        dic[seq] = ''
        getMisSeq(seq)
GFILE.close()
# for fasta file
#os.system('bowtie2 -x /1/public_resources/hg19/hg19 -p 2 -r -N 0 -U ./ngago_off/ngago_off.fastq --no-unal --no-hd -S ./ngago_sam/ngago_off.sam')
# for fastq file
os.system('bowtie2 -x /1/public_resources/hg19/hg19 -p 2 -N 0 -U ./ngago_off.fastq --no-unal --no-hd -S ./ngago_off.sam')
SAMFILE = open('./ngago_off.sam','r')
lnum = 0
for ln in SAMFILE:
    lnum += 1
    ln = ln.strip('\n')
    ln_ls = ln.split('\t')
    if len(ln_ls) == 19:
        md = ln_ls[17]
    elif len(ln_ls) == 20:
        md = ln_ls[18]
    else:
        md = getMd(ln_ls)
    mis_loc=[]
    if md == 'MD:Z:24':
        mis_loc = getMisLoc(ln_ls[0], ln_ls[9])
    else:
        continue
    score = 0
    score = getScore(mis_loc)
    chr_pos = ''
    if int(ln_ls[1]) == 0:
        chr_pos = ln_ls[2]+':+'+str(ln_ls[3])
    elif int(ln_ls[1]) == 16:
        chr_pos = ln_ls[2] + ':-'+str(ln_ls[3])
    #dic[ln_ls[0]].append({ln_ls[9]+':'+str(mis_loc):score}) 
    dic[ln_ls[0]]=dic[ln_ls[0]]+ln_ls[9]+'\t'+str(score)+'\t'+str(mis_loc)+'\t'+chr_pos + ';'
print(len(dic))
SAMFILE.close()
dic_2 = {}
for e in dic:
    score_sum = 100
    #if len(dic[e]) == 0:
        #continue
    dic_ls = dic[e].split(';')
    key_str = e+';'
    for i in dic_ls:
        i_ls = i.split('\t')
        if len(i_ls) < 3:
            continue
        #print(i_ls)
        score_sum -= float(i_ls[1])*100
        key_str = key_str+i+';'
    dic_2[key_str] = score_sum
dic_2_sort = sorted(dic_2.iteritems(), key = lambda asd:asd[1], reverse = True)
RESULT = open('ngago_result.txt','w')
for i in dic_2_sort:
    info_ls = i[0].split(';')
    RESULT.write('Guide Sequence: '+info_ls[0]+'\n')
    RESULT.write('Guide Score: '+str(i[1])+'\n')
    if float(i[1]) == 100.0:
        RESULT.write('----------------------------------------\n')
        continue
    RESULT.write('Sequence\tScore\tMismatches\tLocus\n')
    for i in info_ls[1:]:
        #off_ls = i.split(':')
        RESULT.write(i+'\n')
    RESULT.write('----------------------------------------\n')
RESULT.close()



'''
GFILE = open('ngago_gseq.txt','r')
for line in GFILE:
    line = line.strip('\n')

    SAMFILE = open('./ngago_sam/'+line, 'r')
    RESULT = open('./ngago_result/'+line,'w')
    RESULT.write('Guide Sequence: '+line+'\n')
    #RESULT.write('Guide Score\tSequence\tScore\n')
    dic = {}
    score_sum = 1
    for ln in SAMFILE:
        ln = ln.strip('\n')
        ln_ls = ln.split('\t')
        if len(ln_ls) == 19:
            md = ln_ls[17]
        elif len(ln_ls) == 20:
            md = ln_ls[18]
        else:
            md = getMd(ln_ls)
        mis_loc=[]
        if md == 'MD:Z:24':
            mis_loc = getMisLoc(line, ln_ls[9])
        score = 0
        score = getScore(mis_loc)
        score_sum -= socre
        dic[ln+':'+str(mis_loc)] = score
    dic_sort = sorted(dic.iteritems(), key = lambda asd:asd[1], reverse = True)
    RESULT.write('Guide score: '+str(score_sum)+'\n')
    RESULT.write('-----------------------------------\n')
    RESULT.write('Sequence\tScore\tMismatch Location\n')
    for e in dic_sort:
        sequnence = e[0].split(':')[0]
        seq_mis_loc= e[0].split(':')[1]
        seq_socre = e[1]
        RESULT.write(sequence+'\t'+str(seq_score)+'\t'+str(seq_mis_loc)+'\n')
    RESULT.close()
    SAMFILE.close()

GFILE.close()
'''
