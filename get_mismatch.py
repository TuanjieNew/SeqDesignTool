#!/usr/bin/env python
#fn; get_mismatch.py

#ACTGCAGCGTCATAGTTTTTGAG
import os
import copy



def getMismatch(start,seq,name,end):
    #name = seq
    quality = 'IIIIIIIIIIIIIIIIIIIIII'
    OUTFILE = open('./mis_test.fastq','a')
    ls = list(seq)
    ls_1 = copy.deepcopy(ls)
    ii = start+1
    for i in ls_1[ii:end]:
        if i == 'A':
            ls_1[ii] = 'T'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'G'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'C'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

        if i == 'T':
            ls_1[ii] = 'A'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'G'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'C'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)
        if i == 'G':
            ls_1[ii] = 'T'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'A'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'C'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

        if i == 'C':
            ls_1[ii] = 'T'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'G'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)

            ls_1[ii] = 'A'
            ls_1[21] = 'G'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            ls_1[22] = 'A'
            OUTFILE.write('@'+name+'\n'+''.join(ls_1)+'\n'+'+'+name+'\n'+quality+'\n')
            getMismatch(ii,ls_1,name,end - 1)
        ii+=1


seq = 'GCTGCGTCGTCGTAGTTTTTTGG'
getMismatch(-1, seq, seq, 21)
