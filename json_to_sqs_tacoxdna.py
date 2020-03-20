# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 22:31:31 2020

@author: Shubham
"""

################
# Input arguments in order - input json file, input scaffold sequence file,
# starting vhelix position and starting column position of scaffold
#
# Doesn't account for loops, only skips. Sequence file should be nice
################

import json
import sys

def compl(base):
    if base == 'A':
        base = 'T'
    elif base == 'T':
        base = 'A'
    elif base == 'G':
        base = 'C'
    elif base == 'C':
        base = 'G'
    return base


with open(sys.argv[1], "r") as json_file:
    json_input=json.load(json_file)

with open(sys.argv[2], "r") as seq_file:
    seq=seq_file.read()
    seq=seq.upper()
    seq=list(seq)
    for index,base in enumerate(seq):
        if base.isspace():
            raise ValueError('Sequence contains whitespace ' \
	    +repr(str(base))+' at position '+str(index)+'!')

#start sites
vhs = int(sys.argv[3])
cs = int(sys.argv[4])

#processing json for easier access via vhelix id
ori = {}
preindex = 0
for helix in json_input['vstrands']:
        ori[helix['num']] = helix
        ori[helix['num']]['preindex'] = preindex #associate vhelix with index
        preindex += 1


#initialize empty output matrix populated with some default base for eg. A
m = len(ori)
n = len(ori[0]['scaf'])

mxn = []
for x in range(m):
     sub = []
     for x in range(n):
         sub.append('A')
     mxn.append(sub)

 
    
skiplist = [] #list of skip sites

if ori[vhs]['scaf'][cs][0:2] == [-1,-1] and \
   ori[vhs]['scaf'][cs][2:4] != [-1,-1]: #if it is a start site
       i=0 #sequence counter
       vh = vhs #vhelix id
       c= cs #column id

       while ori[vh]['scaf'][c][2:4] != [-1,-1]: #if it not an end site

           vhref = vh
           preindex = ori[vh]['preindex']
           if ori[vh]['skip'][c] == 0: #if current position is not a skip
               if (vh%2) == 1:
                   mxn[preindex][c] = seq[i]
               else :
                   mxn[preindex][c] = compl(seq[i])
               i += 1
           else :
              skiplist.append([preindex,c])
           vh = ori[vh]['scaf'][c][2]
           c = ori[vhref]['scaf'][c][3]
       #since while loop stops just before the last nucleotide  
       preindex = ori[vh]['preindex']
       if vh%2 == 1:
            mxn[preindex][c] = seq[i]
       else :
            mxn[preindex][c] = compl(seq[i])
       
else : 
    raise ValueError('Not a scaffold start site!')
    
#remove skip sites from mxn
for el in skiplist:
    i = el[0]
    j = el[1]
    del mxn[i][j]
    
        
with open("caca.sqs","w+") as out:
    i=0
    for row in mxn:
        out.write(''.join(map(str, row)))
        if i != (len(mxn)-1):
            out.write('\n')
        i += 1
