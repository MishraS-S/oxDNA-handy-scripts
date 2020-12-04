# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 22:31:31 2020
@author: Shubham
"""
__descr__ = """
Input arguments in order - input json file, input scaffold sequence file,
starting vhelix position and starting column position of scaffold if circular

Doesn't account for loops, only skips. Sequence file should be nice. There should one continuous scaffold
"""
import json
import argparse
from pathlib import Path


def proc_input():
    parser = argparse.ArgumentParser(
        description=__descr__,
    )
    parser.add_argument("-j", "--json",
                        help="input json file",
                        type=str
                        )
    parser.add_argument("-s", "--sequence",
                        help="input scaffold sequence file",
                        type=str
                        )
    parser.add_argument("-v", "--vHelix",
                        help="virtual helix number where the scaffold sequence begins (count the number of rows from the top in cadnano)",
                        type=int,
                        default=-1  # if linear scaffold, this default value will find the start value automatically
                        )
    parser.add_argument("-c", "--column",
                        help="column number where the scaffold sequence begins (count the number of little squares from the left in cadnano, no matter if they are empy or not)",
                        type=int,
                        default=-1  # if linear scaffold, this default value will find the start value automatically
                        )
    parser.add_argument("-o", "--out_name",
                        help="output file name",
                        type= str,
                        default= "caca.sqs"
                        )
    args = parser.parse_args()
    return args


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


def load_data(args):
    def compute_start(json):
        # processing json for easier access via vhelix id
        ori = {}
        preindex = 0
        for helix in json_input['vstrands']:
            ori[helix['num']] = helix
            ori[helix['num']]['preindex'] = preindex  # associate vhelix with index
            preindex += 1
        for vhs in range(len(ori)):
            for cs in range(len(ori[helix['num']])):
                if ori[vhs]['scaf'][cs][0:2] == [-1, -1] and \
                        ori[vhs]['scaf'][cs][2:4] != [-1, -1]:  # if it is a start site
                    return vhs, cs
        return [0,5]  # TODO: test the starting position for linear scaffold

    with open(Path(args.json), "r") as json_file:
        json_input = json.load(json_file)

    with open(Path(args.sequence), "r") as seq_file:
        seq = seq_file.read()
        seq = seq.upper()
        seq = list(seq)
        for index, base in enumerate(seq):
            if base.isspace():
                raise ValueError('Sequence contains whitespace ' \
                                 + repr(str(base)) + ' at position ' + str(index) + '!')

    # start sites
    if args.vHelix != -1 and args.column != -1:
        vhs = args.vHelix
        cs = args.column
    else:
        vhs, cs = compute_start(args.json)

    return [json_input, seq, vhs, cs]


def compute_output(data):  # TODO: make this compatible with circular scaffold
    json_input, seq, vhs, cs = data
    # processing json for easier access via vhelix id
    ori = {}
    preindex = 0
    for helix in json_input['vstrands']:
        ori[helix['num']] = helix
        ori[helix['num']]['preindex'] = preindex  # associate vhelix with index
        preindex += 1

    # initialize empty output matrix populated with some default base for eg. A
    m = len(ori)
    n = len(ori[0]['scaf'])

    mxn = []
    for x in range(m):
        sub = []
        for x in range(n):
            sub.append('R')
        mxn.append(sub)

    skiplist = []  # list of skip sites
    print(ori[1]['scaf'])  # tim edit
    if ori[vhs]['scaf'][cs][0:2] == [-1, -1] and \
            ori[vhs]['scaf'][cs][2:4] != [-1, -1]:  # if it is a start site
        i = 0  # sequence counter
        vh = vhs  # vhelix id
        c = cs  # column id

        while ori[vh]['scaf'][c][2:4] != [-1, -1]:  # if it not an end site

            vhref = vh
            preindex = ori[vh]['preindex']
            if ori[vh]['skip'][c] == 0:  # if current position is not a skip
                if (vh % 2) == 1:
                    mxn[preindex][c] = seq[i]
                else:
                    mxn[preindex][c] = seq[i]  # Tim edit, changed from  compl(seq[i])
                i += 1
            else:
                skiplist.append([preindex, c]) #  Tim: isn't there i += 1 missing?
            vh = ori[vh]['scaf'][c][2]
            c = ori[vhref]['scaf'][c][3]
        # since while loop stops just before the last nucleotide
        preindex = ori[vh]['preindex']
        if vh % 2 == 1:
            mxn[preindex][c] = seq[i]
        else:
            mxn[preindex][c] = seq[i] # Tim edit, changed from compl(seq[i])
    else:
        raise ValueError('Not a scaffold start site!')

    # remove skip sites from mxn
    for el in skiplist:
        i = el[0]
        j = el[1]
        del mxn[i][j]

    return mxn


def write_sequenceFile(name, mxn):
    with open(name, "w+") as out:
        i = 0
        for row in mxn:
            out.write(''.join(map(str, row)))
            if i != (len(mxn) - 1):
                out.write('\n')
            i += 1


def main(args):
    data = load_data(args)
    out = compute_output(data)
    write_sequenceFile(args.out_name, out)


if __name__ == "__main__":
    args = proc_input()
    print(type(args))
    #main(args)
