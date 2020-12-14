# -*- coding: utf-8 -*-
_author__ = "Tim Peisker"
__copyright__ = "Copyright 2020, Dietzlab (TUM)"
__license__ = "None"
__maintainer__ = "Tim Peisker"
__email__ = "tim.peisker@tum.de"
__status__ = "Development"
__descr__ = """
Based on the original script from https://github.com/Mishrito/oxDNA-handy-scripts

Input arguments in order - input json file, input scaffold sequence file,
starting vhelix position and starting column position of scaffold if circular

Doesn't account for loops, only skips. Sequence file should be nice. There should one continuous scaffold

Note from Tim Peisker: Cadnano saves DNA origami structures as .json files. These .json files store each virtual Helix
from Cadnano as a separate helix entry. A Cadnano virtual Helix consists of two rows with little squares (one for
the scaffold and one for the staples). Each square as represented in the .json file has a list attribute of 4 values.
The values are [5'vHelix neighbor, 5'column neighbor, 3'vHelix neighbor, 3'column neighbor]. If there is no connection
in the 3' or 5' direction, the respective values will be -1.
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


def compute_start(json_path):
    with open(Path(json_path), "r") as json_file:
        json_input = json.load(json_file)
    # processing json for easier access via vhelix id
    ori = {}
    preindex = 0
    for helix in json_input['vstrands']:
        ori[helix['num']] = helix
        ori[helix['num']]['preindex'] = preindex  # associate vhelix with index
        preindex += 1
    for vhs in range(len(ori)):
        for cs in range(len(ori[helix['num']]['scaf'])):
            if ori[vhs]['scaf'][cs][0:2] == [-1, -1] and \
                    ori[vhs]['scaf'][cs][2:4] != [-1, -1]:  # if it is a start site
                return [vhs, cs]
    return [-1, -1]  # TODO: test the starting position for linear scaffold


def load_data(args):

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


def compute_output(data):
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
    if cs != -1 and vhs != -1:  # if it is a start site TODO test if this works with a false starting site
        i = 0  # sequence counter
        vh = vhs  # vhelix id
        c = cs  # column id

        while ori[vh]['scaf'][c][2:4] != [-1, -1] and ori[vh]['scaf'][c][2:4] != [vhs, cs]:  # if it not an end site or back to the origin

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
