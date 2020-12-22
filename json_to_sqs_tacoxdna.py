# -*- coding: utf-8 -*-

__descr__ = """
Input arguments in order:
-input json file
-input scaffold sequence file
-starting vhelix position of scaffold if circular
-starting column position of scaffold if circular
-output_file name

Doesn't account for loops, only skips. Sequence file should be nice. There should one continuous scaffold

To find the vhelix and column position values of your start site:
For linear scaffolds you don't need to do anything. The default values of -1 will tell the script to find it
automatically. Just make sure that the beginning of your sequence file matches the 5' starting site. 
If, however, your scaffold is circular, you need to tell the script the beginning of your sequence. To do that, 
open the design in Cadnano and locate the little square that is your start site. For vhelix enter the id found in the
orange circle to the very left of the corresponding row of little squares. For column enter the number of squares your
starting position is away from the leftmost square (including empty squares), e.g. if your sequence starts at the
leftmost square, it will be column = 0.
"""

"""
Cadnano saves the 5' and 3' neighbors of each square in the following structure:
[5'vHelix neighbor, 5'column neighbor, 3'vHelix neighbor, 3'column neighbor]
If there is no connection in the 3' or 5' direction, the respective values will be -1.
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
                        type=str,
                        default="caca.sqs"
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
            has_5p_neighbor = ori[vhs]['scaf'][cs][0:2] != [-1, -1]
            has_3p_neighbor = ori[vhs]['scaf'][cs][2:4] != [-1, -1]
            is_start_site = not has_5p_neighbor and has_3p_neighbor
            if is_start_site:
                return [vhs, cs]
    return [-1, -1]


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
    if cs != -1 and vhs != -1:  # if it is a start site TODO test if this works with a false starting site
        i = 0  # sequence counter
        vh = vhs  # vhelix id
        c = cs  # column id

        while ori[vh]['scaf'][c][2:4] != [-1, -1] and ori[vh]['scaf'][c][2:4] != [vhs,
                                                                                  cs]:  # if it not an end site or back to the origin

            vhref = vh
            preindex = ori[vh]['preindex']
            if ori[vh]['skip'][c] == 0:  # if current position is not a skip
                if (vh % 2) == 1:
                    mxn[preindex][c] = seq[i]
                else:
                    mxn[preindex][c] = seq[i]  # Tim edit, changed from  compl(seq[i])
                i += 1
            else:
                skiplist.append([preindex, c])  # Tim: isn't there i += 1 missing?
            vh = ori[vh]['scaf'][c][2]
            c = ori[vhref]['scaf'][c][3]
        # since while loop stops just before the last nucleotide
        preindex = ori[vh]['preindex']
        if vh % 2 == 1:
            mxn[preindex][c] = seq[i]
        else:
            mxn[preindex][c] = seq[i]  # Tim edit, changed from compl(seq[i])
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


def main():
    args = proc_input()
    data = load_data(args)
    out = compute_output(data)
    write_sequenceFile(args.out_name, out)


if __name__ == "__main__":
    main()
