#!/usr/bin/env python3

import json
from itertools import product
import argparse
import os
import glob

def read_fasta(fasta) : 

    fa_dict = {}
    seq = ''
    header = False
    with open(fasta, 'r') as f : 
        for line in f : 
            if line.startswith(">") : 
                if header : 
                    fa_dict[header] = seq

                    header = line.strip().replace(">", "")
                    seq = ''
                else : 
                    header = line.strip().replace(">", "")
            else : 
                seq += line.strip()
        else : 
            fa_dict[header] = seq
    f.close()

    return fa_dict

def fa_to_json(fa1, fa2) : 

    fa1_dict = read_fasta(fa1)
    fa2_dict = read_fasta(fa2)
    combined_fa_dict = {**fa1_dict, **fa2_dict}

    fa1_names = list(fa1_dict.keys())
    fa2_names = list(fa2_dict.keys())

    cross_product = [ list(pair) for pair in product(fa1_names, fa2_names)]
    filt_cross_prod = [pair for pair in cross_product if len(set(pair)) > 1]

    for i in filt_cross_prod : 
        A_name = i[0]
        B_name = i[1]

        output = f"{A_name}_{B_name}.json"

        A_seq = combined_fa_dict[A_name]
        B_seq = combined_fa_dict[B_name]

        af_input = {
            "name" : A_name, 
            "sequences" : [
                {
                    "protein" : {
                        "id" : "A",
                        "sequence" : A_seq
                    }
                },
                {
                    "protein" : {
                        "id" : "B",
                        "sequence" : B_seq
                    }
                }
            ],
            "modelSeeds" : [
                1
            ], 
            "dialect" : "alphafold3",
            "version" : 1
        }

        with open(output, 'w') as json_file : 
            json.dump(af_input, json_file, indent = 4)
        json_file.close()

def get_args() : 

    """Parse command line parameters from input"""
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-fa1", type = str, required = True)
    parser.add_argument("-fa2", type = str, required = True)

    return parser.parse_args()

def main() : 

    args = get_args()
    fa_to_json(args.fa1, args.fa2)

if __name__ == "__main__" : 

    main()



