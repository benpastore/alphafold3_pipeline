import os
import json
import glob
import sys

directory=sys.argv[1]

files = glob.glob(f"{directory}/seed*sample*/summary_confidences.json")
pae = []
for f in files :
    with open(f, 'r') as my_json : 
        info = json.load(my_json)
        pae_min = info["chain_pair_pae_min"][0][1]
        pae.append(pae_min)
    my_json.close()

sample = sys.argv[2]
lines = 'pae_min\tsample\n'
for p in pae : 
    lines += f"{p}\t{sample}\n"

output=f"{sample}.pae_min.tsv"
op = open(output, 'w')
op.write(lines)
op.close()