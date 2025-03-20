import os
import json
import pandas as pd
import glob

directory=sys.argv[1]

files = glob.glob(f"{directory}/seed*sample*/summary_confidences.json")
for f in files :
    with open(f, 'r') as my_json : 
        info = json.load(my_json)
        pae_min = info["chain_pair_pae_min"][0][1]
        pae.append(pae_min)
f.close()

df = pd.DataFrame({
    'pae_min' : pae, 
})

sample = os.path.basename(directory)
df['sample'] = sample

df.to_csv(f"{os.path.basename(directory)}.pae_min.tsv", sep = "\t", header = True, index = False)