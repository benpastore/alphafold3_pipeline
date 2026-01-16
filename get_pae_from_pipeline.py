#!/usr/bin/env python3

file = "/fs/scratch/PCON0160/00AA_AF3/ad/2704871dacfa16d7a3715ae8e420f2/.command.sh"

lines = "pae_min\tsample\n"
with open(file, 'r') as f : 
    for line in f : 
        if line.startswith("cat") : 
            info = line.strip().split(" ") 
            for i in info : 
                if "/fs/scratch/PCON0160" in i : 
                    with open(i, 'r') as k : 
                        for kline in k : 
                            if not kline.startswith("pae_min") : 
                                kinfo = kline.strip().split("\t")
                                lines += f"{kinfo[0]}\t{kinfo[1]}\n"
                    k.close()
f.close() 

op =  open("prg1_hits_from_screen_against_proteome.tsv", 'w')
op.write(lines)
op.close()

