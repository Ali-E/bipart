import pandas as pd
import sys
import os


os.system('rm filenames.txt')
os.system('sh src_ext_res.sh')

input_file = sys.argv[1]

ext_num = []
df = pd.DataFrame(columns=['#T', '-RT ln Z', 'Z', 'idx'])

pairs = set()
temp_number = 0
names = []

with open('filenames.txt') as f:
    names = f.readlines()


for i, name in enumerate(names):
    if name[:len(input_file)] != input_file:
        continue

    cname = name[:-4]
    files = cname.split("-")
    if len(files) < 2 or files[0][len(input_file):] == files[1][len(input_file):]:
        continue

    pair_number = int(files[0][len(input_file):])//2
    row = pd.read_csv(name[:-1], delimiter='\t')
    row['idx'] = pair_number
    df = df.append(row)


df = df.sort_values(by=['idx'])
df = df.set_index('idx')

df.to_csv("pirna_results.tsv", sep='\t')

