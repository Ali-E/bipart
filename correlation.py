import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import math
import sys
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr



plt.switch_backend('agg')

df_bp = pd.read_csv("pre_computed/bppart_37.csv", delimiter='\t')
df_pi = pd.read_csv("pre_computed/table_37.csv")

df_len = pd.read_csv(sys.argv[1], header=None, names=["len1", "len2"])
df_pi["len1"] = df_len["len1"]
df_pi["len2"] = df_len["len2"]

row = df_pi.iloc[32076]
if row["-RT ln Z"] == -np.inf:
    print("ccccc")

idx_to_del = []
for idx, row in df_pi.iterrows():
    if row["-RT ln Z"] == -np.inf or row["Z"] == np.inf or math.isnan(row["-RT ln Z"]) or \
      math.isnan(row["Z"]) or np.isnan(row["Z"]) or np.isnan(row["-RT ln Z"]):
        idx_to_del.append(idx)


df_pos = pd.read_csv("pre_computed/bpmax_37.csv", delimiter='\t')
df_pos = df_pos[:len(df_pi)]

df_bp = df_bp[:len(df_pi)]


df_pi = df_pi.drop(df_pi.index[idx_to_del])
df_pos = df_pos.drop(df_pos.index[idx_to_del])
df_bp = df_bp.drop(df_bp.index[idx_to_del])


x = df_bp["-log(QI)"]/(df_pi["len1"]+df_pi["len2"])
y = df_pi["-RT ln Z"]/(df_pi["len1"]+df_pi["len2"])
z = -df_pos["full"]/(df_pi["len1"]+df_pi["len2"])

res = pearsonr(x, y)
print("correlation bppart & pirna: ", res)

res = pearsonr(x, z)
print("correlation bppart & bpmax: ", res)

res = pearsonr(y, z)
print("correlation bpmax & pirna: ", res)

res = spearmanr(x, y)
print("rank correlation bppart & pirna: ", res)

res = spearmanr(x, z)
print("rank correlation bppart & bpmax: ", res)

res = spearmanr(y, z)
print("rank correlation bpmax & pirna: ", res)



from matplotlib.ticker import NullFormatter

nullfmt = NullFormatter()      

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

plt.figure(1, figsize=(8, 8))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

axScatter.scatter(x, y, edgecolors='black')
axScatter.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), 'r')
axScatter.set_ylim(top=0.2)

axScatter.set_xlabel(r'$-ln\:QI$', fontsize=18)
axScatter.set_ylabel(r'$-RT\;ln\:Z$', fontsize=16)

axHistx.hist(x, bins=25, edgecolor='black')
axHistx.set_yticks([2000, 4000, 6000])
axHistx.set_yticklabels(['2k', '4k', '6k'])

axHisty.hist(y, bins=25, orientation='horizontal', edgecolor='black')
axHisty.set_xticks([2000, 4000, 6000])
axHisty.set_xticklabels(['2k', '4k', '6k'])

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())
plt.savefig("cor_Pearson_37")

