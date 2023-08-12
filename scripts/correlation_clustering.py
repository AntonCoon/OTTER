#!/usr/bin/env python
# coding: utf-8

import os
import sys
import subprocess as sp

import scripts.settings as settings
import importlib
importlib.reload(settings)

# In[99]:


import pandas
import vcfpy
from pathlib import Path
from collections import defaultdict
import seaborn
from collections import Counter
import numpy
from matplotlib import pyplot as plt
import scipy
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

seaborn.set_style("whitegrid")

# In[353]:

second_patient_path = Path(sys.argv[1])
patient = second_patient_path.stem.lower()
results_folder = second_patient_path / Path('pipeline-out')
os.makedirs(results_folder, exist_ok=True)


allele_frequency_first_time_point = defaultdict(list)
for ts in range(1, 6):
    for repr in [1, 2]:
        reader = vcfpy.Reader.from_path(second_patient_path / f"ts_{ts}_rep{repr}.vcf")

        for record in reader:
            key = (
                record.CHROM,
                record.POS,
                record.REF,
                record.ALT[0].type,
                record.ALT[0].value,
            )
            allele_frequency_first_time_point[key].append((ts, record.INFO["AF"]))


# In[354]:


time_point_amount = []


# In[355]:


for k, v in allele_frequency_first_time_point.items():
    time_point_amount.append(len({tp for tp, _ in v}))


# In[356]:


# In[357]:


not_rear_mut = defaultdict(list)
for k, v in allele_frequency_first_time_point.items():

    amount = len({tp for tp, _ in v})
    name = "_".join(map(str, k))
    if amount >= 2:
        not_rear_mut[amount].append(("_".join(map(str, sorted({tp for tp, _ in v}))), name))


# In[358]:


with open(f"{settings.TMP_FOLDER}/not_rare_{patient.lower()}.csv", "w") as nr_file:
    for amount in [2, 3, 4, 5]:
        for points, name in not_rear_mut[amount]:
            nr_file.write(f"{amount}, {points}, {name}\n")


# In[359]:


all_mutations_are_presented = dict()
for k, v in allele_frequency_first_time_point.items():
    amount = len({tp for tp, _ in v})
    if amount == 5:
        name = "_".join(map(str, k))
        freqs = [[] for _ in range(5)]
        for idx, freq in v:
            freqs[idx - 1].append(freq)
        all_mutations_are_presented[name] = freqs


# In[360]:


ordered_names = list(all_mutations_are_presented.keys())


# In[361]:


a_freqs = [[] for _ in range(len(ordered_names))]
for idx, name in enumerate(ordered_names):
    a_freqs[idx] = numpy.array(list(map(numpy.mean, all_mutations_are_presented[name])))


# In[362]:


a_freqs = numpy.array(a_freqs)


# In[363]:


_min = numpy.min(a_freqs, axis=1)[numpy.newaxis].T
_max = numpy.max(a_freqs, axis=1)[numpy.newaxis].T
a_freqs = (a_freqs - _min) / (_max - _min)


# In[364]:


freq_norm_dict = {ordered_names[idx]: freqs for idx, freqs in enumerate(a_freqs)}


# In[365]:


d = pandas.DataFrame(data=a_freqs.T, columns=ordered_names)


# In[366]:


corr = d.corr()
mask = numpy.triu(numpy.ones_like(corr, dtype=bool))
f, ax = plt.subplots(figsize=(11, 9))
cmap = seaborn.dark_palette("#FFFFFF", as_cmap=True, reverse=True)
seaborn.heatmap(
    corr,
    mask=mask,
    cmap=cmap,
    center=0,
    square=True,
    linewidths=.5,
    cbar_kws={"shrink": .5},
    vmin=corr.to_numpy().min(),
    vmax=corr.to_numpy().max(),
)
plt.savefig(f"{str(results_folder)}/bh_2_correlation.pdf", bbox_inches = 'tight')
plt.close()


# In[369]:


labels = ordered_names
p = len(labels)

plt.figure(figsize=(20,8))
plt.title('Hierarchical Clustering Dendrogram', fontsize=20)
plt.xlabel('mutations', fontsize=16)
plt.ylabel('distance', fontsize=16)

Z = linkage(corr, 'ward')

R = dendrogram(
    Z,
    p=p,  # show only the last p merged
    no_plot=True,
)

temp = {R["leaves"][i]: labels[R["leaves"][i] - 1] for i in range(len(R["leaves"]))}
def llf(xx):
    lab = "{}".format(temp[xx])
    lab = lab.replace('_', ' ').split()
    lab = ' '.join(lab[i] for i in [0, 1, 2, 4])
    return lab


with plt.rc_context({'lines.linewidth': 4.3}):
    x = dendrogram(
        Z,
        truncate_mode='lastp',  # show only the last p merged
        p=p,  # show only the last p merged clusters
        leaf_label_func=llf,
        leaf_rotation=90.,
        leaf_font_size=12.,
        show_contracted=True,  # to get a distribution impression in truncated branches
        # color_threshold = 6,
    )
    plt.savefig(f"{str(results_folder)}/{patient.lower()}_dendrogram.pdf", bbox_inches = 'tight')
plt.close()

pd = pandas
d = pd.DataFrame([x['ivl'], x['leaves_color_list']]).T
d.columns = ['var', 'cluster']
d['cluster'] = d['cluster'].apply(lambda x: int(x[1:]))

# convert to vcf
temp_vcf = settings.TMP_FOLDER + 'corr_clust.vcf'
output_vcf = str(results_folder) + '/corr_clust.vcf'

vcf_columns = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR".split()

with open(output_vcf, 'w') as fout:
    fout.write('\t'.join(vcf_columns) + '\n')
    for i, row in d.iterrows():
        chrom, pos, ref, alt = row['var'].split()
        cluster = row['cluster']
        fout.write('\t'.join([chrom, pos, '.', ref, alt, '20', 'PASS', f'cluster={cluster}', '.', '.']) + '\n')

# run variant calling...
cmd = ['ipython', 'scripts/snpeff_variantcalling.py', output_vcf, str(results_folder)]
print(' '.join(cmd))
sp.run(cmd, check=True)

# skip following cells
sys.exit()

# In[325]:


[idx for idx, name in enumerate(list(temp.values())) if name == "chr8_32685712_C_SNV_T" or name == "chr21_41473456_A_SNV_G"]


# In[327]:


for name in list(temp.values())[71:74]:
    plt.plot(freq_norm_dict[name])


# In[ ]:





# In[264]:


plt.plot(freq_norm_dict["chr15_87963853_T_SNV_A"])
plt.plot(freq_norm_dict["chr1_11121264_C_SNV_T"])


# In[217]:


groups = [list(temp.values())[:10], list(temp.values())[10:17], list(temp.values())[17:]]


# In[218]:


first_group = {
    "mutation": [],
    "frequency": [],
    "point": [],
}

for name in groups[0]:
    allele_freq_line = list(freq_norm_dict[name])
    first_group["mutation"].extend([name] * len(allele_freq_line))
    first_group["frequency"].extend(allele_freq_line)
    first_group["point"].extend(list(range(1, len(allele_freq_line) + 1)))


# In[219]:


plt.figure(figsize=(12,8))
ax = seaborn.lineplot(
    pandas.DataFrame.from_dict(first_group),
    x="point",
    y="frequency",
    hue="mutation",
    style="mutation",
    markers=["o"] * 10
)
seaborn.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))


# In[220]:


second_group = {
    "mutation": [],
    "frequency": [],
    "point": [],
}

for name in groups[1]:
    allele_freq_line = list(freq_norm_dict[name])
    second_group["mutation"].extend([name] * len(allele_freq_line))
    second_group["frequency"].extend(allele_freq_line)
    second_group["point"].extend(list(range(1, len(allele_freq_line) + 1)))


# In[221]:


plt.figure(figsize=(12,8))
ax = seaborn.lineplot(
    pandas.DataFrame.from_dict(second_group),
    x="point",
    y="frequency",
    hue="mutation",
    style="mutation",
    markers=["o"] * len(groups[1])
)
seaborn.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))


# In[222]:


third_group = {
    "mutation": [],
    "frequency": [],
    "point": [],
}

for name in groups[2]:
    allele_freq_line = list(freq_norm_dict[name])
    third_group["mutation"].extend([name] * len(allele_freq_line))
    third_group["frequency"].extend(allele_freq_line)
    third_group["point"].extend(list(range(1, len(allele_freq_line) + 1)))


# In[223]:


plt.figure(figsize=(12,8))
ax = seaborn.lineplot(
    pandas.DataFrame.from_dict(third_group),
    x="point",
    y="frequency",
    hue="mutation",
    style="mutation",
    markers=["o"] * len(groups[2])
)
seaborn.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))


# In[224]:


plt.figure(figsize=(12,8))
ax = seaborn.lineplot(
    pandas.DataFrame.from_dict(first_group),
    x="point",
    y="frequency",
)
ax = seaborn.lineplot(
    pandas.DataFrame.from_dict(second_group),
    x="point",
    y="frequency",
)
ax = seaborn.lineplot(
    pandas.DataFrame.from_dict(third_group),
    x="point",
    y="frequency",
)


# In[ ]:




