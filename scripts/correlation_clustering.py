import os
import sys
import subprocess as sp

import scripts.settings as settings
import importlib
importlib.reload(settings)


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


second_patient_path = Path(sys.argv[1])
patient = second_patient_path.stem.lower()
results_folder = sys.argv[2]
if results_folder[-1] != '/':
    results_folder += '/'


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



time_point_amount = []



for k, v in allele_frequency_first_time_point.items():
    time_point_amount.append(len({tp for tp, _ in v}))


not_rear_mut = defaultdict(list)
for k, v in allele_frequency_first_time_point.items():

    amount = len({tp for tp, _ in v})
    name = "_".join(map(str, k))
    if amount >= 2:
        not_rear_mut[amount].append(("_".join(map(str, sorted({tp for tp, _ in v}))), name))




with open(f"{settings.TMP_FOLDER}/not_rare_{patient.lower()}.csv", "w") as nr_file:
    for amount in [2, 3, 4, 5]:
        for points, name in not_rear_mut[amount]:
            nr_file.write(f"{amount}, {points}, {name}\n")




all_mutations_are_presented = dict()
for k, v in allele_frequency_first_time_point.items():
    amount = len({tp for tp, _ in v})
    if amount == 5:
        name = "_".join(map(str, k))
        freqs = [[] for _ in range(5)]
        for idx, freq in v:
            freqs[idx - 1].append(freq)
        all_mutations_are_presented[name] = freqs




ordered_names = list(all_mutations_are_presented.keys())




a_freqs = [[] for _ in range(len(ordered_names))]
for idx, name in enumerate(ordered_names):
    a_freqs[idx] = numpy.array(list(map(numpy.mean, all_mutations_are_presented[name])))

if not a_freqs:
    print("Not enough mutations found for correlation clustering analysis...")
    sys.exit(0)



a_freqs = numpy.array(a_freqs)




_min = numpy.min(a_freqs, axis=1)[numpy.newaxis].T
_max = numpy.max(a_freqs, axis=1)[numpy.newaxis].T
a_freqs = (a_freqs - _min) / (_max - _min)




freq_norm_dict = {ordered_names[idx]: freqs for idx, freqs in enumerate(a_freqs)}




d = pandas.DataFrame(data=a_freqs.T, columns=ordered_names)



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

