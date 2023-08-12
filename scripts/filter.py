import sys
import subprocess as sp
import os

import scripts.settings as settings
import importlib
importlib.reload(settings)

verbose = True
LEVEL = sys.argv[1]
assert LEVEL in ('--soft', '--hard', '--medium')
LEVEL = LEVEL[2:]

patient_folder = sys.argv[2]
if patient_folder[-1] != '/':
    patient_folder += '/'

# remove junk from previous runs...
sp.run(f'rm -rf {settings.TMP_FOLDER}*', shell=True, check=True)

# check how many timestamps we have
timestamps = []
for fn in os.listdir(patient_folder):
    if not fn.endswith('.vcf'):
        continue
    ts = int(fn.split('_')[1])
    timestamps.append(ts)
timestamps = list(set(timestamps))
timestamps.sort()

if LEVEL in ('hard', 'medium'):
    # select variants present in both replicates...
    present_in_both_rep_paths = []
    for ts in timestamps:
        rep1 = patient_folder + f"ts_{ts}_rep1.vcf"
        rep2 = patient_folder + f"ts_{ts}_rep2.vcf"

        both_rep_pres_path = settings.TMP_FOLDER + f'both_rep_ts{ts}.tsv'

        cmd = f'scripts/print_present_in_two_vcfs {rep1} {rep2} > {both_rep_pres_path}'
        if verbose:
            print(cmd)
        sp.run(cmd, shell=True, check=True)
        present_in_both_rep_paths.append(both_rep_pres_path)

    if LEVEL == 'hard':
        # from previously selected variants,
        # select only variants that are present in all time points
        cur_tmp_file = present_in_both_rep_paths[0]
        for i in range(1, len(present_in_both_rep_paths)):
            new_tmp_file = str(settings.TMP_FOLDER) + f'/merged_{i+1}.tsv'

            cmd = f'comm -1 -2 {cur_tmp_file} {present_in_both_rep_paths[i]} > {new_tmp_file}'
            if verbose:
                print(cmd)
            sp.run(cmd, shell=True, check=True)
            cur_tmp_file = new_tmp_file
    elif LEVEL == 'medium':
        # from previously selected variants,
        # select only those present in at least 2 time points...
        cur_tmp_file = settings.TMP_FOLDER + 'variants_in_at_least_2.tsv'
        open(cur_tmp_file, 'w').close()

        # add all pairwise intersections, then do sort | uniq
        for i in range(len(present_in_both_rep_paths)):
            for j in range(i + 1, len(present_in_both_rep_paths)):
                sp.run(f'comm -1 -2 {present_in_both_rep_paths[i]} {present_in_both_rep_paths[j]} >> {cur_tmp_file}', shell=True, check=True)

        present_in_at_least_2_path = settings.TMP_FOLDER + 'variants_in_at_least_2_sort_uniq.tsv'
        sp.run(f'sort {cur_tmp_file} | uniq > {present_in_at_least_2_path}', shell=True, check=True)
        cur_tmp_file = present_in_at_least_2_path

elif LEVEL == 'soft':
    cur_tmp_file = settings.TMP_FOLDER + 'variants_in_at_least_2.tsv'
    open(cur_tmp_file, 'w').close()

    filenames = []
    for ts in timestamps:
        for rep in (1, 2):
            filenames.append(patient_folder + f"ts_{ts}_rep{rep}.vcf")
    for i in range(len(filenames)):
        for j in range(i + 1, len(filenames)):
            f1 = filenames[i]
            f2 = filenames[j]
            sp.run(f'scripts/print_present_in_two_vcfs {f1} {f2} >> {cur_tmp_file}', shell=True, check=True)

    present_in_at_least_2_path = settings.TMP_FOLDER + 'variants_in_at_least_2_sort_uniq.tsv'
    sp.run(f'sort {cur_tmp_file} | uniq > {present_in_at_least_2_path}', shell=True, check=True)
    cur_tmp_file = present_in_at_least_2_path

# convert result to vcf
vcf_output = str(settings.TMP_FOLDER) + f'filtered_{LEVEL}.vcf'
vcf_columns = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR".split()
with open(cur_tmp_file, 'r') as fin, open(vcf_output, 'w') as fout:
    fout.write('\t'.join(vcf_columns) + '\n')
    for line in fin:
        line = line.split()
        vcf_fields = [line[0], line[1], '.', line[2], line[3], 20, 'PASS', '.', '.', '.']
        vcf_fields = [str(x) for x in vcf_fields]
        fout.write('\t'.join(vcf_fields) + '\n')

output_folder = patient_folder + 'pipeline-out/'
os.makedirs(output_folder, exist_ok=True)

# sort vcf and place it in output folder
sorted_vcf = output_folder + f'filtered_{LEVEL}.vcf'
if verbose:
    print('Sort filtered vcf...')
sp.run(f'head -n 1 {vcf_output} > {sorted_vcf}', check=True, shell=True)
sp.run(f'tail -n+2 {vcf_output} | sort --key=1,1 --key=2n >> {sorted_vcf}', check=True, shell=True)
