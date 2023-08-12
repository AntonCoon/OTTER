import sys
import os
import subprocess as sp

import scripts.cli as cli
import importlib
importlib.reload(cli)

if len(sys.argv) != 3:
    cli.print_usage()
    sys.exit(0)

fltr = sys.argv[1]
assert (fltr == '--filter-soft' or fltr == '--filter-medium' or fltr == '--filter-hard')

folder = sys.argv[2]
assert os.path.isdir(folder)
if folder[-1] != '/':
    folder += '/'
output_folder = folder + 'pipeline-out/'



filtering = True
cor_clust = True

if filtering:
    print('Step 1: Filter vcf...')
    if fltr == '--filter-soft':
        cmd = ['ipython', 'scripts/filter.py', '--', '--soft', folder]
        filtered_vcf = output_folder + 'filtered_soft.vcf'
    elif fltr == '--filter-medium':
        cmd = ['ipython', 'scripts/filter.py', '--', '--medium', folder]
        filtered_vcf = output_folder + 'filtered_medium.vcf'
    else:
        cmd = ['ipython', 'scripts/filter.py', '--', '--hard', folder]
        filtered_vcf = output_folder + 'filtered_hard.vcf'
    sp.run(cmd, check=True)


    print('Step 2: Annotate filtered vcf...')
    cmd = ['ipython', 'scripts/snpeff_variantcalling.py', filtered_vcf, output_folder]
    sp.run(cmd, check=True)

if cor_clust:
    print('Step 3: Correlation clustering...')
    cmd = ['ipython', 'scripts/correlation_clustering.py', folder]
    sp.run(cmd, check=True)
