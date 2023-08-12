import sys
import subprocess as sp
from pathlib import Path

import importlib
import scripts.settings as settings
importlib.reload(settings)

inp_vcf = sys.argv[1]
output_folder = sys.argv[2]
if output_folder[-1] != '/':
    output_folder += '/'

annotated_vcf = inp_vcf[:-4] + '_snpEff-ann.vcf'

cmd = f"java -Xmx8g -jar bin/snpEff/snpEff.jar GRCh38.105 -onlyProtein -no-intergenic -no-downstream -no-intron -no-upstream {inp_vcf} > {annotated_vcf}"
if settings.VERBOSE:
    print(cmd)
sp.run(cmd, shell=True, check=True)

bn = Path(inp_vcf).stem

sp.run(['mv', 'snpEff_genes.txt', output_folder + bn + '_snpEff_genes.txt'], check=True)
sp.run(['mv', 'snpEff_summary.html', output_folder + bn + '_snpEff_summary.html'], check=True)
