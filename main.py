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

print('Step 1: Filter vcf...')
if fltr == '--filter-soft':
    sys.exit(1)
elif fltr == '--filter-medium':
    cmd = ['ipython', 'scripts/filter.py', '--', '--medium', folder]
    # sys.exit(1)
else:
    cmd = ['ipython', 'scripts/filter.py', '--', '--hard', folder]

sp.run(cmd, check=True)
