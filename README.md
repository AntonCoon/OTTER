# appl - A PiPeLine

Installation:
```bash
git clone git@github.com:AntonCoon/OncoTracker.git
cd OncoTracker && make
```

Usage:
``
python main.py --filter-soft|--filter-medium|--filter-hard| <path_to_folder_with_vcf_files>
```
...created a folder with results in `path_to_folder_with_vcf_files`.

Available options are:
- `--filter-soft` leaves all the variants which occur in at least 2 files;
- `--filter-medium` leaves only variants which are present in both replicas in at least 2 timepoints;
- `--filter-hard` leaves only variants which are present in all of the subject's files.


---

Requirements:
- jre
- python modules:
    - ipython
    - pandas
    - scipy
    - vcfpy
    - seaborn

Alternatively, can be deployed within a Docker/Podman container:
```bash
docker build -t otter:1 -f .
docker run -it --rm --mount type=bind,src="$(pwd)",target=/pipeline otter:1

# then run from the contaier
python3 main.py --filter-medium example_data/BH_2/
```
