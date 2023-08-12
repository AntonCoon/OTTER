# appl - A PiPeLine

Installationg:
```bash
git clone git@github.com:AntonCoon/OncoTracker.git
cd OncoTracker && make
```

Usage (subject to a change):
```bash
python main.py --filter-soft|--filter-medium|--filter-hard| <path_to_folder_with_vcf_files>
```
Options:
- `--filter-soft` leaves all the variants which occured in at least 2 files;
- `--filter-medium` leaves only variants which are present in both replicas in at least 2 timepoints;
- `--filter-hard` leaves only variants which are present in all of the subject's files.

Creates a folder with results in `path_to_folder_with_vcf_files`.
