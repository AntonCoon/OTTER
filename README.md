[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# OTTER 🦦: Oncology Traces TrackER
A Biohack 2023 project - Analysis of multiple time points of cfDNA from plasma of patients with oncological diagnoses

🦦 consits of filtering, clusterization and variant calling.

## Installation
```bash
git clone git@github.com:AntonCoon/OncoTracker.git
cd OncoTracker && make
```

---

## Requirements
- jre;
- python modules:
    - ipython;
    - pandas;
    - scipy;
    - vcfpy;
    - seaborn.

Alternatively, after installation, a pipeline can be executen within a Docker/Podman container:
```bash
docker build -t otter:1 .
docker run -it --rm --mount type=bind,src="$(pwd)",target=/pipeline otter:1

# then run from the contaier
python3 main.py --filter-medium example_data/BH_2/
```

## Usage

```bash
python main.py --filter-soft|--filter-medium|--filter-hard| <path/to/folder/with/vcf/files>
```
A folder with results be created in `path/to/folder/with/vcf/files`.

Available options are:
- `--filter-soft` leaves all the variants which occur in at least 2 files;
- `--filter-medium` leaves only variants which are present in both replicas in at least 2 timepoints;
- `--filter-hard` leaves only variants which are present in all of the subject's files.

### Results description

Results consist of
- `filtered_LEVEL.vcf` - filtered variants;
- `corr_cust.vcf` - variants filtered and clustered by correlated in variant allele frequency dynamics (here only variants present in all timepoints are considered);
- [SnpEff](https://github.com/pcingola/SnpEff) output:
    - `{filtered_LEVEL,corr_clust}_snpEff-ann.vcf` - annotated variants;
    - `{filtered_LEVEL,corr_clust}_snpEff_genes.txt` - [gene counts summary](https://pcingola.github.io/SnpEff/snpeff/outputsummary/#gene-counts-summary-snpeff_genestxt);
    - `{filtered_LEVEL,corr_clust}_snpEff_summary.html` -  [SnpEff output summary](https://pcingola.github.io/SnpEff/snpeff/outputsummary/#html-summary-snpeff_summaryhtml);
- `*correlation.pdf` - a heatmap of pairwise correlation used for clustering;
- `*dendrogram.pdf` - a dengrogram with variants colored by cluster cluster.
