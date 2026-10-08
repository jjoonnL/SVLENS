# SVLENS

Code for *Gene-level analysis of structural variants reveals a lead-variant-centered signal architecture and phenotype-structured recurrence*.

## Project structure

```text
SVLENS/
├── scripts/             # Six analysis stages, in execution order
├── src/acat.py          # Shared ACAT implementation
├── metadata/            # Trait categories, families, and ratio flags
├── examples/            # Required input filenames and columns
├── tests/               # Focused code tests
└── requirements.txt     # Python dependencies
```

## Setup

1. Clone the repository and enter it:

   ```bash
   git clone https://github.com/jjoonnL/SVLENS.git
   cd SVLENS
   ```

2. Create and activate a Conda environment (recommended):

   ```bash
   conda create -n svassoc python=3.14 pip -y
   conda activate svassoc
   ```

3. Install the dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

## Required inputs

Download the source data separately. They are not included in this repository.

| Input | Provider | Expected location |
| --- | --- | --- |
| UK Biobank WGS SV association summary statistics | [deCODE](https://www.decode.com/summarydata/) — UK Biobank Whole-Genome Sequencing Consortium entry | `<DATA_ROOT>/SV_association/` |
| GENCODE v49 basic GTF (GRCh38) | [GENCODE release 49](https://www.gencodegenes.org/human/release_49.html) | `<DATA_ROOT>/Gene_ref/gencode.v49.basic.annotation.gtf.gz` |

Place the NFE quantitative-trait `.txt.gz` files directly in `SV_association/`, retaining their original filenames. The analyzed traits and phenotype groupings are in [`metadata/traits.tsv`](metadata/traits.tsv). See [`examples/input_schema.md`](examples/input_schema.md) for filenames, columns, and coordinate conventions. Keep input data and generated results outside this repository.

## Run the analysis

Run the numbered scripts in order. Step 1 takes explicit `--gtf` and `--out` paths. Steps 2–6 take `--project-root <DATA_ROOT>`. Use `--help` for other options.

| Step | Script | Purpose |
| --- | --- | --- |
| 1 | `01_reference.py` | Build the GENCODE gene table. |
| 2 | `02_sv_acat.py` | Annotate SVs and perform SV-type-stratified gene-level ACAT. |
| 3 | `03_lead_architecture.py` | Identify lead SVs, removal classes, and residual-SV patterns. |
| 4 | `04_recurrence.py` | Analyze exact lead-SV recurrence and phenotype-category enrichment. |
| 5 | `05_sv_features.py` | Compare recurrent-SV features. |
| 6 | `06_maf1_sensitivity.py` | Repeat key SV analyses at MAF <1%. |

From the repository directory with the Conda environment active, replace `/path/to/data` with the directory containing the downloaded inputs. To run all steps:

```bash
(
  set -e
  data_root=/path/to/data
  python scripts/01_reference.py --gtf "$data_root/Gene_ref/gencode.v49.basic.annotation.gtf.gz" --out "$data_root/results/gene_table.parquet"
  python scripts/02_sv_acat.py --project-root "$data_root"
  python scripts/03_lead_architecture.py --project-root "$data_root"
  python scripts/04_recurrence.py --project-root "$data_root"
  python scripts/05_sv_features.py --project-root "$data_root"
  python scripts/06_maf1_sensitivity.py --project-root "$data_root"
)
```

To run one step at a time, set the data directory once and copy each step below in the same terminal:

```bash
data_root=/path/to/data
```

Step 1:

```bash
python scripts/01_reference.py \
  --gtf "$data_root/Gene_ref/gencode.v49.basic.annotation.gtf.gz" \
  --out "$data_root/results/gene_table.parquet"
```

Step 2:

```bash
python scripts/02_sv_acat.py --project-root "$data_root"
```

Step 3:

```bash
python scripts/03_lead_architecture.py --project-root "$data_root"
```

Step 4:

```bash
python scripts/04_recurrence.py --project-root "$data_root"
```

Step 5:

```bash
python scripts/05_sv_features.py --project-root "$data_root"
```

Step 6:

```bash
python scripts/06_maf1_sensitivity.py --project-root "$data_root"
```

## Outputs

All files are written under `<DATA_ROOT>/results/`:

```text
results/
├── gene_table.parquet                 # GENCODE gene reference
├── traits/<TRAIT>/                   # Per-trait SV–gene and ACAT results
├── main/                             # Six principal results
│   ├── gene_trait_associations.csv   # Significant non-ratio associations
│   ├── lead_sv_summary.csv            # Exact lead-SV recurrence and locus flags
│   ├── architecture_summary.csv       # Lead-removal classes and null comparisons
│   ├── category_enrichment.csv        # Trait- and family-level enrichment
│   ├── sv_feature_tests.csv           # Primary and origin-exclusion tests
│   └── maf1_summary.csv               # MAF <1% comparison with the primary analysis
└── work/                             # Files needed between steps and supporting tables
```

Each `traits/<TRAIT>/` folder contains `sv_annotated.parquet`, `sv_weighted.parquet`, and `acat_gene.parquet`. The `work/` folder holds the all-trait association table, lead-removal detail, and supporting feature and MAF <1% tables. The six `main/` files are the starting point for inspecting the results.
