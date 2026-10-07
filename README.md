# SV association architecture

Code for *Gene-level analysis of structural variants reveals a lead-variant-centered signal architecture and phenotype-structured recurrence*.

## Project structure

```text
sv-association-architecture/
├── scripts/             # Six analysis stages, in execution order
├── src/acat.py          # Shared ACAT implementation
├── metadata/            # Trait categories, families, and ratio flags
├── examples/            # Required input filenames and columns
├── docs/                # Manuscript result-to-code provenance
├── tests/               # Focused code tests
└── requirements.txt     # Python dependencies
```

## Setup

1. Clone the repository and enter it:

   ```bash
   git clone https://github.com/jjoonnL/sv-association-architecture.git
   cd sv-association-architecture
   ```

2. Create and activate a Conda environment (recommended):

   ```bash
   conda create -n sv-architecture python=3.14 pip -y
   conda activate sv-architecture
   ```

3. Install the dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

## Required inputs

Download the source data separately; they are not included in this repository.

| Input | Provider | Expected location |
| --- | --- | --- |
| UK Biobank WGS SV association summary statistics | [deCODE](https://www.decode.com/summarydata/) — UK Biobank Whole-Genome Sequencing Consortium entry | `<DATA_ROOT>/SV_association/NFE/qtbig.A/` |
| GENCODE v49 basic GTF (GRCh38) | [GENCODE release 49](https://www.gencodegenes.org/human/release_49.html) | `<DATA_ROOT>/Gene_ref/gencode.v49.basic.annotation.gtf.gz` |

The analyzed traits and phenotype groupings are in [`metadata/traits.tsv`](metadata/traits.tsv). See [`examples/input_schema.md`](examples/input_schema.md) for filenames, columns, and coordinate conventions. Keep input data and generated results outside this repository.

## Run the analysis

Run the numbered scripts in order. Step 1 takes explicit `--gtf` and `--out` paths; steps 2–6 take `--project-root <DATA_ROOT>`. Use `--help` for other options.

| Step | Script | Purpose |
| --- | --- | --- |
| 1 | `01_reference.py` | Build the GENCODE gene table. |
| 2 | `02_sv_acat.py` | Annotate SVs and perform SV-type-stratified gene-level ACAT. |
| 3 | `03_lead_architecture.py` | Identify lead SVs, removal classes, and residual-SV patterns. |
| 4 | `04_recurrence.py` | Analyze exact lead-SV recurrence and phenotype-category enrichment. |
| 5 | `05_sv_features.py` | Compare recurrent-SV features. |
| 6 | `06_maf1_sensitivity.py` | Repeat key SV analyses at MAF <1%. |

Use `<DATA_ROOT>/results/gene_table.parquet` as the output for step 1. Steps 2–6 write under `<DATA_ROOT>/results/`. The [provenance map](docs/provenance.md) lists principal outputs and their manuscript figure/table connections. Figure layout and Supplementary Table formatting are not part of this code release.

The complete pipeline has not been rerun from this distribution copy. No source summary statistics or generated results are redistributed here.
