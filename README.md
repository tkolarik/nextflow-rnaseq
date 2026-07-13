# RNA-seq Quantification Pipeline on Nextflow/AWS

Production-style RNA-seq quantification workflow built with Nextflow DSL2, Docker, and an AWS Batch profile.

The pipeline takes paired-end FASTQ files through read trimming/QC, transcript-level quantification, MultiQC reporting, and count-matrix generation. The `test` profile uses small real paired-end fixtures from the nf-core RNA-seq test dataset; the `live` profile runs a public ENA airway RNA-seq sample against GENCODE human transcripts.

```text
samplesheet.csv -> FASTQ stage/download -> fastp -> Salmon index -> Salmon quant -> count matrix
                                              \                         /
                                               -------- MultiQC --------
```

## Why This Repo Exists

This repository is designed to demonstrate practical workflow engineering for genomics roles: reproducible inputs, containerized tools, resumable execution, cloud configuration, and documented outputs. It is intentionally small enough to run locally while preserving the same execution pattern used for larger RNA-seq studies on AWS Batch.

## Pipeline Steps

| Step | Tool | Purpose |
| --- | --- | --- |
| FASTQ staging | `curl`/Nextflow `path` staging | Pull paired FASTQ files from public URLs or stage local paths from the sample sheet |
| Read trimming/QC | `fastp` | Adapter/quality trimming plus per-sample HTML/JSON QC reports |
| Transcriptome indexing | `salmon index` | Build a Salmon transcriptome index from a FASTA URL |
| Quantification | `salmon quant` | Estimate transcript-level abundance from paired FASTQs |
| Count matrix | `pandas` helper | Merge per-sample `quant.sf` outputs into one count matrix |
| QC aggregation | `MultiQC` | Aggregate fastp and Salmon logs into one report |

## Repository Layout

```text
.
├── main.nf                          # Nextflow DSL2 workflow
├── nextflow.config                  # Local, test, live, and AWS Batch profiles
├── conf/awsbatch.config             # AWS Batch executor configuration
├── Dockerfile                       # Reproducible runtime image
├── envs/pipeline.yml                # Conda/micromamba tool specification
├── bin/
│   ├── make_count_matrix.py         # Quantification aggregation helper
│   └── assert_transcript_counts.py  # Smoke-test output assertion
├── tests/main.nf.test               # nf-test pipeline smoke test
├── .github/workflows/ci.yml         # CI: build image, run smoke test + nf-test
└── data/
    ├── samplesheet.csv              # Two-sample nf-core demo sheet
    ├── test_samplesheet.csv         # One-sample nf-core smoke-test sheet
    └── airway_samplesheet.csv       # Live public ENA airway RNA-seq sample
```

## Requirements

- Nextflow `>=23.10.0`
- Docker for local execution
- AWS Batch, S3, ECR, and IAM permissions for cloud execution

## Quick Start

Build the container:

```bash
docker build -t rna-seq-nextflow-aws:0.1.0 .
```

Run the one-sample smoke test:

```bash
nextflow run main.nf -profile test
```

Run the two-sample demo:

```bash
nextflow run main.nf -profile local --samplesheet data/samplesheet.csv --outdir results/demo
```

Run the live public airway RNA-seq sample:

```bash
nextflow run main.nf -profile live
```

The live profile downloads `SRR1039508` from ENA and the GENCODE v44 human transcriptome from EBI. Expect roughly 2.5 GB of FASTQ input plus the reference download.

Resume after an interruption:

```bash
nextflow run main.nf -profile local -resume
```

## Inputs

The sample sheet is a CSV with one row per paired-end sample:

```csv
sample_id,condition,fastq_1,fastq_2
demo_A,treated,https://example.org/sample_A_R1.fastq.gz,https://example.org/sample_A_R2.fastq.gz
demo_B,control,https://example.org/sample_B_R1.fastq.gz,https://example.org/sample_B_R2.fastq.gz
```

Columns:

| Column | Required | Description |
| --- | --- | --- |
| `sample_id` | Yes | Unique sample identifier used in output names |
| `condition` | No | Metadata column carried into `sample_metadata.tsv` |
| `fastq_1` | Yes | URL or local path to read 1 FASTQ/FASTQ.gz |
| `fastq_2` | Yes | URL or local path to read 2 FASTQ/FASTQ.gz |

The default transcriptome is the compact nf-core prokaryotic reference used by the smoke tests. For real analysis, provide a transcriptome FASTA matching the organism and annotation version:

```bash
nextflow run main.nf \
  -profile local \
  --samplesheet data/samplesheet.csv \
  --transcriptome https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz \
  --outdir results/gencode_demo
```

The checked-in live sample sheet contains verified public ENA FASTQ URLs:

```csv
sample_id,condition,fastq_1,fastq_2
SRR1039508,airway,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz
```

## Outputs

```text
results/
├── fastq/                  # Staged raw FASTQs
├── fastp/                  # fastp HTML/JSON reports and trimmed FASTQs
├── reference/              # Salmon index
├── salmon/                 # Per-sample Salmon quantification directories
├── counts/
│   ├── transcript_counts.tsv # Merged NumReads matrix by transcript ID
│   └── sample_metadata.tsv # Sample metadata used for the matrix
├── multiqc/
│   └── multiqc_report.html # Aggregated QC report
└── pipeline_info/          # Nextflow trace, timeline, DAG, execution report
```

`transcript_counts.tsv` uses Salmon `NumReads` values and keeps transcript identifiers in the first column. In a production gene-level workflow, the same pattern can be extended with `tximport` and a transcript-to-gene map.

## Testing

Run the smoke test and output assertions locally:

```bash
make test
```

The CI workflow in `.github/workflows/ci.yml` builds the Docker image, runs `nextflow run main.nf -profile test`, asserts the transcript count matrix shape/content, installs nf-test, and runs the nf-test suite.

## AWS Batch Profile

The `awsbatch` profile is defined in `conf/awsbatch.config`. It expects a prebuilt image in ECR or another registry readable by AWS Batch, plus an S3 work directory.

Build and push the image, replacing account, region, and repository names:

```bash
aws ecr create-repository --repository-name rna-seq-nextflow-aws
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin <account>.dkr.ecr.us-east-1.amazonaws.com
docker build -t rna-seq-nextflow-aws:0.1.0 .
docker tag rna-seq-nextflow-aws:0.1.0 <account>.dkr.ecr.us-east-1.amazonaws.com/rna-seq-nextflow-aws:0.1.0
docker push <account>.dkr.ecr.us-east-1.amazonaws.com/rna-seq-nextflow-aws:0.1.0
```

Run on AWS Batch:

```bash
export AWS_REGION=us-east-1
export AWS_BATCH_QUEUE=nextflow-rnaseq-queue
export NXF_WORK=s3://my-nextflow-work/rnaseq/work
export PIPELINE_OUTDIR=s3://my-nextflow-results/rnaseq/demo
export PIPELINE_CONTAINER=<account>.dkr.ecr.us-east-1.amazonaws.com/rna-seq-nextflow-aws:0.1.0

nextflow run main.nf -profile awsbatch --samplesheet s3://my-inputs/samplesheet.csv
```

AWS resources are intentionally not created by this repo because teams vary in how they provision Batch compute environments, job queues, IAM roles, and VPC networking. The profile is ready to plug into an existing Batch environment.

## Engineering Notes

- Nextflow `publishDir` separates durable results from transient work directories.
- Process labels centralize CPU, memory, and runtime settings for local and cloud execution.
- Retry behavior is configured at the process level in `nextflow.config` instead of being attached to one process ad hoc.
- The container pins bioinformatics tool versions through `envs/pipeline.yml`.
- `-resume` works across local and cloud runs when the same work directory is retained.
- Public ENA and EBI URLs are used in the live profile so reviewers can run the pipeline without requesting protected data.
- The demo builds a Salmon index from the transcriptome FASTA for portability. In a larger production deployment, the reference index should be prebuilt and cached in S3/EFS to avoid repeated GENCODE downloads and indexing.

## Common Commands

Inspect the execution DAG without running all tasks:

```bash
nextflow run main.nf -profile test -preview
```

Clean local runtime artifacts:

```bash
rm -rf work results .nextflow .nextflow.log*
```

## License

MIT. See `LICENSE`.
