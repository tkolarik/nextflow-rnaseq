.PHONY: build test run clean

IMAGE ?= rna-seq-nextflow-aws:0.1.0
NEXTFLOW ?= nextflow

build:
	docker build -t $(IMAGE) .

test:
	$(NEXTFLOW) run main.nf -profile test
	python3 bin/assert_transcript_counts.py results/test/counts/transcript_counts.tsv prok_R1

run:
	$(NEXTFLOW) run main.nf -profile local --samplesheet data/samplesheet.csv --outdir results/demo

clean:
	rm -rf work results .nextflow .nextflow.log*
