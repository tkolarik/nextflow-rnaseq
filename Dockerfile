FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="rna-seq-nextflow-aws"
LABEL org.opencontainers.image.description="Tools for a compact Nextflow RNA-seq quantification pipeline"

COPY envs/pipeline.yml /tmp/pipeline.yml

RUN micromamba install -y -n base -f /tmp/pipeline.yml && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"

USER root
COPY bin/make_count_matrix.py /usr/local/bin/make_count_matrix.py
RUN chmod +x /usr/local/bin/make_count_matrix.py
USER mambauser

WORKDIR /work
