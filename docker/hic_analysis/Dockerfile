# minimal Debian base
FROM debian:trixie-slim

# build-time only (ARG does not persist into the running container)
ARG DEBIAN_FRONTEND=noninteractive
ARG MICROMAMBA_VERSION=2.9.0-0
ARG IMAGE_VERSION=0.2.1

LABEL org.opencontainers.image.title="hic_analysis" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.description="Hi-C analysis image for the GEP2 Snakemake pipeline"

# environment lives in /opt/conda; micromamba's package cache in /opt/micromamba
ENV CONDA_DIR=/opt/conda
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH="$CONDA_DIR/bin:$PATH"

# fail a RUN if any command in a pipe fails, not just the last one
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# 1. System dependencies and locales
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates curl procps locales && \
    sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en

# 2. micromamba as a standalone static binary
RUN curl -fsSL -o /usr/local/bin/micromamba \
        "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" && \
    chmod +x /usr/local/bin/micromamba && \
    micromamba --version

# 3. Conda environment:
#    save fully resolved environment (exact package URLs) inside the image.
COPY hic_analysis_docker_env.yml /opt/hic_analysis_docker_env.yml
RUN micromamba config set channel_priority flexible && \
    micromamba create -y -p "$CONDA_DIR" -f /opt/hic_analysis_docker_env.yml && \
    micromamba env export -p "$CONDA_DIR" --explicit > /opt/hic_analysis_docker_env.lock.txt && \
    micromamba clean -a -y

# 4. Verify: every executable exists and the Python packages import.
#    Reports all missing tools before failing.
RUN missing=0; \
    for t in micromamba python pigz seqtk seqkit chromap pairtools cooler \
             hicPlotMatrix hicConvertFormat yahs juicer agp_to_fasta sdust longdust \
             PretextMap PretextGraph PretextSnapshot bwa-mem2 minibwa minimap2 \
             samtools mosdepth bedtools diamond tidk; do \
        command -v "$t" >/dev/null || { echo "MISSING: $t"; missing=1; }; \
    done; \
    [ "$missing" -eq 0 ] && \
    python -c "import cooler, pairtools, hicexplorer, matplotlib, numpy" && \
    echo "All tools present."

# Informational only: print versions to the build log (never fails the build)
RUN set +o pipefail; \
    for cmd in "chromap --version" "pairtools --version" "cooler --version" \
               "hicPlotMatrix --version" "samtools --version" "minimap2 --version" \
               "bedtools --version" "diamond version" "tidk --version"; do \
        echo "--- $cmd: $($cmd 2>&1 | head -1)"; \
    done

CMD ["/bin/bash"]
