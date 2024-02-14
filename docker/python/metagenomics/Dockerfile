FROM python:3.9-buster

LABEL maintainer="Gregory Concepcion <gconcepcion@pacificbiosciences.com>"

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG TOOL_REPO_HASH
ARG NUMPY_VERSION
ARG PANDAS_VERSION
ARG SEABORN_VERSION
ARG BIOPYTHON_VERSION

ENV IMAGE_NAME="${IMAGE_NAME}" \
    IMAGE_TAG="${IMAGE_TAG}" \
    NUMPY_VERSION="${NUMPY_VERSION}" \
    PANDAS_VERSION="${PANDAS_VERSION}" \
    SEABORN_VERSION="${SEABORN_VERSION}" \
    BIOPYTHON_VERSION="${BIOPYTHON_VERSION}" \
    PATH="${PATH}:/opt/scripts"


RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        wget \
        git && \
    rm -rf /var/lib/apt/lists/* && \
    git clone --recurse-submodules https://github.com/PacificBiosciences/pb-metagenomics-tools.git /opt/pb-metagenomics-tools && \
    cd /opt/pb-metagenomics-tools && \
    git checkout "${TOOL_REPO_HASH}" && \
    cp -r /opt/pb-metagenomics-tools/HiFi-MAG-Pipeline/scripts /opt/scripts && \
	sed -i 's|filepath = os.path.join(bin_dir, "{}.fa".format(bin))|filepath = os.path.join("bins", "{}.fa".format(bin))|' /opt/scripts/Filter-Checkm2-Bins.py && \
    chmod +x /opt/scripts/* && \ 
	pip install numpy=="${NUMPY_VERSION}" pandas=="${PANDAS_VERSION}" seaborn=="${SEABORN_VERSION}" biopython=="${BIOPYTHON_VERSION}"

COPY scripts/* /opt/scripts/
RUN chmod +x /opt/scripts/Fasta-Make-Long-Seq-Bins.py
