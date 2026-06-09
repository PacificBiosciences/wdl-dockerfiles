FROM ubuntu:26.04

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

ARG UV_PYTHON
ENV UV_PYTHON "${UV_PYTHON}"
RUN uv venv /opt/venv --python "${UV_PYTHON}"
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN ln -s $(uv python find "${UV_PYTHON}") /usr/local/bin/python3
RUN ln -s $(uv python find "${UV_PYTHON}") /usr/local/bin/python

RUN apt-get -qq update \
	&& apt-get -qq upgrade \
	&& apt-get -qq install --no-install-recommends \
		curl \
		ca-certificates \
		git \
		wget \
		build-essential \
		zlib1g-dev \
		liblzma-dev \
		libbz2-dev \
		libncurses5-dev \
		libcurl4-openssl-dev \
		libdeflate-dev \
		jq \
		gawk \
		parallel \
	&& rm -rf /var/lib/apt/lists/*

ARG HTSLIB_VERSION
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf htslib-${HTSLIB_VERSION}.tar.bz2 --directory /opt \
	&& rm htslib-${HTSLIB_VERSION}.tar.bz2
RUN cd /opt/htslib-${HTSLIB_VERSION} \
	&& make -j \
	&& make install

ARG BCFTOOLS_VERSION
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN cd /opt/bcftools-${BCFTOOLS_VERSION} \
	&& make -j \
	&& make install

ARG SAMTOOLS_VERSION
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} \
	&& make -j \
	&& make install

ARG BEDTOOLS_VERSION
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static \
	-O /usr/local/bin/bedtools \
	&& chmod +x /usr/local/bin/bedtools

ARG NUMPY_VERSION
ARG PANDAS_VERSION
ARG MATPLOTLIB_VERSION
ARG SEABORN_VERSION
ARG PYSAM_VERSION
ARG VCFPY_VERSION
ARG BIOPYTHON_VERSION
RUN uv pip install \
    "numpy~=${NUMPY_VERSION}" \
    "pandas~=${PANDAS_VERSION}" \
    "matplotlib~=${MATPLOTLIB_VERSION}" \
    "seaborn~=${SEABORN_VERSION}" \
    "pysam~=${PYSAM_VERSION}" \
    "biopython~=${BIOPYTHON_VERSION}"
