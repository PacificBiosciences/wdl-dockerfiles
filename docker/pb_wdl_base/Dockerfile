FROM ubuntu:26.04

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG IMAGE_DESCRIPTION
ARG IMAGE_VENDOR
ARG IMAGE_URL
ARG IMAGE_MAINTAINER

LABEL org.opencontainers.image.title="${IMAGE_NAME}" \
      org.opencontainers.image.description="${IMAGE_DESCRIPTION}" \
      org.opencontainers.image.authors="${IMAGE_MAINTAINER}" \
      org.opencontainers.image.vendor="${IMAGE_VENDOR}" \
      org.opencontainers.image.url="${IMAGE_URL}" \
      org.opencontainers.image.version="${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq -y upgrade \
	&& apt-get -qq -y install --no-install-recommends \
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
		libssl-dev \
		libgsl-dev \
		jq \
		gawk \
		parallel \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

ENV UV_PYTHON_INSTALL_DIR=/opt
ARG UV_PYTHON_VERSION
ENV UV_PYTHON_VERSION="${UV_PYTHON_VERSION}"
ENV VIRTUAL_ENV="/opt/venv"
ENV PATH="/opt/venv/bin:$PATH"
RUN uv venv "${VIRTUAL_ENV}" --python "${UV_PYTHON_VERSION}" \
	&& ln -s $(uv python find "${UV_PYTHON_VERSION}") /usr/local/bin/python3 \
	&& ln -s $(uv python find "${UV_PYTHON_VERSION}") /usr/local/bin/python \
	&& chmod -R o+rw "${VIRTUAL_ENV}" /opt/ \
	&& touch /.python_history && chmod a+rw /.python_history

ARG HTSLIB_VERSION
RUN wget "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
	&& tar --no-same-owner -jxvf "htslib-${HTSLIB_VERSION}.tar.bz2" --directory /opt \
	&& rm "htslib-${HTSLIB_VERSION}.tar.bz2" \
	&& cd "/opt/htslib-${HTSLIB_VERSION}" \
	&& ./configure \
		--enable-libcurl \
		--enable-gcs \
		--enable-s3 \
		--with-libdeflate \
	&& make -j \
	&& make install \
	&& ldconfig \
	&& make clean

ARG BCFTOOLS_VERSION
RUN wget "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
	&& tar --no-same-owner -jxvf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" --directory /opt \
	&& rm "bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
	&& cd "/opt/bcftools-${BCFTOOLS_VERSION}" \
	&& ./configure --with-htslib=system --enable-libgsl \
	&& make -j \
	&& make install \
	&& make clean
ENV BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

ARG SAMTOOLS_VERSION
RUN wget "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
	&& tar --no-same-owner -jxvf "samtools-${SAMTOOLS_VERSION}.tar.bz2" --directory /opt \
	&& rm "samtools-${SAMTOOLS_VERSION}.tar.bz2" \
	&& cd "/opt/samtools-${SAMTOOLS_VERSION}" \
	&& ./configure --with-htslib=system \
	&& make -j \
	&& make install \
	&& make clean

ARG BEDTOOLS_VERSION
RUN wget "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static" \
		-O /usr/local/bin/bedtools \
	&& chmod +x /usr/local/bin/bedtools

ARG NUMPY_VERSION
ENV NUMPY_VERSION="${NUMPY_VERSION}"
ARG PANDAS_VERSION
ENV PANDAS_VERSION="${PANDAS_VERSION}"
ARG MATPLOTLIB_VERSION
ENV MATPLOTLIB_VERSION="${MATPLOTLIB_VERSION}"
ARG SEABORN_VERSION
ENV SEABORN_VERSION="${SEABORN_VERSION}"
ARG PYSAM_VERSION
ENV PYSAM_VERSION="${PYSAM_VERSION}"
ARG BIOPYTHON_VERSION
ENV BIOPYTHON_VERSION="${BIOPYTHON_VERSION}"
RUN uv pip install \
	"numpy~=${NUMPY_VERSION}" \
	"pandas~=${PANDAS_VERSION}" \
	"matplotlib~=${MATPLOTLIB_VERSION}" \
	"seaborn~=${SEABORN_VERSION}" \
	"pysam~=${PYSAM_VERSION}" \
	"biopython~=${BIOPYTHON_VERSION}"

RUN htsfile --version \
	&& samtools --version \
	&& bcftools --version \
	&& bcftools plugin --list \
	&& bedtools --version
