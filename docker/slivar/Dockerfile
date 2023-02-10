FROM python:3.9-buster

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		build-essential \
		zlib1g-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
		tabix \
	&& rm -rf /var/lib/apt/lists/*

ARG SLIVAR_VERSION
RUN wget https://github.com/brentp/slivar/releases/download/v${SLIVAR_VERSION}/slivar --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/slivar
RUN wget https://github.com/brentp/slivar/releases/download/v${SLIVAR_VERSION}/pslivar --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/pslivar

ARG BCFTOOLS_VERSION
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN cd /opt/bcftools-${BCFTOOLS_VERSION} \
	&& make \
	&& make install

ARG VCFPY_VERSION
ARG PYSAM_VERSION
RUN python3 -m pip install vcfpy==${VCFPY_VERSION} pysam==${PYSAM_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts
