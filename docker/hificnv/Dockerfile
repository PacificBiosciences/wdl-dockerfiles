FROM ubuntu:focal

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
	&& rm -rf /var/lib/apt/lists/*

ARG HIFICNV_VERSION
RUN wget https://github.com/PacificBiosciences/HiFiCNV/releases/download/v${HIFICNV_VERSION}/hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu/hificnv /usr/local/bin/

ARG BCFTOOLS_VERSION
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN cd /opt/bcftools-${BCFTOOLS_VERSION} \
	&& make \
	&& make install
