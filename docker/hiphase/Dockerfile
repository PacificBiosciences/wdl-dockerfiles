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
		libncurses5-dev \
        liblzma-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
	&& rm -rf /var/lib/apt/lists/*

ARG HIPHASE_VERSION
RUN wget https://github.com/PacificBiosciences/HiPhase/releases/download/v${HIPHASE_VERSION}/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar -zxvf hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu/hiphase /usr/local/bin/

ARG BCFTOOLS_VERSION
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN cd /opt/bcftools-${BCFTOOLS_VERSION} \
	&& make \
	&& make install

ARG SAMTOOLS_VERSION
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} \
	&& make \
	&& make install
