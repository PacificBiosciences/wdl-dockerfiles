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
		python3-dev \
		python3-pip \
	&& rm -rf /var/lib/apt/lists/*

ARG TRGT_VERSION
RUN mkdir /opt/trgt
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-linux_x86_64.gz \
	&& gunzip -c trgt-v${TRGT_VERSION}-linux_x86_64.gz > /opt/trgt/trgt \
	&& chmod +x /opt/trgt/trgt
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trvz-v${TRGT_VERSION}-linux_x86_64.gz \
	&& gunzip -c trvz-v${TRGT_VERSION}-linux_x86_64.gz > /opt/trgt/trvz \
	&& chmod +x /opt/trgt/trvz
ENV PATH "${PATH}":/opt/trgt

ARG SAMTOOLS_VERSION
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} \
	&& make \
	&& make install

ARG BCFTOOLS_VERSION
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN cd /opt/bcftools-${BCFTOOLS_VERSION} \
	&& make \
	&& make install

ARG PYSAM_VERSION
RUN python3 -m pip install pysam==${PYSAM_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts
