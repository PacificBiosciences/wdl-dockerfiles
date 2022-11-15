FROM ubuntu:xenial

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -y update && apt-get -y install \
		build-essential \
		wget \
		tar \
		autoconf \
		gzip \
		bzip2 \
		lbzip2 \
		libbz2-dev \
		libcurl4-openssl-dev \
		liblzma-dev \
		ncurses-dev \
		libncurses5-dev \
		zlib1g-dev \
		libbz2-dev \
		jq \
		curl \
	&& rm -rf /var/lib/apt/lists/*

ARG SAMTOOLS_VERSION
ENV SAMTOOLS_VERSION "${SAMTOOLS_VERSION}"

RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2

WORKDIR /opt/samtools-${SAMTOOLS_VERSION}
RUN autoheader\
	&& autoconf \
	&& ./configure --prefix=/usr/local \
	&& make \
	&& make install

ARG BCFTOOLS_VERSION
ENV BCFTOOLS_VERSION "${BCFTOOLS_VERSION}"

RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar -xvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2

WORKDIR /opt/bcftools-${BCFTOOLS_VERSION}
RUN autoheader \
	&& autoconf \
	&& ./configure --prefix=/usr/local \
	&& make \
	&& make install

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ENV TRGT_VERSION "${IMAGE_TAG}"
ENV TRVZ_VERSION "${IMAGE_TAG}"

RUN mkdir /opt/trgt
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/${TRGT_VERSION}/trgt-${TRGT_VERSION}-linux_x86_64.gz \
	&& gunzip -c trgt-${TRGT_VERSION}-linux_x86_64.gz > /opt/trgt/trgt \
	&& chmod +x /opt/trgt/trgt \
	&& rm trgt-${TRGT_VERSION}-linux_x86_64.gz
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/${TRVZ_VERSION}/trvz-${TRVZ_VERSION}-linux_x86_64.gz \
	&& gunzip -c trvz-${TRVZ_VERSION}-linux_x86_64.gz > /opt/trgt/trvz \
	&& chmod +x /opt/trgt/trvz \
	&& rm trvz-${TRVZ_VERSION}-linux_x86_64.gz
ENV PATH $PATH:/opt/trgt
