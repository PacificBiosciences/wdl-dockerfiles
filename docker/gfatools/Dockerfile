FROM ubuntu:focal

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		git \
		wget \
		build-essential \
		zlib1g-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
	&& rm -rf /var/lib/apt/lists/*

ARG GFATOOLS_VERSION
RUN wget https://github.com/lh3/gfatools/archive/refs/tags/v${GFATOOLS_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf v${GFATOOLS_VERSION}.tar.gz --directory /opt \
	&& rm v${GFATOOLS_VERSION}.tar.gz
RUN cd /opt/gfatools-${GFATOOLS_VERSION} \
	&& make \
	&& ln -s /opt/gfatools-${GFATOOLS_VERSION}/gfatools /usr/local/bin/

ARG HTSLIB_VERSION
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf htslib-${HTSLIB_VERSION}.tar.bz2 --directory /opt \
	&& rm htslib-${HTSLIB_VERSION}.tar.bz2
RUN cd /opt/htslib-${HTSLIB_VERSION} \
	&& make \
	&& make install

ARG K8_VERSION
RUN wget https://github.com/attractivechaos/k8/releases/download/${K8_VERSION}/k8-${K8_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf k8-${K8_VERSION}.tar.bz2 --directory /opt \
	&& rm k8-${K8_VERSION}.tar.bz2
RUN cd /opt/k8-${K8_VERSION} \
	&& ln -s /opt/k8-${K8_VERSION}/k8-Linux /usr/local/bin/k8

ARG CALN50_GIT_HASH
WORKDIR /opt
RUN git clone https://github.com/lh3/calN50.git
RUN cd calN50 \
	&& git checkout "${CALN50_GIT_HASH}" \
	&& ln -s /opt/calN50/calN50.js /usr/local/bin/
