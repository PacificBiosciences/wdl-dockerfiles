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
	&& rm -rf /var/lib/apt/lists/*

ARG HTSBOX_VERSION_HASH
RUN wget https://github.com/lh3/htsbox/archive/${HTSBOX_VERSION_HASH}.tar.gz \
	&& tar --no-same-owner -zxvf ${HTSBOX_VERSION_HASH}.tar.gz --directory /opt \
	&& rm ${HTSBOX_VERSION_HASH}.tar.gz
RUN cd /opt/htsbox-${HTSBOX_VERSION_HASH} \
	&& make \
	&& ln -s /opt/htsbox-${HTSBOX_VERSION_HASH}/htsbox /usr/local/bin/
