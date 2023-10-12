FROM ubuntu:focal

MAINTAINER Billy Rowell <wrowell@pacb.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		build-essential \
	&& rm -rf /var/lib/apt/lists/*

ARG BEDTOOLS_VERSION
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static \
	-O /usr/local/bin/bedtools \
	&& chmod +x /usr/local/bin/bedtools

ARG DATAMASH_VERSION
RUN wget https://ftp.gnu.org/gnu/datamash/datamash-${DATAMASH_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf datamash-${DATAMASH_VERSION}.tar.gz --directory /opt \
	&& rm datamash-${DATAMASH_VERSION}.tar.gz
RUN cd /opt/datamash-${DATAMASH_VERSION} \
	&& ./configure \
	&& make \
	&& make check \
	&& make install
