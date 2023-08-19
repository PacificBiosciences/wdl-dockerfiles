FROM python:3.9-buster

MAINTAINER Billy Rowell <wrowell@pacificbiosciences.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
	&& rm -rf /var/lib/apt/lists/*

ARG DATAMASH_VERSION
RUN wget https://ftp.gnu.org/gnu/datamash/datamash-${DATAMASH_VERSION}.tar.gz \
	&& tar -zxvf datamash-${DATAMASH_VERSION}.tar.gz --directory /opt \
	&& rm datamash-${DATAMASH_VERSION}.tar.gz
RUN cd /opt/datamash-${DATAMASH_VERSION} \
	&& ./configure \
	&& make \
	&& make check \
	&& make install

ARG PYSAM_VERSION
# Issue with the use_2to3 library required to build pysam for setuptools > 58
RUN python3 -m pip install setuptools==58
RUN python3 -m pip install pysam==${PYSAM_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts
