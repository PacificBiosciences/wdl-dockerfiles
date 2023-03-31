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

ARG HIFIASM_VERSION
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${HIFIASM_VERSION}.tar.gz \
	&& tar -zxvf ${HIFIASM_VERSION}.tar.gz --directory /opt \
	&& rm ${HIFIASM_VERSION}.tar.gz
RUN cd /opt/hifiasm-${HIFIASM_VERSION} \
	&& make \
	&& ln -s /opt/hifiasm-${HIFIASM_VERSION}/hifiasm /usr/local/bin/
