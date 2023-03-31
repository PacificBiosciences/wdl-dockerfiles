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

ARG YAK_VERSION
RUN wget https://github.com/lh3/yak/archive/refs/tags/v${YAK_VERSION}.tar.gz \
	&& tar -zxvf v${YAK_VERSION}.tar.gz --directory /opt \
	&& rm v${YAK_VERSION}.tar.gz
RUN cd /opt/yak-${YAK_VERSION} \
	&& make \
	&& ln -s /opt/yak-${YAK_VERSION}/yak /usr/local/bin/
