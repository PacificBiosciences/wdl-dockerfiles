FROM ubuntu:focal

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
	&& rm -rf /var/lib/apt/lists/*

# Need to explicitly grant executable to non-root users
ARG LIMA_VERSION
RUN wget https://github.com/PacificBiosciences/barcoding/releases/download/v${LIMA_VERSION}/lima.tar.gz \
	&& mkdir /opt/lima \
	&& tar -zxvf lima.tar.gz --directory /opt/lima \
	&& rm lima.tar.gz \
	&& chmod +x /opt/lima/*

ENV PATH ${PATH}:/opt/lima
