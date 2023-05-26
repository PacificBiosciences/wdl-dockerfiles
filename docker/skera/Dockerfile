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

ARG SKERA_VERSION
RUN wget https://github.com/PacificBiosciences/skera/releases/download/v${SKERA_VERSION}/skera --directory-prefix /opt/skera \
	&& chmod +x /opt/skera/skera \
	&& ln -s /opt/skera/skera /usr/local/bin
