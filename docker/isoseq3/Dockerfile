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

ARG ISOSEQ3_VERSION
RUN wget https://github.com/PacificBiosciences/IsoSeq/releases/download/v${ISOSEQ3_VERSION}/isoseq3 --directory-prefix /opt/isoseq3 \
	&& chmod +x /opt/isoseq3/isoseq3 \
	&& ln -s /opt/isoseq3/isoseq3 /usr/local/bin
