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

ARG ISOSEQ_VERSION
RUN wget https://github.com/PacificBiosciences/IsoSeq/releases/download/v${ISOSEQ_VERSION}/isoseq --directory-prefix /opt/isoseq \
	&& chmod +x /opt/isoseq/isoseq \
	&& ln -s /opt/isoseq/isoseq /usr/local/bin
