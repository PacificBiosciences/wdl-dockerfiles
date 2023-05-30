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

ARG PIGEON_VERSION
RUN wget https://github.com/PacificBiosciences/pigeon/releases/download/v${PIGEON_VERSION}/pigeon --directory-prefix /opt/pigeon \
	&& chmod +x /opt/pigeon/pigeon \
	&& ln -s /opt/pigeon/pigeon /usr/local/bin
