FROM ubuntu:focal

MAINTAINER William Rowell <wrowell@pacb.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		python3-dev \
		python3-pip \
	&& rm -rf /var/lib/apt/lists/*

ARG JELLYFISH_VERSION
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-linux \
	&& mv jellyfish-linux /usr/local/bin/jellyfish \
	&& chmod +x /usr/local/bin/jellyfish

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}:/opt/scripts/"
