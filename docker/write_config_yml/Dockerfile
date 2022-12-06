FROM ubuntu:focal

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		python3-dev \
		python3-pip \
	&& rm -rf /var/lib/apt/lists/*

ENV PYYAML_VERSION 6.0
RUN python3 -m pip install pyyaml==${PYYAML_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}:/opt/scripts/"
