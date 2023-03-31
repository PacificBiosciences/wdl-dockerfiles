FROM python:3.9-buster

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
	&& rm -rf /var/lib/apt/lists/*

ARG PYYAML_VERSION
RUN python3 -m pip install pyyaml==${PYYAML_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts
