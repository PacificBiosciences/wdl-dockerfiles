FROM python:3.10-buster

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG WHATSHAP_VERSION
RUN python3 -m pip install whatshap==${WHATSHAP_VERSION}
