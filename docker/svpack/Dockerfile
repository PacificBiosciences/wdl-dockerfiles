FROM python:3.9-buster

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG SVPACK_GIT_HASH
WORKDIR /opt
RUN git clone https://github.com/PacificBiosciences/svpack.git \
	&& cd svpack \
	&& git checkout "${SVPACK_GIT_HASH}" \
	&& chmod +x /opt/svpack/svpack

ARG PYSAM_VERSION
# Issue with the use_2to3 library required to build pysam for setuptools > 58
RUN python3 -m pip install setuptools==58
RUN python3 -m pip install pysam==${PYSAM_VERSION}

ENV PATH "${PATH}":/opt/svpack
