FROM python:3.9-buster

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG PYYAML_VERSION
ARG NUMPY_VERSION
ARG SCIPY_VERSION
ARG SCIKIT_LEARN_VERSION
ARG PANDAS_VERSION
ARG PYSAM_VERSION
ARG MAPPY_VERSION
# Issue with the use_2to3 library required to build pysam for setuptools > 58
# RUN python3 -m pip install setuptools==58
RUN python3 -m pip install pyyaml==${PYYAML_VERSION} \
	numpy==${NUMPY_VERSION} \
	scipy==${SCIPY_VERSION} \
	scikit-learn==${SCIKIT_LEARN_VERSION} \
	pandas==${PANDAS_VERSION} \
	pysam==${PYSAM_VERSION} \
	mappy==${MAPPY_VERSION}

ARG PANGU_VERSION
WORKDIR /opt
RUN git clone https://github.com/PacificBiosciences/pangu.git \
	&& cd pangu \
	&& git checkout "v${PANGU_VERSION}" \
	&& python3 -m pip install .
