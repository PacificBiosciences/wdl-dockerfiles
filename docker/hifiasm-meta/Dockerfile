FROM ubuntu:focal

MAINTAINER Karen Fang <karen@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		build-essential \
		zlib1g-dev \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG HIFIASM_META_GIT_HASH
ARG HIFIASM_META_VERSION
ARG HIFIASM_VERSION
WORKDIR /opt
RUN git clone https://github.com/xfengnefx/hifiasm-meta.git \
	&& cd hifiasm-meta \
	&& git checkout "${HIFIASM_META_GIT_HASH}" \
	&& make

ENV PATH "${PATH}":/opt/hifiasm-meta