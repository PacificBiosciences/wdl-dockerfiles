FROM ubuntu:jammy

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		unzip \
		python2-dev \
		python3-dev \
		python3-pip \
		zlib1g-dev \
		libncurses5-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
		libgtest-dev \
		pkg-config \
		git \
		cmake \
	&& rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/python python /usr/bin/python2 2

# Need bash for curly brace pattern replace for boost version
SHELL ["/bin/bash", "-c"]
ARG BOOST_VERSION
RUN wget https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/boost_${BOOST_VERSION//./_}.tar.gz \
	&& tar -zxvf boost_${BOOST_VERSION//./_}.tar.gz --directory /opt \
	&& rm boost_${BOOST_VERSION//./_}.tar.gz
RUN cd /opt/boost_${BOOST_VERSION//./_} \
	&& ./bootstrap.sh \
	&& ./b2 \
	&& ./b2 headers \
	&& ./b2 install

ARG MESON_VERSION
RUN python3 -m pip install meson==${MESON_VERSION} cram

ARG HTSLIB_VERSION
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& tar -jxvf htslib-${HTSLIB_VERSION}.tar.bz2 --directory /opt \
	&& rm htslib-${HTSLIB_VERSION}.tar.bz2
RUN cd /opt/htslib-${HTSLIB_VERSION} \
	&& make \
	&& make install

ARG SAMTOOLS_VERSION
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} \
	&& make \
	&& make install

ARG NINJA_VERSION
RUN mkdir /opt/ninja \
	&& wget https://github.com/ninja-build/ninja/releases/download/v${NINJA_VERSION}/ninja-linux.zip \
	&& unzip -d /opt/ninja ninja-linux.zip \
	&& rm ninja-linux.zip \
	&& ln -s /opt/ninja/ninja /usr/local/bin/ninja

ARG PBCOPPER_VERSION
RUN wget https://github.com/PacificBiosciences/pbcopper/archive/refs/tags/v${PBCOPPER_VERSION}.tar.gz \
	&& tar -zxvf v${PBCOPPER_VERSION}.tar.gz --directory /opt \
	&& rm v${PBCOPPER_VERSION}.tar.gz
RUN cd /opt/pbcopper-${PBCOPPER_VERSION} \
	&& meson build . \
	&& cd build \
	&& meson install

ARG PBBAM_VERSION
RUN wget https://github.com/PacificBiosciences/pbbam/archive/refs/tags/v${PBBAM_VERSION}.tar.gz \
	&& tar -zxvf v${PBBAM_VERSION}.tar.gz --directory /opt \
	&& rm v${PBBAM_VERSION}.tar.gz
RUN cd /opt/pbbam-${PBBAM_VERSION} \
	&& meson build . \
	&& meson install -C build

ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/lib/x86_64-linux-gnu/
