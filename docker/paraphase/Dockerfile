FROM ubuntu:focal

MAINTAINER William Rowell <wrowell@pacificbiosciences.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
        python3 \
        python3-pip \
	&& rm -rf /var/lib/apt/lists/*

ARG MINIMAP2_VERSION
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf minimap2-${MINIMAP2_VERSION}.tar.bz2 --directory /opt \
	&& rm minimap2-${MINIMAP2_VERSION}.tar.bz2
RUN cd /opt/minimap2-${MINIMAP2_VERSION} \
	&& make \
	&& ln -s /opt/minimap2-${MINIMAP2_VERSION}/minimap2 /usr/local/bin/

ARG SAMTOOLS_VERSION
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} \
	&& make \
	&& make install

ARG PARAPHASE_VERSION
RUN python3 -m pip install paraphase==${PARAPHASE_VERSION}
