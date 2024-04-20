FROM python:3.11-bookworm

MAINTAINER Billy Rowell <wrowell@pacificbiosciences.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
		&& apt-get -qq install \
		libjemalloc2 \
		procps \
		git \
		build-essential \
		zlib1g-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
		tabix \
		&& rm -rf /var/lib/apt/lists/*

ARG SENTIEON_VERSION
ENV SENTIEON_VERSION "${SENTIEON_VERSION}"
RUN mkdir -p /opt/sentieon/ \
	&& wget "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" \
	&& tar --no-same-owner -zxf sentieon-genomics-${SENTIEON_VERSION}.tar.gz --directory /opt/sentieon/ \
	&& rm sentieon-genomics-${SENTIEON_VERSION}.tar.gz

ARG SENTIEON_PACBIO_MODEL_VERSION
ENV SENTIEON_PACBIO_MODEL_VERSION "${SENTIEON_PACBIO_MODEL_VERSION}"
RUN mkdir -p /opt/sentieon/model \
	&& wget "https://s3.amazonaws.com/sentieon-release/other/${SENTIEON_PACBIO_MODEL_VERSION}.bundle" \
	&& mv "${SENTIEON_PACBIO_MODEL_VERSION}.bundle" /opt/sentieon/model/
ENV MODEL_BUNDLE="/opt/sentieon/model/${SENTIEON_PACBIO_MODEL_VERSION}.bundle"

ARG BCFTOOLS_VERSION
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN cd /opt/bcftools-${BCFTOOLS_VERSION} \
	&& make \
	&& make install

ARG BEDTOOLS_VERSION
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools.static \
	-O /usr/local/bin/bedtools \
	&& chmod +x /usr/local/bin/bedtools

ARG SENTIEON_SCRIPTS_GIT_HASH
ENV SENTIEON_SCRIPTS_GIT_HASH "${SENTIEON_SCRIPTS_GIT_HASH}"
WORKDIR /opt
RUN git clone https://github.com/Sentieon/sentieon-scripts.git \
	&& cd sentieon-scripts \
	&& git checkout "${SENTIEON_SCRIPTS_GIT_HASH}" \
	&& wget https://raw.githubusercontent.com/Sentieon/sentieon-cli/main/sentieon_cli/scripts/gvcf_combine.py \
	&& mv gvcf_combine.py /opt/sentieon-scripts/dnascope_LongRead/ \
	&& chmod +x /opt/sentieon-scripts/dnascope_LongRead/*.py \
	&& chmod +x /opt/sentieon-scripts/dnascope_LongRead/*.sh
ENV PATH ${PATH}:/opt/sentieon-scripts/dnascope_LongRead 

ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-$SENTIEON_VERSION
ENV PATH $SENTIEON_INSTALL_DIR/bin/:$PATH
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000

CMD ["/bin/bash"]