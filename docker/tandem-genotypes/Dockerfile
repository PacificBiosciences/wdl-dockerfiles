FROM python:3.9-buster

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		gfortran \
		libpcre2-dev \
		openjdk-11-jre \
	&& rm -rf /var/lib/apt/lists/*

# Needed for visualization with tandem genotypes
ARG R_VERSION
RUN R_VERSION_BASE=$(echo "${R_VERSION}" | cut -d '.' -f 1); wget https://cran.r-project.org/src/base/R-${R_VERSION_BASE}/R-${R_VERSION}.tar.gz \
	&& tar -zxvf R-${R_VERSION}.tar.gz --directory /opt \
	&& rm R-${R_VERSION}.tar.gz
RUN cd /opt/R-${R_VERSION} \
	&& ./configure \
	&& make \
	&& make install

ARG TANDEM_GENOTYPES_VERSION
RUN wget https://github.com/mcfrith/tandem-genotypes/archive/refs/tags/${TANDEM_GENOTYPES_VERSION}.tar.gz \
	&& tar -zxvf ${TANDEM_GENOTYPES_VERSION}.tar.gz --directory /opt \
	&& rm ${TANDEM_GENOTYPES_VERSION}.tar.gz
RUN cd /opt/tandem-genotypes-${TANDEM_GENOTYPES_VERSION} \
	&& python3 -m pip install .

ARG PYSAM_VERSION
# Issue with the use_2to3 library required to build pysam for setuptools > 58
RUN python3 -m pip install setuptools==58
RUN python3 -m pip install pysam==${PYSAM_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts
