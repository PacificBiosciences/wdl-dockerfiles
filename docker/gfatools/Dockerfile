FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG GFATOOLS_VERSION
RUN wget https://github.com/lh3/gfatools/archive/refs/tags/v${GFATOOLS_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf v${GFATOOLS_VERSION}.tar.gz --directory /opt \
	&& rm v${GFATOOLS_VERSION}.tar.gz
RUN cd /opt/gfatools-${GFATOOLS_VERSION} \
	&& make -j \
	&& ln -s /opt/gfatools-${GFATOOLS_VERSION}/gfatools /usr/local/bin/

ARG K8_VERSION
RUN wget https://github.com/attractivechaos/k8/releases/download/v${K8_VERSION}/k8-${K8_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf k8-${K8_VERSION}.tar.bz2 --directory /opt \
	&& rm k8-${K8_VERSION}.tar.bz2
RUN cd /opt/k8-${K8_VERSION} \
	&& ln -s /opt/k8-${K8_VERSION}/k8-x86_64-Linux /usr/local/bin/k8

ARG CALN50_GIT_HASH
WORKDIR /opt
RUN git clone https://github.com/lh3/calN50.git
RUN cd calN50 \
	&& git checkout "${CALN50_GIT_HASH}" \
	&& ln -s /opt/calN50/calN50.js /usr/local/bin/
