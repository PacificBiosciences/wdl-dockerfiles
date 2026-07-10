FROM quay.io/pacbio/pb_wdl_base@sha256:fa70d211e884a7348b53ff0120eed1f8070f053dd4a1c4b0bfeec74e1a42c5b0

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG IMAGE_DESCRIPTION
ARG IMAGE_VENDOR
ARG IMAGE_URL
ARG IMAGE_MAINTAINER

LABEL org.opencontainers.image.title="${IMAGE_NAME}" \
	org.opencontainers.image.description="${IMAGE_DESCRIPTION}" \
	org.opencontainers.image.authors="${IMAGE_MAINTAINER}" \
	org.opencontainers.image.vendor="${IMAGE_VENDOR}" \
	org.opencontainers.image.url="${IMAGE_URL}" \
	org.opencontainers.image.version="${IMAGE_TAG}"

ARG GFATOOLS_VERSION
RUN wget https://github.com/lh3/gfatools/archive/refs/tags/v${GFATOOLS_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf v${GFATOOLS_VERSION}.tar.gz --directory /opt \
	&& rm v${GFATOOLS_VERSION}.tar.gz
RUN cd /opt/gfatools-${GFATOOLS_VERSION} \
	&& make -j \
	&& ln -s /opt/gfatools-${GFATOOLS_VERSION}/gfatools /usr/local/bin/ \
	&& make clean

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
