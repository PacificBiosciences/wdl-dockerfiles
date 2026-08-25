FROM quay.io/pacbio/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae

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

ARG SLIVAR_VERSION
RUN wget https://github.com/brentp/slivar/releases/download/v${SLIVAR_VERSION}/slivar --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/slivar
RUN wget https://github.com/brentp/slivar/releases/download/v${SLIVAR_VERSION}/pslivar --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/pslivar
COPY scripts/* /opt/scripts/
RUN chmod +x /opt/scripts/*
RUN wget https://raw.githubusercontent.com/brentp/slivar/v${SLIVAR_VERSION}/js/slivar-functions.js --directory-prefix /opt/scripts/
ENV PATH="${PATH}:/opt/scripts"

RUN slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+'
