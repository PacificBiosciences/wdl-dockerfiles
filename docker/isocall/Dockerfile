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

ARG ISOCALL_VERSION
RUN wget https://github.com/PacificBiosciences/isocall/releases/download/${ISOCALL_VERSION}/isocall-v${ISOCALL_VERSION}-linux-x86_64.gz \
	&& gunzip -c isocall-v${ISOCALL_VERSION}-linux-x86_64.gz > /usr/local/bin/isocall \
	&& chmod +x /usr/local/bin/isocall

RUN isocall --version
