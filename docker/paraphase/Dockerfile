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

ARG PARAPHASE_VERSION
RUN wget https://github.com/PacificBiosciences/paraphase/releases/download/v${PARAPHASE_VERSION}/paraphase-v${PARAPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& tar --no-same-owner -zxvf paraphase-v${PARAPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt \
	&& rm paraphase-v${PARAPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& ln -s /opt/paraphase-v${PARAPHASE_VERSION}-x86_64-unknown-linux-gnu/paraphase /usr/local/bin/

RUN paraphase --version
