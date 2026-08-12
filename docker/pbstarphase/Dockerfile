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

ARG PBSTARPHASE_VERSION
ARG CPIC_BUILD_VERSION
RUN wget https://github.com/PacificBiosciences/pb-StarPhase/archive/refs/tags/v${PBSTARPHASE_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf v${PBSTARPHASE_VERSION}.tar.gz --directory /opt \
	&& rm v${PBSTARPHASE_VERSION}.tar.gz \
	&& ln -s /opt/pb-StarPhase-${PBSTARPHASE_VERSION}/data/${CPIC_BUILD_VERSION} /opt/pbstarphase_db.json.gz
RUN wget https://github.com/PacificBiosciences/pb-StarPhase/releases/download/v${PBSTARPHASE_VERSION}/pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& tar --no-same-owner -zxvf pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt \
	&& rm pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& mv /opt/pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu/pbstarphase /usr/local/bin \
	&& chmod +x /usr/local/bin/pbstarphase

RUN pbstarphase --version
