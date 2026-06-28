FROM quay.io/pacbio/pb_wdl_base@sha256:1fa8b22abd2e32afe18cb2ff41dd0426d813a0450515efc0e626cf35f5658af3

LABEL maintainer="Jocelyne Bruand <jbruand@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG SKERA_VERSION
RUN wget https://github.com/PacificBiosciences/skera/releases/download/v${SKERA_VERSION}/skera -O /usr/local/bin/skera \
	&& chmod +x /usr/local/bin/skera
