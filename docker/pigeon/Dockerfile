FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Jocelyne Bruand <jbruand@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG PIGEON_VERSION
RUN wget https://github.com/PacificBiosciences/pigeon/releases/download/v${PIGEON_VERSION}/pigeon \
	-O /usr/local/bin/pigeon \
	&& chmod +x /usr/local/bin/pigeon
