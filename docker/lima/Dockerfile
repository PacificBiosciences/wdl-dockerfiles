FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG LIMA_VERSION
RUN wget https://github.com/PacificBiosciences/barcoding/releases/download/v${LIMA_VERSION}/lima.tar.gz \
	&& tar --no-same-owner -zxvf lima.tar.gz --directory /usr/local/bin \
	&& chmod +x /usr/local/bin/lima
