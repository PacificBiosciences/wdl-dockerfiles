FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG HIPHASE_VERSION
RUN wget https://github.com/PacificBiosciences/HiPhase/releases/download/v${HIPHASE_VERSION}/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu/hiphase /usr/local/bin/
