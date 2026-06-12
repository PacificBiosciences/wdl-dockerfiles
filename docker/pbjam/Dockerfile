FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Matt Holt <mholt@pacb.com>, Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt update && apt install -y --no-install-recommends \
		fonts-dejavu \
		&& apt clean \
		&& rm -rf /var/lib/apt/lists/*

ARG PBJAM_VERSION
RUN wget https://github.com/PacificBiosciences/pbjam/releases/download/v${PBJAM_VERSION}/pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu/pbjam /usr/local/bin/
