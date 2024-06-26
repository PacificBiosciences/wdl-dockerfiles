FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG TRGT_VERSION
RUN mkdir /opt/trgt
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& tar --no-same-owner -zxvf  trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& rm trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& mv trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /usr/local/bin/ \
	&& chmod +x /usr/local/bin/trgt

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}:/opt/scripts/"
