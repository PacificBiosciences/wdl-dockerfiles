FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG JELLYFISH_VERSION
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v${JELLYFISH_VERSION}/jellyfish-${JELLYFISH_VERSION}.tar.gz \
  && tar --no-same-owner -zxvf jellyfish-${JELLYFISH_VERSION}.tar.gz --directory /opt \
	&& rm jellyfish-${JELLYFISH_VERSION}.tar.gz

RUN cd /opt/jellyfish-${JELLYFISH_VERSION} \
  && ./configure \
  && make -j \
  && make install

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}:/opt/scripts/"
