FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

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
