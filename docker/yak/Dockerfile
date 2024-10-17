FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG YAK_VERSION
RUN wget https://github.com/lh3/yak/archive/refs/tags/v${YAK_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf v${YAK_VERSION}.tar.gz --directory /opt \
	&& rm v${YAK_VERSION}.tar.gz
RUN cd /opt/yak-${YAK_VERSION} \
	&& make -j \
	&& ln -s /opt/yak-${YAK_VERSION}/yak /usr/local/bin/
