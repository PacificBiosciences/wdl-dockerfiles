FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG MINIMAP2_VERSION
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 \
	&& tar --no-same-owner -jxvf minimap2-${MINIMAP2_VERSION}.tar.bz2 --directory /opt \
	&& rm minimap2-${MINIMAP2_VERSION}.tar.bz2
RUN cd /opt/minimap2-${MINIMAP2_VERSION} \
	&& make -j \
	&& ln -s /opt/minimap2-${MINIMAP2_VERSION}/minimap2 /usr/local/bin/
