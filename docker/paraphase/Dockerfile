FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

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

ARG PARAPHASE_VERSION
RUN python3 -m pip install paraphase==${PARAPHASE_VERSION}
