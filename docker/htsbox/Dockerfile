FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG HTSBOX_VERSION_HASH
RUN wget https://github.com/lh3/htsbox/archive/${HTSBOX_VERSION_HASH}.tar.gz \
	&& tar --no-same-owner -zxvf ${HTSBOX_VERSION_HASH}.tar.gz --directory /opt \
	&& rm ${HTSBOX_VERSION_HASH}.tar.gz
RUN cd /opt/htsbox-${HTSBOX_VERSION_HASH} \
	&& make -j \
	&& ln -s /opt/htsbox-${HTSBOX_VERSION_HASH}/htsbox /usr/local/bin/
