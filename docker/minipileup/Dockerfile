FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG MINIPILEUP_VERSION_HASH
RUN wget https://github.com/lh3/minipileup/archive/${MINIPILEUP_VERSION_HASH}.tar.gz \
	&& tar --no-same-owner -zxvf ${MINIPILEUP_VERSION_HASH}.tar.gz --directory /opt \
	&& rm ${MINIPILEUP_VERSION_HASH}.tar.gz
RUN cd /opt/minipileup-${MINIPILEUP_VERSION_HASH} \
	&& make -j \
	&& ln -s /opt/minipileup-${MINIPILEUP_VERSION_HASH}/minipileup /usr/local/bin/
