FROM pb_wdl_base:build1

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
	&& make \
	&& ln -s /opt/htsbox-${HTSBOX_VERSION_HASH}/htsbox /usr/local/bin/
