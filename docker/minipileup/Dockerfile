FROM pb_wdl_base:build1

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
	&& make \
	&& ln -s /opt/minipileup-${MINIPILEUP_VERSION_HASH}/minipileup /usr/local/bin/
