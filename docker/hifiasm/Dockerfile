FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG HIFIASM_VERSION
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${HIFIASM_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf ${HIFIASM_VERSION}.tar.gz --directory /opt \
	&& rm ${HIFIASM_VERSION}.tar.gz
RUN cd /opt/hifiasm-${HIFIASM_VERSION} \
	&& make \
	&& ln -s /opt/hifiasm-${HIFIASM_VERSION}/hifiasm /usr/local/bin/
