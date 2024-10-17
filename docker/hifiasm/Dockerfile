FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

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
	&& make -j \
	&& ln -s /opt/hifiasm-${HIFIASM_VERSION}/hifiasm /usr/local/bin/
