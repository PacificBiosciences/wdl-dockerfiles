FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG HIFIASM_META_VERSION
WORKDIR /opt
RUN wget https://github.com/xfengnefx/hifiasm-meta/archive/refs/tags/hamtv${HIFIASM_META_VERSION}.tar.gz \
	&& tar -zxvf hamtv${HIFIASM_META_VERSION}.tar.gz --directory /opt \
	&& rm hamtv${HIFIASM_META_VERSION}.tar.gz
RUN cd /opt/hifiasm-meta-hamtv${HIFIASM_META_VERSION} \
	&& make -j \
	&& ln -s /opt/hifiasm-meta-hamtv${HIFIASM_META_VERSION}/hifiasm_meta /usr/local/bin/

ARG GFATOOLS_VERSION
RUN wget https://github.com/lh3/gfatools/archive/refs/tags/v${GFATOOLS_VERSION}.tar.gz \
	&& tar -zxvf v${GFATOOLS_VERSION}.tar.gz --directory /opt \
	&& rm v${GFATOOLS_VERSION}.tar.gz
RUN cd /opt/gfatools-${GFATOOLS_VERSION} \
	&& make -j \
	&& ln -s /opt/gfatools-${GFATOOLS_VERSION}/gfatools /usr/local/bin/
