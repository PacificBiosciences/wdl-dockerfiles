FROM pb_wdl_base:build1

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
	&& make \
	&& ln -s /opt/yak-${YAK_VERSION}/yak /usr/local/bin/
