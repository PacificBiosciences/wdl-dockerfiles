FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG OWL_VERSION
RUN wget https://github.com/PacificBiosciences/owl/releases/download/v${OWL_VERSION}/owl-${OWL_VERSION}-linux-x86_64.tar.gz \
	&& tar --no-same-owner -zxvf owl-${OWL_VERSION}-linux-x86_64.tar.gz --directory /opt \
	&& rm owl-${OWL_VERSION}-linux-x86_64.tar.gz \
	&& mv /opt/owl /usr/local/bin \
	&& chmod +x /usr/local/bin/owl