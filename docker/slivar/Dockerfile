FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG SLIVAR_VERSION
RUN wget https://github.com/brentp/slivar/releases/download/v${SLIVAR_VERSION}/slivar --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/slivar
RUN wget https://github.com/brentp/slivar/releases/download/v${SLIVAR_VERSION}/pslivar --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/pslivar

COPY scripts/* /opt/scripts/
RUN chmod +x /opt/scripts/*

RUN wget https://raw.githubusercontent.com/brentp/slivar/v${SLIVAR_VERSION}/js/slivar-functions.js --directory-prefix /opt/scripts/

ENV PATH "${PATH}:/opt/scripts"
