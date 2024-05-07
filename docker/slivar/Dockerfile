FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

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
