FROM quay.io/pacbio/pb_wdl_base@sha256:ef224ec6b77c959baf24bb3c06189acfc62300af92e1e74d311e281735749c56

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG LIMA_VERSION
RUN wget https://github.com/PacificBiosciences/barcoding/releases/download/v${LIMA_VERSION}/lima.tar.gz \
	&& tar --no-same-owner -zxvf lima.tar.gz --directory /usr/local/bin \
	&& chmod +x /usr/local/bin/lima
