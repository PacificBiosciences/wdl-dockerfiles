FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG PBSV_VERSION
RUN wget https://github.com/PacificBiosciences/pbsv/releases/download/v${PBSV_VERSION}/pbsv --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/pbsv
