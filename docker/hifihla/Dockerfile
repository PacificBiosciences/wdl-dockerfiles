FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG HIFIHLA_VERSION
RUN wget https://github.com/PacificBiosciences/hifihla/releases/download/v${HIFIHLA_VERSION}/hifihla && \
	mv hifihla /usr/local/bin/hifihla && \
	chmod +x /usr/local/bin/hifihla
