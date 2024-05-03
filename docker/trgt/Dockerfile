FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG TRGT_VERSION
RUN mkdir /opt/trgt
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& gunzip -c trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz > /usr/local/bin/trgt \
	&& chmod +x /usr/local/bin/trgt
