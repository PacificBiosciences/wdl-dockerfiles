FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG HIFICNV_VERSION
RUN wget https://github.com/PacificBiosciences/HiFiCNV/releases/download/v${HIFICNV_VERSION}/hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/hificnv-v${HIFICNV_VERSION}-x86_64-unknown-linux-gnu/hificnv /usr/local/bin/
