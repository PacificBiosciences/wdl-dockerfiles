FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG SAWFISH_VERSION
RUN wget https://github.com/PacificBiosciences/sawfish/releases/download/v${SAWFISH_VERSION}/sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu/bin/sawfish /usr/local/bin/ && \
	chmod +x /usr/local/bin/sawfish
