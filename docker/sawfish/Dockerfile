FROM quay.io/pacbio/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG IMAGE_DESCRIPTION
ARG IMAGE_VENDOR
ARG IMAGE_URL
ARG IMAGE_MAINTAINER

LABEL org.opencontainers.image.title="${IMAGE_NAME}" \
	org.opencontainers.image.description="${IMAGE_DESCRIPTION}" \
	org.opencontainers.image.authors="${IMAGE_MAINTAINER}" \
	org.opencontainers.image.vendor="${IMAGE_VENDOR}" \
	org.opencontainers.image.url="${IMAGE_URL}" \
	org.opencontainers.image.version="${IMAGE_TAG}"

ARG SAWFISH_VERSION
RUN wget https://github.com/PacificBiosciences/sawfish/releases/download/v${SAWFISH_VERSION}/sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/sawfish-v${SAWFISH_VERSION}-x86_64-unknown-linux-gnu/bin/sawfish /usr/local/bin/ && \
	chmod +x /usr/local/bin/sawfish

ARG SAWSHARK_VERSION
RUN wget https://github.com/PacificBiosciences/sawshark/releases/download/v${SAWSHARK_VERSION}/sawshark-v${SAWSHARK_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf sawshark-v${SAWSHARK_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm sawshark-v${SAWSHARK_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/sawshark-v${SAWSHARK_VERSION}-x86_64-unknown-linux-gnu/bin/sawshark /usr/local/bin/ && \
	chmod +x /usr/local/bin/sawshark

RUN sawfish --version \
	&& sawshark --version
