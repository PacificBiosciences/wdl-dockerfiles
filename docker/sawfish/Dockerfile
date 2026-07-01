FROM quay.io/pacbio/pb_wdl_base@sha256:fa70d211e884a7348b53ff0120eed1f8070f053dd4a1c4b0bfeec74e1a42c5b0

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
