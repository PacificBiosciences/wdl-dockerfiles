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

# install fonts for plots
RUN apt update && apt install -y --no-install-recommends \
		fonts-dejavu \
		&& apt clean \
		&& rm -rf /var/lib/apt/lists/*

ARG PBJAM_VERSION
RUN wget https://github.com/PacificBiosciences/pbjam/releases/download/v${PBJAM_VERSION}/pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/pbjam-v${PBJAM_VERSION}-x86_64-unknown-linux-gnu/pbjam /usr/local/bin/

RUN pbjam --version
