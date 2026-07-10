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

ARG TRGT_VERSION
RUN mkdir /opt/trgt
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& tar --no-same-owner -zxvf  trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& rm trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& mv trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu/trgt /usr/local/bin/ \
	&& chmod +x /usr/local/bin/trgt

COPY scripts/* /opt/scripts/
ENV PATH="${PATH}:/opt/scripts/"

RUN trgt --version
