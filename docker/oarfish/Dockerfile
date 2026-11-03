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

RUN apt-get -qq update \
	&& apt-get -qq -y install --no-install-recommends \
		xz-utils \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

ARG OARFISH_VERSION
RUN wget https://github.com/COMBINE-lab/oarfish/releases/download/v${OARFISH_VERSION}/oarfish-x86_64-unknown-linux-gnu.tar.xz \
	&& wget https://github.com/COMBINE-lab/oarfish/releases/download/v${OARFISH_VERSION}/oarfish-x86_64-unknown-linux-gnu.tar.xz.sha256 \
	&& sha256sum -c oarfish-x86_64-unknown-linux-gnu.tar.xz.sha256 \
	&& tar --no-same-owner -xJf oarfish-x86_64-unknown-linux-gnu.tar.xz --directory /opt \
	&& mv /opt/oarfish-x86_64-unknown-linux-gnu /opt/oarfish-v${OARFISH_VERSION} \
	&& ln -s /opt/oarfish-v${OARFISH_VERSION}/oarfish /usr/local/bin/oarfish \
	&& chmod +x /usr/local/bin/oarfish \
	&& rm oarfish-x86_64-unknown-linux-gnu.tar.xz oarfish-x86_64-unknown-linux-gnu.tar.xz.sha256

RUN oarfish --version
