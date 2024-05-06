FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG PBSTARPHASE_VERSION
ARG CPIC_BUILD
RUN wget https://github.com/PacificBiosciences/pb-StarPhase/archive/refs/tags/v${PBSTARPHASE_VERSION}.tar.gz \
	&& tar --no-same-owner -zxvf v${PBSTARPHASE_VERSION}.tar.gz --directory /opt \
	&& rm v${PBSTARPHASE_VERSION}.tar.gz \
	&& ln -s /opt/pb-StarPhase-${PBSTARPHASE_VERSION}/data/${CPIC_BUILD} /opt/pbstarphase_db.json.gz \
	&& wget https://github.com/PacificBiosciences/pb-StarPhase/releases/download/v${PBSTARPHASE_VERSION}/pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& tar --no-same-owner -zxvf pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt \
	&& rm pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& mv /opt/pbstarphase-v${PBSTARPHASE_VERSION}-x86_64-unknown-linux-gnu/pbstarphase /usr/local/bin \
  && chmod +x /usr/local/bin/pbstarphase
