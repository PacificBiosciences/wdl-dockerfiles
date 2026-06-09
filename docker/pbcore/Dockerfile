FROM quay.io/pacbio/pb_wdl_base@sha256:7e7822f6b540f911d46e2e40dec321a63bed9cae6e51846a4b9c3cae5e682e16
LABEL maintainer="Jocelyne Bruand <jbruand@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME="${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG="${IMAGE_TAG}"

ARG PBCORE_VERSION
ARG SETUPTOOLS_VERSION

ENV DEBIAN_FRONTEND=noninteractive

RUN uv pip install \
    setuptools==${SETUPTOOLS_VERSION} \
    git+https://github.com/PacificBiosciences/pbcore.git@${PBCORE_VERSION}
