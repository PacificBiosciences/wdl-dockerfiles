FROM quay.io/pacbio/pb_wdl_base@sha256:ef224ec6b77c959baf24bb3c06189acfc62300af92e1e74d311e281735749c56
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
