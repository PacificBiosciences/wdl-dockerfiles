FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

COPY scripts/* /opt/scripts/
RUN chmod +x /opt/scripts/*

ENV PATH="${PATH}:/opt/scripts"