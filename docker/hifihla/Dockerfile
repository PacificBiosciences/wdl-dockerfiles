FROM python:3.9-buster

LABEL maintainer="Gregory Concepcion <gconcepcion@pacificbiosciences.com>"

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG HIFIHLA_VER

ENV IMAGE_NAME="${IMAGE_NAME}" \
    IMAGE_TAG="${IMAGE_TAG}" \
    PATH="${PATH}:/opt/scripts"


RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        wget \
        git && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /opt/bin && \
    cd /opt/bin && \
    wget https://github.com/PacificBiosciences/hifihla/releases/download/${HIFIHLA_VER}/hifihla && \
    chmod +x hifihla

ENV PATH /opt/bin:$PATH
