FROM python:3.9-buster

LABEL maintainer="Gregory Concepcion <gconcepcion@pacificbiosciences.com>"

ARG IMAGE_NAME
ARG IMAGE_TAG

ENV IMAGE_NAME="${IMAGE_NAME}" \
    IMAGE_TAG="${IMAGE_TAG}" \
    PATH="${PATH}:/opt/scripts"


RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        tar \
        wget \
        git && \
    rm -rf /var/lib/apt/lists/* && \
    cd /opt &&\
    wget https://github.com/PacificBiosciences/pb-StarPhase/archive/refs/tags/v0.8.1.tar.gz && \
    tar zxvf v0.8.1.tar.gz && \
    ln -s /opt/pb-StarPhase-0.8.1 /opt/pbstarphase && \
    mkdir /opt/bin && \
    cd /opt/bin && \
    wget https://github.com/PacificBiosciences/pb-StarPhase/releases/download/v0.8.1/pbstarphase-v0.8.1-x86_64-unknown-linux-gnu.tar.gz && \
    tar zxvf pbstarphase-v0.8.1-x86_64-unknown-linux-gnu.tar.gz && \
    chmod +x pbstarphase-v0.8.1-x86_64-unknown-linux-gnu/pbstarphase && \
    ln -s /opt/bin/pbstarphase-v0.8.1-x86_64-unknown-linux-gnu/pbstarphase /opt/bin
ENV PATH /opt/bin:$PATH
