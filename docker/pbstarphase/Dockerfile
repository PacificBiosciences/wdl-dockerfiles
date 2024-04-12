FROM python:3.9-buster

LABEL maintainer="Gregory Concepcion <gconcepcion@pacificbiosciences.com>"

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG TOOL_VER

ENV IMAGE_NAME="${IMAGE_NAME}" \
    IMAGE_TAG="${IMAGE_TAG}" \
    TOOL_VER="${TOOL_VER}" \
    PATH="${PATH}:/opt/scripts"


RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        tar \
        wget \
        git && \
    rm -rf /var/lib/apt/lists/* && \
    cd /opt &&\
    ## Unpack the DB
    wget https://github.com/PacificBiosciences/pb-StarPhase/archive/refs/tags/${TOOL_VER}.tar.gz && \
    tar zxvf ${TOOL_VER}.tar.gz && \
    ln -s pbStarPhase-0.9.0/data/v0.9.0/cpic_20240404.json.gz ./pbstarphase_db.json.gz && \
    mkdir /opt/bin && \
    cd /opt/bin && \
    ## Install binary
    wget https://github.com/PacificBiosciences/pb-StarPhase/releases/download/${TOOL_VER}/pbstarphase-${TOOL_VER}-x86_64-unknown-linux-gnu.tar.gz && \
    tar zxvf pbstarphase-${TOOL_VER}-x86_64-unknown-linux-gnu.tar.gz && \
    chmod +x pbstarphase-${TOOL_VER}-x86_64-unknown-linux-gnu/pbstarphase && \
    ln -s /opt/bin/pbstarphase-${TOOL_VER}-x86_64-unknown-linux-gnu/pbstarphase /opt/bin
ENV PATH /opt/bin:$PATH
